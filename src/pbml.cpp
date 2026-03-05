/*
===============
Intro
===============
This is the implementation of PBML (Positional Boyer-Moore-Li): A method to enumerate SMEMs (Set Maximal Exact Matches) from external 
query haplotypes against a reference haplotype panel represented as PBWT using Li's forward-backward 
algorithm combined with Boyer-Moore skip logic, originally proposed by T. Gagie as BML (Boyer-Moore-Li) for BWT,
here adapted to PBWT. For forward-backward; PBML uses longest common prefix (LCP) for forward and longest common 
suffix (LCS) queries for backward matching. 
To materialize the LCP and LCS, PBML needs two PBWTs: a forward PBWT built left-to-right, and a reverse PBWT 
built right-to-left. Both are stored using RLE (run-length encoded) data structures. LCP queries 
on the forward PBWT (extending right) and LCS queries on the reverse PBWT (extending left) to do forward-backward 
using rank queries for interval tracking. Boyer-Moore skipping logic is applied to ensure a site is not
traversed again and again.
PBML has a space complexity of O(R), R being the number of runs in a PBWT.
We use φ (phi) operations to recover prefix array indices of SMEMs from the run-length encoded structures.
The phi_with_hint mechanism is inspired by the MOVI architecture (Zakeri et al.):
    https://www.cell.com/iscience/fulltext/S2589-0042(24)02691-9
which builds on the Move structure by Nishimoto and Tabei:
    https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ICALP.2021.101
PBML Features:
  - Supports *kL*-SMEMs: SMEMs of minimum length *L* occurring at least *k* times in the panel, key in IBD analysis.
  - Index serialization for fast repeated queries, single index can support all query parameterization.
All persistent data structures are RLE-compressed. The temporary all_columns array (used during construction
to build the reverse PBWT from memory) is freed after construction and not required for querying.
================================================================================
*/
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
#include <unordered_map>
#include <memory>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <omp.h>

// ============================================================
// SECTION: Utility — Peek haplotype count from VCF/BCF header
// ============================================================
static size_t peekHaplotypeCount(const std::string &filename)
{
    htsFile *fp = hts_open(filename.c_str(), "rb");
    if (!fp)
        throw std::runtime_error("Failed to open VCF for peeking: " + filename);
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr)
    {
        hts_close(fp);
        throw std::runtime_error("Failed to read header for peeking: " + filename);
    }
    size_t n_samples = bcf_hdr_nsamples(hdr);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    return n_samples * 2;
}

// Peek haplotype count from a PBML index file header
static size_t peekIndexHaplotypeCount(const std::string &filename)
{
    std::ifstream in(filename, std::ios::binary);
    if (!in)
        throw std::runtime_error("Failed to open index for peeking: " + filename);
    uint32_t magic, version;
    in.read(reinterpret_cast<char *>(&magic), sizeof(magic));
    in.read(reinterpret_cast<char *>(&version), sizeof(version));
    if (magic != 0x4C4D4250)
        throw std::runtime_error("Invalid index file (bad magic number)");
    size_t n_sites, n_haplotypes;
    in.read(reinterpret_cast<char *>(&n_sites), sizeof(n_sites));
    in.read(reinterpret_cast<char *>(&n_haplotypes), sizeof(n_haplotypes));
    in.close();
    return n_haplotypes;
}

// ============================================================
// SECTION: Command Line Parsing
// ============================================================
struct CLIOptions {
    std::string mode;                        // "index", "query", "run"
    std::string panel_file;                  // -p, --panel
    std::string query_file;                  // -q, --query
    std::string index_file;                  // -i, --index
    std::string output_file = "smems.tsv";   // -o, --output
    int L = 1;                               // -L, --length
    int k = 1;                               // -k, --min-occ
    int threads = 0;                         // -t, --threads (0 = auto)
    bool verbose = false;                    // -v, --verbose
    bool help = false;                       // -h, --help
};
void printUsage(const char* progname) {
    std::cerr << "PBML - Positional Boyer-Moore-Li for kL-SMEM finding\n\n"
              << "Usage: " << progname << " <mode> [options]\n\n"
              << "Modes:\n"
              << "  index    Build and save index from panel\n"
              << "  query    Query from saved index\n"
              << "  run      Build and query (no index saved)\n\n"
              << "Options:\n"
              << "  -p, --panel <file>     Panel VCF/BCF file (index/run modes)\n"
              << "  -q, --query <file>     Query VCF/BCF file (query/run modes)\n"
              << "  -i, --index <file>     Index file (query: input, index: output)\n"
              << "  -o, --output <file>    Output file [default: smems.tsv]\n"
              << "  -L, --length <int>     Minimum SMEM length [default: 1]\n"
              << "  -k, --min-occ <int>    Minimum SMEM occurrences in panel [default: 1]\n"
              << "  -t, --threads <int>    Number of threads [default: all]\n"
              << "  -v, --verbose          Verbose output\n"
              << "  -h, --help             Show this help\n\n"
              << "Examples:\n"
              << "  " << progname << " index -p panel.bcf -i panel.pbml\n"
              << "  " << progname << " query -i panel.pbml -q query.bcf -L 5 -o smems.tsv\n"
              << "  " << progname << " run -p panel.bcf -q query.bcf -L 5 -k 1\n";
}
CLIOptions parseArgs(int argc, char* argv[]) {
    CLIOptions opts;
    if (argc < 2) {
        opts.help = true;
        return opts;
    }
    std::string first_arg = argv[1];
    if (first_arg == "-h" || first_arg == "--help") {
        opts.help = true;
        return opts;
    }
    opts.mode = first_arg;
    if (opts.mode != "index" && opts.mode != "query" && opts.mode != "run") {
        throw std::runtime_error("Unknown mode: " + opts.mode + ". Use 'index', 'query', or 'run'.");
    }
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            opts.help = true;
            return opts;
        }
        else if (arg == "-v" || arg == "--verbose") {
            opts.verbose = true;
        }
        else if (arg == "-p" || arg == "--panel") {
            if (++i >= argc) throw std::runtime_error("Missing value for " + arg);
            opts.panel_file = argv[i];
        }
        else if (arg == "-q" || arg == "--query") {
            if (++i >= argc) throw std::runtime_error("Missing value for " + arg);
            opts.query_file = argv[i];
        }
        else if (arg == "-i" || arg == "--index") {
            if (++i >= argc) throw std::runtime_error("Missing value for " + arg);
            opts.index_file = argv[i];
        }
        else if (arg == "-o" || arg == "--output") {
            if (++i >= argc) throw std::runtime_error("Missing value for " + arg);
            opts.output_file = argv[i];
        }
        else if (arg == "-L" || arg == "--length") {
            if (++i >= argc) throw std::runtime_error("Missing value for " + arg);
            opts.L = std::stoi(argv[i]);
        }
        else if (arg == "-k" || arg == "--min-occ") {
            if (++i >= argc) throw std::runtime_error("Missing value for " + arg);
            opts.k = std::stoi(argv[i]);
        }
        else if (arg == "-t" || arg == "--threads") {
            if (++i >= argc) throw std::runtime_error("Missing value for " + arg);
            opts.threads = std::stoi(argv[i]);
        }
        else {
            throw std::runtime_error("Unknown option: " + arg);
        }
    }
    return opts;
}
void validateOptions(const CLIOptions& opts) {
    if (opts.mode == "index") {
        if (opts.panel_file.empty())
            throw std::runtime_error("index mode requires -p/--panel");
        if (opts.index_file.empty())
            throw std::runtime_error("index mode requires -i/--index for output");
    }
    else if (opts.mode == "query") {
        if (opts.index_file.empty())
            throw std::runtime_error("query mode requires -i/--index");
        if (opts.query_file.empty())
            throw std::runtime_error("query mode requires -q/--query");
    }
    else if (opts.mode == "run") {
        if (opts.panel_file.empty())
            throw std::runtime_error("run mode requires -p/--panel");
        if (opts.query_file.empty())
            throw std::runtime_error("run mode requires -q/--query");
    }
    if (opts.L <= 0)
        throw std::runtime_error("L must be positive");
    if (opts.k <= 0)
        throw std::runtime_error("k must be positive");
}

// ============================================================
// SECTION: PBWT Column Class (templated on HapT)
// ============================================================
template <typename HapT>
class PBWTColumn
{
public:
    // Index file constants
    static constexpr uint32_t PBML_MAGIC = 0x4C4D4250;  // "PBML" in little-endian
    static constexpr uint32_t PBML_VERSION = 2;          // v2: variable-width haplotype type
    static constexpr uint8_t HAP_BYTES = sizeof(HapT);

    // Return type for LCP function
    struct BMLreturn
    {
        unsigned int lce;
        std::pair<HapT, HapT> interval;
        int interval_col;
        HapT interval_end_hap;
    };

    // φ hint for each successor entry
    struct PhiInfo
    {
        HapT hap;
        uint32_t idx;
    };

    // Successor data structure for phi operations
    struct SuccessorInfo
    {
        uint32_t *sites;
        HapT *predecessors;
        uint32_t *pred_hints;
        size_t count;
        size_t capacity;
        SuccessorInfo() : sites(nullptr), predecessors(nullptr), pred_hints(nullptr), count(0), capacity(0) {}
    };

    // Combined rank + bit result for forward PBWT queries
    struct RankBitResult
    {
        int rank_start;
        int rank_end;
        bool bit_end;
    };

private:
    // RLE structures for forward PBWT
    std::unique_ptr<HapT[]> runLens_fwd;
    std::unique_ptr<char[]> startBits_fwd;
    std::unique_ptr<uint32_t[]> colPtrs_fwd;
    std::unique_ptr<HapT[]> colCs_fwd;
    size_t total_runs_fwd;

    // RLE structures for reverse PBWT
    std::unique_ptr<HapT[]> runLens_rev;
    std::unique_ptr<char[]> startBits_rev;
    std::unique_ptr<uint32_t[]> colPtrs_rev;
    std::unique_ptr<HapT[]> colCs_rev;
    size_t total_runs_rev;

    // Essential structures
    size_t n_sites;
    size_t n_haplotypes;
    std::vector<sdsl::bit_vector> queries;
    std::vector<sdsl::bit_vector> all_columns;
    std::unique_ptr<HapT[]> run_begin_positions;
    std::unique_ptr<PhiInfo[]> run_phi_info;
    std::unique_ptr<HapT[]> end_prefs;

    SuccessorInfo *hap_successor;
    uint32_t *hap_successor_sites;
    HapT *hap_successor_preds;
    uint32_t *hap_successor_hints;

    size_t run_capacity;
    size_t run_count;
    bool verbose;

public:
    // Default constructor
    PBWTColumn()
        : total_runs_fwd(0), total_runs_rev(0),
          n_sites(0), n_haplotypes(0), hap_successor(nullptr),
          hap_successor_sites(nullptr), hap_successor_preds(nullptr),
          hap_successor_hints(nullptr), run_capacity(0), run_count(0),
          verbose(false)
    {
    }

    // Constructor for "run" mode (build + query)
    PBWTColumn(const std::string &panel_file, const std::string &query_file,
               const int L, const unsigned int k, const std::string &output_file,
               bool verbose_flag = false)
        : total_runs_fwd(0), total_runs_rev(0),
          n_sites(0), n_haplotypes(0), hap_successor(nullptr),
          hap_successor_sites(nullptr), hap_successor_preds(nullptr),
          hap_successor_hints(nullptr), run_capacity(0), run_count(0),
          verbose(verbose_flag)
    {
        if (L <= 0)
            throw std::invalid_argument("L must be positive");
        auto start_total = std::chrono::high_resolution_clock::now();
        buildFromPanel(panel_file);
        processQueryFile(query_file);
        processBML(static_cast<unsigned>(L), k, output_file);
        auto end_total = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> total_time = end_total - start_total;
        std::cout << "Total time: " << total_time.count() / 1000.0 << "s\n";
    }

    ~PBWTColumn()
    {
        if (hap_successor)
            delete[] hap_successor;
        if (hap_successor_sites)
            delete[] hap_successor_sites;
        if (hap_successor_preds)
            delete[] hap_successor_preds;
        if (hap_successor_hints)
            delete[] hap_successor_hints;
    }

    // ============================================================
    // SECTION: Public Interface
    // ============================================================
    void setVerbose(bool v) { verbose = v; }

    void buildFromPanel(const std::string &filename)
    {
        auto start = std::chrono::high_resolution_clock::now();
        readAndStoreAllColumns(filename);
        if (n_haplotypes > static_cast<size_t>(std::numeric_limits<HapT>::max()))
            throw std::runtime_error("Panel has " + std::to_string(n_haplotypes) +
                                     " haplotypes but HapT can only hold up to " +
                                     std::to_string(std::numeric_limits<HapT>::max()));
        end_prefs = std::make_unique<HapT[]>(n_sites);
        colPtrs_fwd = std::make_unique<uint32_t[]>(n_sites + 2);
        colCs_fwd = std::make_unique<HapT[]>(n_sites);
        startBits_fwd = std::make_unique<char[]>(n_sites);
        colPtrs_rev = std::make_unique<uint32_t[]>(n_sites + 2);
        colCs_rev = std::make_unique<HapT[]>(n_sites);
        startBits_rev = std::make_unique<char[]>(n_sites);
        buildForwardPBWT();
        buildReversePBWT();
        all_columns.clear();
        all_columns.shrink_to_fit();
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;
        std::cout << "Built PBWTs of " << n_haplotypes << " haplotypes and "
                  << n_sites << " variable sites in " << elapsed.count() / 1000.0 << "s\n";
        std::cout << "Forward PBWT: " << total_runs_fwd << " runs\n";
        std::cout << "Reverse PBWT: " << total_runs_rev << " runs\n";
    }

    void processQueryFile(const std::string &query_file)
    {
        htsFile *qfp = hts_open(query_file.c_str(), "rb");
        if (!qfp)
            throw std::runtime_error("Failed to open query VCF: " + query_file);
        bcf_hdr_t *qhdr = bcf_hdr_read(qfp);
        if (!qhdr)
        {
            hts_close(qfp);
            throw std::runtime_error("Failed to read query header.");
        }
        bcf1_t *qrec = bcf_init();
        if (!qrec)
        {
            bcf_hdr_destroy(qhdr);
            hts_close(qfp);
            throw std::runtime_error("Failed to init query record.");
        }
        std::vector<std::vector<char>> variants_by_site;
        int32_t *gt_arr = nullptr;
        int ngt_arr = 0;
        while (bcf_read(qfp, qhdr, qrec) >= 0)
        {
            bcf_unpack(qrec, BCF_UN_ALL);
            int ngt = bcf_get_genotypes(qhdr, qrec, &gt_arr, &ngt_arr);
            if (ngt <= 0)
                continue;
            int nsmpl = bcf_hdr_nsamples(qhdr);
            int max_ploidy = ngt / nsmpl;
            size_t num_haplotypes = static_cast<size_t>(nsmpl) * static_cast<size_t>(max_ploidy);
            std::vector<char> variant_row(num_haplotypes, 0);
            for (int i = 0; i < nsmpl; ++i)
            {
                int32_t *ptr = gt_arr + i * max_ploidy;
                for (int j = 0; j < max_ploidy; ++j)
                {
                    if (ptr[j] == bcf_int32_vector_end)
                        break;
                    if (bcf_gt_is_missing(ptr[j]))
                    {
                        throw std::runtime_error("Missing allele found in query VCF");
                    }
                    if (bcf_gt_allele(ptr[j]) == 1)
                        variant_row[i * max_ploidy + j] = 1;
                }
            }
            variants_by_site.push_back(std::move(variant_row));
        }
        if (!variants_by_site.empty())
        {
            size_t num_haplotypes = variants_by_site[0].size();
            queries.resize(num_haplotypes);
            for (size_t hap_idx = 0; hap_idx < num_haplotypes; ++hap_idx)
            {
                sdsl::bit_vector hap_bv(variants_by_site.size(), 0);
                for (size_t var_idx = 0; var_idx < variants_by_site.size(); ++var_idx)
                {
                    if (variants_by_site[var_idx][hap_idx])
                        hap_bv[var_idx] = 1;
                }
                queries[hap_idx] = std::move(hap_bv);
            }
        }
        if (gt_arr)
            free(gt_arr);
        bcf_destroy(qrec);
        bcf_hdr_destroy(qhdr);
        hts_close(qfp);
        std::cout << "Loaded " << queries.size() << " query haplotypes\n";
    }

    void processBML(unsigned int L, unsigned int k, const std::string &output_file)
    {
        if (queries.empty())
        {
            std::cerr << "ERROR: No queries loaded!\n";
            return;
        }
        auto start = std::chrono::high_resolution_clock::now();
        const size_t num_threads = static_cast<size_t>(std::max(1, omp_get_max_threads()));
        int n = static_cast<int>(queries.size());
        int B = 2;
        int chunk_size = (1 > n / (num_threads + B) ? 1 : n / (num_threads + B));

        if (k > 1)
        {
            // Streaming mode: per-thread temp files, μ-PBWT-style format
            // Format per SMEM: start, length, [hap1 hap2 hap3 ...]
            std::vector<std::string> tmp_paths(num_threads);
            std::vector<FILE *> tmp_files(num_threads, nullptr);
            for (size_t t = 0; t < num_threads; ++t)
            {
                tmp_paths[t] = output_file + ".tmp." + std::to_string(t);
                tmp_files[t] = std::fopen(tmp_paths[t].c_str(), "wb");
                if (!tmp_files[t])
                    throw std::runtime_error("Failed to open temp file: " + tmp_paths[t]);
                std::setvbuf(tmp_files[t], nullptr, _IOFBF, 256 * 1024);
            }

#pragma omp parallel for schedule(dynamic, chunk_size)
            for (int i = 0; i < n; ++i)
            {
                int thread_id = omp_get_thread_num();
                BML_streaming(queries[i], L, k, static_cast<size_t>(i), tmp_files[thread_id]);
            }

            // Concatenate temp files into final output
            std::ofstream out(output_file, std::ios::binary);
            out.rdbuf()->pubsetbuf(nullptr, 65536);
            std::vector<char> copy_buf(256 * 1024);
            for (size_t t = 0; t < num_threads; ++t)
            {
                std::fclose(tmp_files[t]);
                FILE *in = std::fopen(tmp_paths[t].c_str(), "rb");
                if (in)
                {
                    size_t bytes_read;
                    while ((bytes_read = std::fread(copy_buf.data(), 1, copy_buf.size(), in)) > 0)
                        out.write(copy_buf.data(), static_cast<std::streamsize>(bytes_read));
                    std::fclose(in);
                }
                std::remove(tmp_paths[t].c_str());
            }
            out.close();
        }
        else
        {
            // Buffered mode: k=1, original format (one line per haplotype)
            std::vector<std::string> thread_outputs(num_threads);
            for (auto &output : thread_outputs)
                output.reserve(512 * 1024);

#pragma omp parallel for schedule(dynamic, chunk_size)
            for (int i = 0; i < n; ++i)
            {
                int thread_id = omp_get_thread_num();
                BML_buffered(queries[i], L, k, static_cast<size_t>(i), thread_outputs[thread_id]);
            }

            std::ofstream out(output_file, std::ios::binary);
            out.rdbuf()->pubsetbuf(nullptr, 65536);
            for (const auto &thread_output : thread_outputs)
            {
                if (!thread_output.empty())
                    out.write(thread_output.data(), static_cast<std::streamsize>(thread_output.size()));
            }
            out.close();
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;
        std::cout << "Queried " << queries.size() << " haplotypes in " << elapsed.count() / 1000.0 << "s\n";
    }

    // ============================================================
    // SECTION: Index Serialization
    // ============================================================
    void saveIndex(const std::string& filename) const {
        std::ofstream out(filename, std::ios::binary);
        if (!out)
            throw std::runtime_error("Failed to open index file for writing: " + filename);
        // Header
        out.write(reinterpret_cast<const char*>(&PBML_MAGIC), sizeof(PBML_MAGIC));
        out.write(reinterpret_cast<const char*>(&PBML_VERSION), sizeof(PBML_VERSION));
        // Metadata
        out.write(reinterpret_cast<const char*>(&n_sites), sizeof(n_sites));
        out.write(reinterpret_cast<const char*>(&n_haplotypes), sizeof(n_haplotypes));
        out.write(reinterpret_cast<const char*>(&total_runs_fwd), sizeof(total_runs_fwd));
        out.write(reinterpret_cast<const char*>(&total_runs_rev), sizeof(total_runs_rev));
        // Store HapT size so loader can verify
        uint8_t hap_bytes = HAP_BYTES;
        out.write(reinterpret_cast<const char*>(&hap_bytes), sizeof(hap_bytes));
        // Compute total successors
        size_t total_successors = 0;
        for (size_t i = 0; i < n_haplotypes; ++i)
            total_successors += hap_successor[i].count;
        out.write(reinterpret_cast<const char*>(&total_successors), sizeof(total_successors));
        // Forward PBWT
        out.write(reinterpret_cast<const char*>(runLens_fwd.get()), total_runs_fwd * sizeof(HapT));
        out.write(reinterpret_cast<const char*>(startBits_fwd.get()), n_sites * sizeof(char));
        out.write(reinterpret_cast<const char*>(colPtrs_fwd.get()), (n_sites + 2) * sizeof(uint32_t));
        out.write(reinterpret_cast<const char*>(colCs_fwd.get()), n_sites * sizeof(HapT));
        // Reverse PBWT
        out.write(reinterpret_cast<const char*>(runLens_rev.get()), total_runs_rev * sizeof(HapT));
        out.write(reinterpret_cast<const char*>(startBits_rev.get()), n_sites * sizeof(char));
        out.write(reinterpret_cast<const char*>(colPtrs_rev.get()), (n_sites + 2) * sizeof(uint32_t));
        out.write(reinterpret_cast<const char*>(colCs_rev.get()), n_sites * sizeof(HapT));
        // Phi structures
        out.write(reinterpret_cast<const char*>(run_begin_positions.get()), total_runs_fwd * sizeof(HapT));
        // Write PhiInfo field by field to avoid struct padding issues
        for (size_t i = 0; i < total_runs_fwd; ++i) {
            out.write(reinterpret_cast<const char*>(&run_phi_info[i].hap), sizeof(HapT));
            out.write(reinterpret_cast<const char*>(&run_phi_info[i].idx), sizeof(uint32_t));
        }
        out.write(reinterpret_cast<const char*>(end_prefs.get()), n_sites * sizeof(HapT));
        // Successor arrays - interleaved per-haplotype format
        for (size_t i = 0; i < n_haplotypes; ++i) {
            size_t count = hap_successor[i].count;
            out.write(reinterpret_cast<const char*>(&count), sizeof(count));
            out.write(reinterpret_cast<const char*>(hap_successor[i].sites), count * sizeof(uint32_t));
            out.write(reinterpret_cast<const char*>(hap_successor[i].predecessors), count * sizeof(HapT));
            out.write(reinterpret_cast<const char*>(hap_successor[i].pred_hints), count * sizeof(uint32_t));
        }
        size_t file_size = out.tellp();
        out.close();
        std::cout << "Saved index to " << filename << " (" << file_size / (1024.0 * 1024.0) << " MB)"
                  << " [HapT=" << (sizeof(HapT) * 8) << "bit]\n";
    }

    void loadIndex(const std::string& filename) {
        auto start = std::chrono::high_resolution_clock::now();
        std::ifstream in(filename, std::ios::binary);
        if (!in)
            throw std::runtime_error("Failed to open index file: " + filename);
        // Header
        uint32_t magic, version;
        in.read(reinterpret_cast<char*>(&magic), sizeof(magic));
        in.read(reinterpret_cast<char*>(&version), sizeof(version));
        if (magic != PBML_MAGIC)
            throw std::runtime_error("Invalid index file (bad magic number)");
        if (version != PBML_VERSION)
            throw std::runtime_error("Unsupported index version: " + std::to_string(version));
        // Metadata
        size_t total_successors;
        in.read(reinterpret_cast<char*>(&n_sites), sizeof(n_sites));
        in.read(reinterpret_cast<char*>(&n_haplotypes), sizeof(n_haplotypes));
        in.read(reinterpret_cast<char*>(&total_runs_fwd), sizeof(total_runs_fwd));
        in.read(reinterpret_cast<char*>(&total_runs_rev), sizeof(total_runs_rev));
        // Verify HapT matches
        uint8_t stored_hap_bytes;
        in.read(reinterpret_cast<char*>(&stored_hap_bytes), sizeof(stored_hap_bytes));
        if (stored_hap_bytes != HAP_BYTES)
            throw std::runtime_error("Index was built with " + std::to_string(stored_hap_bytes * 8) +
                                     "-bit haplotype type but loaded with " +
                                     std::to_string(HAP_BYTES * 8) + "-bit type");
        in.read(reinterpret_cast<char*>(&total_successors), sizeof(total_successors));
        // Allocate arrays
        runLens_fwd = std::make_unique<HapT[]>(total_runs_fwd);
        startBits_fwd = std::make_unique<char[]>(n_sites);
        colPtrs_fwd = std::make_unique<uint32_t[]>(n_sites + 2);
        colCs_fwd = std::make_unique<HapT[]>(n_sites);
        runLens_rev = std::make_unique<HapT[]>(total_runs_rev);
        startBits_rev = std::make_unique<char[]>(n_sites);
        colPtrs_rev = std::make_unique<uint32_t[]>(n_sites + 2);
        colCs_rev = std::make_unique<HapT[]>(n_sites);
        run_begin_positions = std::make_unique<HapT[]>(total_runs_fwd);
        run_phi_info = std::make_unique<PhiInfo[]>(total_runs_fwd);
        end_prefs = std::make_unique<HapT[]>(n_sites);
        // Forward PBWT
        in.read(reinterpret_cast<char*>(runLens_fwd.get()), total_runs_fwd * sizeof(HapT));
        in.read(reinterpret_cast<char*>(startBits_fwd.get()), n_sites * sizeof(char));
        in.read(reinterpret_cast<char*>(colPtrs_fwd.get()), (n_sites + 2) * sizeof(uint32_t));
        in.read(reinterpret_cast<char*>(colCs_fwd.get()), n_sites * sizeof(HapT));
        // Reverse PBWT
        in.read(reinterpret_cast<char*>(runLens_rev.get()), total_runs_rev * sizeof(HapT));
        in.read(reinterpret_cast<char*>(startBits_rev.get()), n_sites * sizeof(char));
        in.read(reinterpret_cast<char*>(colPtrs_rev.get()), (n_sites + 2) * sizeof(uint32_t));
        in.read(reinterpret_cast<char*>(colCs_rev.get()), n_sites * sizeof(HapT));
        // Phi structures
        in.read(reinterpret_cast<char*>(run_begin_positions.get()), total_runs_fwd * sizeof(HapT));
        for (size_t i = 0; i < total_runs_fwd; ++i) {
            in.read(reinterpret_cast<char*>(&run_phi_info[i].hap), sizeof(HapT));
            in.read(reinterpret_cast<char*>(&run_phi_info[i].idx), sizeof(uint32_t));
        }
        in.read(reinterpret_cast<char*>(end_prefs.get()), n_sites * sizeof(HapT));
        // Successor arrays
        hap_successor = new SuccessorInfo[n_haplotypes];
        hap_successor_sites = new uint32_t[total_successors];
        hap_successor_preds = new HapT[total_successors];
        hap_successor_hints = new uint32_t[total_successors];
        size_t offset = 0;
        for (size_t i = 0; i < n_haplotypes; ++i) {
            size_t count;
            in.read(reinterpret_cast<char*>(&count), sizeof(count));
            hap_successor[i].sites = hap_successor_sites + offset;
            hap_successor[i].predecessors = hap_successor_preds + offset;
            hap_successor[i].pred_hints = hap_successor_hints + offset;
            hap_successor[i].count = count;
            hap_successor[i].capacity = count;
            in.read(reinterpret_cast<char*>(hap_successor[i].sites), count * sizeof(uint32_t));
            in.read(reinterpret_cast<char*>(hap_successor[i].predecessors), count * sizeof(HapT));
            in.read(reinterpret_cast<char*>(hap_successor[i].pred_hints), count * sizeof(uint32_t));
            offset += count;
        }
        in.close();
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;
        std::cout << "Loaded index: " << n_haplotypes << " haplotypes, "
                  << n_sites << " sites in " << elapsed.count() / 1000.0 << "s"
                  << " [HapT=" << (sizeof(HapT) * 8) << "bit]\n";
        std::cout << "Forward PBWT: " << total_runs_fwd << " runs\n";
        std::cout << "Reverse PBWT: " << total_runs_rev << " runs\n";
    }

private:
    // ============================================================
    // SECTION: PBWT Construction
    // ============================================================

    void buildForwardPBWT()
    {
        std::unique_ptr<HapT[]> pref = std::make_unique<HapT[]>(n_haplotypes);
        std::iota(pref.get(), pref.get() + n_haplotypes, static_cast<HapT>(0));
        std::unique_ptr<char[]> pbwt_column = std::make_unique<char[]>(n_haplotypes);
        std::unique_ptr<HapT[]> temp = std::make_unique<HapT[]>(n_haplotypes);

        // ---- Pass 1: Inline RLE encoding + successor counting ----

        std::vector<size_t> successor_counts(n_haplotypes, 0);
        std::vector<HapT> runLens_vec;
        runLens_vec.reserve(n_sites * 16);

        colPtrs_fwd[0] = 0;

        std::unique_ptr<HapT[]> last_pref;
        sdsl::bit_vector last_bv;

        for (size_t site_idx = 0; site_idx < n_sites; ++site_idx)
        {
            const sdsl::bit_vector &original_column = all_columns[site_idx];

            HapT zero_count = 0;
            for (size_t i = 0; i < n_haplotypes; ++i)
            {
                char val = original_column[pref[i]] ? 1 : 0;
                pbwt_column[i] = val;
                zero_count += (val == 0);
            }

            startBits_fwd[site_idx] = pbwt_column[0];
            colCs_fwd[site_idx] = zero_count;

            successor_counts[pref[0]]++;
            char current_bit = pbwt_column[0];
            HapT run_len = 1;

            for (size_t i = 1; i < n_haplotypes; ++i)
            {
                if (pbwt_column[i] == current_bit)
                {
                    run_len++;
                }
                else
                {
                    runLens_vec.push_back(run_len);
                    successor_counts[pref[i]]++;
                    current_bit = pbwt_column[i];
                    run_len = 1;
                }
            }
            runLens_vec.push_back(run_len);

            colPtrs_fwd[site_idx + 1] = static_cast<uint32_t>(runLens_vec.size());

            updatePBWTState(pref.get(), pbwt_column.get(), temp.get());

            if (site_idx == n_sites - 1)
            {
                last_pref = std::make_unique<HapT[]>(n_haplotypes);
                std::memcpy(last_pref.get(), pref.get(), n_haplotypes * sizeof(HapT));
                last_bv.resize(n_haplotypes);
                for (size_t i = 0; i < n_haplotypes; ++i)
                    last_bv[i] = static_cast<bool>(pbwt_column[i]);
            }
        }

        // Virtual final column RLE + successor counting
        {
            char current_bit = last_bv[0] ? 1 : 0;
            HapT run_len = 1;
            successor_counts[last_pref[0]]++;

            for (size_t idx = 1; idx < n_haplotypes; ++idx)
            {
                char bit = last_bv[idx] ? 1 : 0;
                if (bit == current_bit)
                {
                    run_len++;
                }
                else
                {
                    runLens_vec.push_back(run_len);
                    current_bit = bit;
                    run_len = 1;
                }
                successor_counts[last_pref[idx]]++;
            }
            runLens_vec.push_back(run_len);
            colPtrs_fwd[n_sites + 1] = static_cast<uint32_t>(runLens_vec.size());
        }

        // Move to final flat array
        total_runs_fwd = runLens_vec.size();
        runLens_fwd = std::make_unique<HapT[]>(total_runs_fwd);
        std::memcpy(runLens_fwd.get(), runLens_vec.data(), total_runs_fwd * sizeof(HapT));

        // Allocate successor storage
        size_t total_successors = 0;
        for (size_t i = 0; i < n_haplotypes; ++i)
            total_successors += successor_counts[i];

        hap_successor = new SuccessorInfo[n_haplotypes];
        hap_successor_sites = new uint32_t[total_successors];
        hap_successor_preds = new HapT[total_successors];
        hap_successor_hints = new uint32_t[total_successors];

        size_t offset = 0;
        for (size_t i = 0; i < n_haplotypes; ++i)
        {
            hap_successor[i].sites = hap_successor_sites + offset;
            hap_successor[i].predecessors = hap_successor_preds + offset;
            hap_successor[i].pred_hints = hap_successor_hints + offset;
            hap_successor[i].count = 0;
            hap_successor[i].capacity = successor_counts[i];
            offset += successor_counts[i];
        }

        // Allocate phi structures
        run_capacity = total_runs_fwd;
        run_begin_positions = std::make_unique<HapT[]>(run_capacity);
        run_phi_info = std::make_unique<PhiInfo[]>(run_capacity);

        // ---- Pass 2: Populate phi structures and successor arrays ----

        std::iota(pref.get(), pref.get() + n_haplotypes, static_cast<HapT>(0));
        run_count = 0;

        for (size_t site_idx = 0; site_idx < n_sites; ++site_idx)
        {
            const sdsl::bit_vector &original_column = all_columns[site_idx];
            sdsl::bit_vector bv(n_haplotypes, 0);

            for (size_t i = 0; i < n_haplotypes; ++i)
            {
                bool val = original_column[pref[i]];
                pbwt_column[i] = val;
                bv[i] = val;
            }

            run_begin_positions[run_count] = 0;
            run_phi_info[run_count] = {pref[0], static_cast<uint32_t>(hap_successor[pref[0]].count)};
            run_count++;

            for (size_t idx = 1; idx < n_haplotypes; ++idx)
            {
                if (bv[idx] != bv[idx - 1])
                {
                    HapT current_hap = pref[idx];
                    HapT pred_hap = pref[idx - 1];
                    auto &succ = hap_successor[current_hap];

                    if (succ.count >= succ.capacity)
                    {
                        throw std::runtime_error("Successor count exceeded capacity!");
                    }

                    run_begin_positions[run_count] = static_cast<HapT>(idx);
                    run_phi_info[run_count] = {current_hap, static_cast<uint32_t>(succ.count)};
                    run_count++;

                    succ.sites[succ.count] = static_cast<uint32_t>(site_idx);
                    succ.predecessors[succ.count] = pred_hap;
                    succ.pred_hints[succ.count] = 0;
                    succ.count++;
                }
            }

            end_prefs[site_idx] = pref[n_haplotypes - 1];
            updatePBWTState(pref.get(), pbwt_column.get(), temp.get());

            if (site_idx == n_sites - 1)
            {
                last_pref = std::make_unique<HapT[]>(n_haplotypes);
                std::memcpy(last_pref.get(), pref.get(), n_haplotypes * sizeof(HapT));
                last_bv.resize(n_haplotypes);
                for (size_t i = 0; i < n_haplotypes; ++i)
                    last_bv[i] = bv[i];
            }
        }

        // Virtual final column — phi structures + successor entries
        run_begin_positions[run_count] = 0;
        run_phi_info[run_count] = {last_pref[0], static_cast<uint32_t>(hap_successor[last_pref[0]].count)};
        run_count++;

        for (size_t idx = 1; idx < n_haplotypes; ++idx)
        {
            HapT current_hap = last_pref[idx];
            HapT pred_hap = last_pref[idx - 1];
            auto &succ = hap_successor[current_hap];

            if (last_bv[idx] != last_bv[idx - 1])
            {
                run_begin_positions[run_count] = static_cast<HapT>(idx);
                run_phi_info[run_count] = {current_hap, static_cast<uint32_t>(succ.count)};
                run_count++;
            }

            if (succ.count >= succ.capacity)
            {
                throw std::runtime_error("Successor count exceeded capacity in final column!");
            }

            succ.sites[succ.count] = static_cast<uint32_t>(n_sites);
            succ.predecessors[succ.count] = pred_hap;
            succ.pred_hints[succ.count] = 0;
            succ.count++;
        }

        // ---- Hint computation: O(R × log(R/M)) via binary search ----

        for (size_t hap = 0; hap < n_haplotypes; ++hap)
        {
            auto &succ = hap_successor[hap];
            for (size_t i = 0; i < succ.count; ++i)
            {
                uint32_t site = succ.sites[i];
                HapT pred = succ.predecessors[i];
                const auto &pred_succ = hap_successor[pred];

                size_t left = 0;
                size_t right = pred_succ.count;
                while (left < right)
                {
                    size_t mid = left + (right - left) / 2;
                    if (pred_succ.sites[mid] < site)
                        left = mid + 1;
                    else
                        right = mid;
                }
                succ.pred_hints[i] = static_cast<uint32_t>((left < pred_succ.count) ? left : pred_succ.count - 1);
            }
        }
    }

    void buildReversePBWT()
    {
        std::unique_ptr<HapT[]> pref = std::make_unique<HapT[]>(n_haplotypes);
        std::iota(pref.get(), pref.get() + n_haplotypes, static_cast<HapT>(0));
        std::unique_ptr<char[]> pbwt_column = std::make_unique<char[]>(n_haplotypes);
        std::unique_ptr<HapT[]> temp = std::make_unique<HapT[]>(n_haplotypes);

        std::vector<HapT> runLens_vec;
        runLens_vec.reserve(n_sites * 16);
        std::vector<uint32_t> num_runs_per_site(n_sites, 0);

        if (!all_columns.empty())
        {
            for (int i = static_cast<int>(n_sites) - 1; i >= 0; --i)
            {
                size_t site_idx = static_cast<size_t>(i);
                const sdsl::bit_vector &original_column = all_columns[site_idx];

                HapT zero_count = 0;
                for (size_t j = 0; j < n_haplotypes; ++j)
                {
                    char val = original_column[pref[j]] ? 1 : 0;
                    pbwt_column[j] = val;
                    zero_count += (val == 0);
                }

                startBits_rev[site_idx] = pbwt_column[0];
                colCs_rev[site_idx] = zero_count;

                uint32_t col_runs = 0;
                char current_bit = pbwt_column[0];
                HapT run_len = 1;
                for (size_t j = 1; j < n_haplotypes; ++j)
                {
                    if (pbwt_column[j] == current_bit)
                    {
                        run_len++;
                    }
                    else
                    {
                        runLens_vec.push_back(run_len);
                        col_runs++;
                        current_bit = pbwt_column[j];
                        run_len = 1;
                    }
                }
                runLens_vec.push_back(run_len);
                col_runs++;
                num_runs_per_site[site_idx] = col_runs;

                updatePBWTState(pref.get(), pbwt_column.get(), temp.get());
            }
        }

        // Rearrange from processing order (right-to-left) to site order (left-to-right)
        total_runs_rev = runLens_vec.size();
        runLens_rev = std::make_unique<HapT[]>(total_runs_rev);

        colPtrs_rev[0] = 0;
        for (size_t site_idx = 0; site_idx < n_sites; ++site_idx)
        {
            colPtrs_rev[site_idx + 1] = colPtrs_rev[site_idx] + num_runs_per_site[site_idx];
        }
        colPtrs_rev[n_sites + 1] = static_cast<uint32_t>(total_runs_rev);

        size_t src_offset = runLens_vec.size();
        for (size_t site_idx = 0; site_idx < n_sites; ++site_idx)
        {
            uint32_t nr = num_runs_per_site[site_idx];
            src_offset -= nr;
            std::memcpy(runLens_rev.get() + colPtrs_rev[site_idx],
                        runLens_vec.data() + src_offset,
                        nr * sizeof(HapT));
        }
    }

    void readAndStoreAllColumns(const std::string &filename)
    {
        htsFile *fp = hts_open(filename.c_str(), "rb");
        if (!fp)
            throw std::runtime_error("Failed to open panel VCF: " + filename);
        bcf_hdr_t *hdr = bcf_hdr_read(fp);
        if (!hdr)
        {
            hts_close(fp);
            throw std::runtime_error("Failed to read panel header.");
        }
        bcf1_t *rec = bcf_init();
        if (!rec)
        {
            bcf_hdr_destroy(hdr);
            hts_close(fp);
            throw std::runtime_error("Failed to initialize BCF record.");
        }
        const size_t n_samples = bcf_hdr_nsamples(hdr);
        if (n_samples == 0)
        {
            bcf_destroy(rec);
            bcf_hdr_destroy(hdr);
            hts_close(fp);
            throw std::runtime_error("No samples present in panel VCF.");
        }
        n_haplotypes = n_samples * 2;
        all_columns.reserve(100000);
        int32_t *gt_arr = nullptr;
        int ngt_arr = 0;
        while (bcf_read(fp, hdr, rec) >= 0)
        {
            bcf_unpack(rec, BCF_UN_FMT);
            int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
            if (ngt <= 0)
                continue;
            const int per_sample = ngt / static_cast<int>(n_samples);
            sdsl::bit_vector column_bv(n_haplotypes, 0);
            for (size_t sample = 0; sample < n_samples; ++sample)
            {
                int32_t *ptr = gt_arr + sample * per_sample;
                if (per_sample > 0 && bcf_gt_is_missing(ptr[0]) == 0 && bcf_gt_allele(ptr[0]) == 1)
                    column_bv[sample * 2] = 1;
                if (per_sample > 1 && bcf_gt_is_missing(ptr[1]) == 0 && bcf_gt_allele(ptr[1]) == 1)
                    column_bv[sample * 2 + 1] = 1;
            }
            all_columns.emplace_back(std::move(column_bv));
        }
        n_sites = all_columns.size();
        if (gt_arr)
            free(gt_arr);
        bcf_destroy(rec);
        bcf_hdr_destroy(hdr);
        hts_close(fp);
    }

    void updatePBWTState(HapT *pref, const char *column, HapT *temp) const
    {
        size_t count0 = 0;
        for (size_t i = 0; i < n_haplotypes; ++i)
        {
            count0 += (column[i] == 0);
        }

        size_t idx0 = 0, idx1 = count0;
        for (size_t i = 0; i < n_haplotypes; ++i)
        {
            if (column[i] == 0)
            {
                temp[idx0++] = pref[i];
            }
            else
            {
                temp[idx1++] = pref[i];
            }
        }

        std::memcpy(pref, temp, n_haplotypes * sizeof(HapT));
    }

    // ============================================================
    // SECTION: RLE Rank Queries
    // ============================================================

    RankBitResult rank_and_bit_rle_forward(int i_start, int i_end, int j) const
    {
        RankBitResult result = {0, 0, false};
        if (j < 0 || j >= static_cast<int>(n_sites))
            return result;
        bool need_rank_start = (i_start >= 0);
        if (i_start < 0)
            i_start = 0;
        if (i_end >= static_cast<int>(n_haplotypes))
            i_end = static_cast<int>(n_haplotypes) - 1;
        long runLensSum = 0;
        char runBit = startBits_fwd[j];
        long p = static_cast<long>(colPtrs_fwd[j]);
        long p_end = static_cast<long>(colPtrs_fwd[j + 1]);
        int rank = 0;
        bool found_start = !need_rank_start;
        bool found_end = false;
        while (p < p_end)
        {
            long run_len = static_cast<long>(runLens_fwd[p]);
            long next_sum = runLensSum + run_len;
            if (!found_start && i_start < next_sum)
            {
                result.rank_start = rank;
                if (runBit == 1)
                    result.rank_start += (i_start - runLensSum + 1);
                found_start = true;
            }
            if (!found_end && i_end < next_sum)
            {
                result.rank_end = rank;
                if (runBit == 1)
                    result.rank_end += (i_end - runLensSum + 1);
                result.bit_end = (runBit == 1);
                found_end = true;
                break;
            }
            if (runBit == 1)
                rank += static_cast<int>(run_len);
            runLensSum = next_sum;
            runBit = 1 - runBit;
            p++;
        }
        return result;
    }

    std::pair<int, int> rank_pair_rle_reverse(int i_start, int i_end, int j) const
    {
        if (j < 0 || j >= static_cast<int>(n_sites))
            return {0, 0};
        bool need_rank_start = (i_start >= 0);
        if (i_start < 0)
            i_start = 0;
        if (i_end >= static_cast<int>(n_haplotypes))
            i_end = static_cast<int>(n_haplotypes) - 1;
        long runLensSum = 0;
        char runBit = startBits_rev[j];
        long p = static_cast<long>(colPtrs_rev[j]);
        long p_end = static_cast<long>(colPtrs_rev[j + 1]);
        int rank = 0;
        int rank_start = 0;
        int rank_end = 0;
        bool found_start = !need_rank_start;
        bool found_end = false;
        while (p < p_end)
        {
            long run_len = static_cast<long>(runLens_rev[p]);
            long next_sum = runLensSum + run_len;
            if (!found_start && i_start < next_sum)
            {
                rank_start = rank;
                if (runBit == 1)
                    rank_start += (i_start - runLensSum + 1);
                found_start = true;
            }
            if (!found_end && i_end < next_sum)
            {
                rank_end = rank;
                if (runBit == 1)
                    rank_end += (i_end - runLensSum + 1);
                found_end = true;
                break;
            }
            if (runBit == 1)
                rank += static_cast<int>(run_len);
            runLensSum = next_sum;
            runBit = 1 - runBit;
            p++;
        }
        return {rank_start, rank_end};
    }

    // ============================================================
    // SECTION: Phi Operations
    // ============================================================

    std::pair<HapT, uint32_t> phi_with_hint(HapT hap, uint32_t col, uint32_t hint) const
    {
        const auto &succ = hap_successor[hap];
        if (succ.count == 0)
            return {hap, 0};
        uint32_t idx = hint;
        if (idx >= succ.count)
            idx = static_cast<uint32_t>(succ.count - 1);
        while (idx > 0 && succ.sites[idx - 1] >= col)
            --idx;
        while (idx + 1 < succ.count && succ.sites[idx] < col)
            ++idx;
        return {succ.predecessors[idx], idx};
    }

    std::vector<HapT> getOriginalIndicesForInterval(size_t col, HapT start, HapT end, HapT interval_end_hap)
    {
        std::vector<HapT> result(end - start + 1);
        HapT current_hap = interval_end_hap;
        HapT prev_calling_hap = interval_end_hap;
        const auto &initial_succ = hap_successor[interval_end_hap];
        uint32_t prev_result_idx = 0;
        if (initial_succ.count > 0)
        {
            size_t left = 0, right = initial_succ.count;
            while (left < right)
            {
                size_t mid = left + (right - left) / 2;
                if (initial_succ.sites[mid] < col)
                    left = mid + 1;
                else
                    right = mid;
            }
            prev_result_idx = (left < initial_succ.count) ? static_cast<uint32_t>(left) : static_cast<uint32_t>(initial_succ.count - 1);
        }
        for (HapT pos = end;; --pos)
        {
            result[pos - start] = current_hap;
            if (pos == start)
                break;
            HapT calling_hap = current_hap;
            uint32_t hint = hap_successor[prev_calling_hap].pred_hints[prev_result_idx];
            std::pair<HapT, uint32_t> result_pair = phi_with_hint(current_hap, static_cast<uint32_t>(col), hint);
            current_hap = result_pair.first;
            prev_calling_hap = calling_hap;
            prev_result_idx = result_pair.second;
        }
        return result;
    }

    // ============================================================
    // SECTION: LCS/LCP Queries
    // ============================================================

    int LCS(const sdsl::bit_vector &query, int col_idx, unsigned int k = 1)
    {
        if (n_sites == 0 || col_idx < 0 || col_idx >= static_cast<int>(n_sites))
            return 0;
        unsigned int start = 0;
        unsigned int end = static_cast<unsigned int>(n_haplotypes - 1);
        unsigned int lcs = 0;
        for (lcs = 0; col_idx - static_cast<int>(lcs) >= 0; ++lcs)
        {
            int query_pos = col_idx - static_cast<int>(lcs);
            if (query_pos < 0)
                break;
            bool query_bit = query[query_pos];
            auto [rank_start, rank_end] = rank_pair_rle_reverse(
                static_cast<int>(start) - 1,
                static_cast<int>(end),
                query_pos);
            unsigned int new_start, new_end;
            if (query_bit == 0)
            {
                new_start = start - rank_start;
                new_end = (end + 1) - rank_end;
            }
            else
            {
                new_start = colCs_rev[query_pos] + rank_start;
                new_end = colCs_rev[query_pos] + rank_end;
            }
            if (new_start + k > new_end)
                break;
            start = new_start;
            end = new_end - 1u;
        }
        return static_cast<int>(lcs);
    }

    BMLreturn LCP(const sdsl::bit_vector &query, int col_idx, unsigned int k = 1)
    {
        if (n_sites == 0 || col_idx < 0 || col_idx >= static_cast<int>(n_sites))
            return {0, {0, 0}, col_idx, 0};
        unsigned int start = 0;
        unsigned int end = static_cast<unsigned int>(n_haplotypes - 1);
        unsigned int lcp = 0;
        std::pair<unsigned int, unsigned int> interval = {start, end};
        int last_valid_col = col_idx - 1;
        HapT last_pref = end_prefs[col_idx];
        for (lcp = 0; col_idx + static_cast<int>(lcp) < static_cast<int>(n_sites); ++lcp)
        {
            int current_col = col_idx + static_cast<int>(lcp);
            bool query_bit = query[current_col];
            RankBitResult rle_result = rank_and_bit_rle_forward(
                static_cast<int>(start) - 1,
                static_cast<int>(end),
                current_col);
            int rank_start = rle_result.rank_start;
            int rank_end = rle_result.rank_end;
            bool last_bit = rle_result.bit_end;
            unsigned int new_start, new_end;
            if (query_bit == 0)
            {
                new_start = start - rank_start;
                new_end = (end + 1) - rank_end;
            }
            else
            {
                new_start = colCs_fwd[current_col] + rank_start;
                new_end = colCs_fwd[current_col] + rank_end;
            }
            if (new_start + k > new_end)
                break;
            if (last_bit != query_bit)
            {
                uint32_t run_idx = colPtrs_fwd[current_col];
                uint32_t num_runs_in_col = colPtrs_fwd[current_col + 1] - run_idx;
                size_t left = 0;
                size_t right = num_runs_in_col;
                while (left < right)
                {
                    size_t mid = left + (right - left) / 2;
                    if (run_begin_positions[run_idx + mid] <= end)
                        left = mid + 1;
                    else
                        right = mid;
                }
                size_t r = (left > 0) ? left - 1 : 0;
                size_t global_run_idx = run_idx + r;
                PhiInfo phi_info = run_phi_info[global_run_idx];
                const auto &succ = hap_successor[phi_info.hap];
                last_pref = succ.predecessors[phi_info.idx];
            }
            start = new_start;
            end = new_end - 1u;
            interval = {start, end};
            last_valid_col = current_col;
        }
        return {lcp, {static_cast<HapT>(interval.first), static_cast<HapT>(interval.second)},
                last_valid_col, last_pref};
    }

    // ============================================================
    // SECTION: BML Algorithm
    // ============================================================

    /*
        BML_buffered: Original output mode for k=1.
        One line per (query, panel_haplotype) pair, buffered in memory.
        Format: query_idx\thap_idx\tstart\tend\tlength\n
    */
    void BML_buffered(const sdsl::bit_vector &query, unsigned int L, unsigned int k, size_t query_index, std::string &output)
    {
        const int m = static_cast<int>(query.size());
        const int L_int = static_cast<int>(L);
        int left = 0;
        char buffer[128];
        while (left + L_int - 1 < m)
        {
            const int col = left + L_int - 1;
            const int lcs_result = LCS(query, col, k);
            if (lcs_result < L_int)
            {
                const int partial_match_len = lcs_result;
                const int jump = std::max(1, (partial_match_len > 0) ? (L_int - partial_match_len) : L_int);
                left += jump;
                continue;
            }
            const int smem_start = col - lcs_result + 1;
            BMLreturn lcp_result = LCP(query, smem_start, k);
            const int smem_length = static_cast<int>(lcp_result.lce);
            const int smem_end = smem_start + smem_length - 1;
            if (smem_length >= L_int)
            {
                std::vector<HapT> original_indices = getOriginalIndicesForInterval(
                    smem_end + 1,
                    lcp_result.interval.first,
                    lcp_result.interval.second,
                    lcp_result.interval_end_hap);
                for (size_t idx = 0; idx < original_indices.size(); ++idx)
                {
                    HapT original_row = original_indices[idx];
                    int len = std::snprintf(buffer, sizeof(buffer), "%zu\t%u\t%d\t%d\t%d\n",
                                            query_index, static_cast<unsigned>(original_row),
                                            smem_start, smem_end, smem_length);
                    if (len > 0)
                        output.append(buffer, static_cast<size_t>(len));
                }
            }
            left = smem_end - L_int + 2;
        }
    }

    /*
        BML_streaming: Streaming output for k>1.
        One line per SMEM with all haplotypes listed, written directly to FILE*.
        Format (matches μ-PBWT): query_idx\nstart, length, [hap1 hap2 hap3 ]\n
        Eliminates memory buildup from large output at high k values.
    */
    void BML_streaming(const sdsl::bit_vector &query, unsigned int L, unsigned int k, size_t query_index, FILE *fp)
    {
        const int m = static_cast<int>(query.size());
        const int L_int = static_cast<int>(L);
        int left = 0;
        bool header_written = false;
        while (left + L_int - 1 < m)
        {
            const int col = left + L_int - 1;
            const int lcs_result = LCS(query, col, k);
            if (lcs_result < L_int)
            {
                const int partial_match_len = lcs_result;
                const int jump = std::max(1, (partial_match_len > 0) ? (L_int - partial_match_len) : L_int);
                left += jump;
                continue;
            }
            const int smem_start = col - lcs_result + 1;
            BMLreturn lcp_result = LCP(query, smem_start, k);
            const int smem_length = static_cast<int>(lcp_result.lce);
            const int smem_end = smem_start + smem_length - 1;
            if (smem_length >= L_int)
            {
                // Write query header once per query
                if (!header_written)
                {
                    std::fprintf(fp, "%zu\n", query_index);
                    header_written = true;
                }

                std::vector<HapT> original_indices = getOriginalIndicesForInterval(
                    smem_end + 1,
                    lcp_result.interval.first,
                    lcp_result.interval.second,
                    lcp_result.interval_end_hap);

                std::fprintf(fp, "%d, %d, [", smem_start, smem_length);
                for (size_t idx = 0; idx < original_indices.size(); ++idx)
                {
                    std::fprintf(fp, "%u ", static_cast<unsigned>(original_indices[idx]));
                }
                std::fprintf(fp, "]\n");
            }
            left = smem_end - L_int + 2;
        }
    }
};

// ============================================================
// SECTION: Runtime Dispatch
// ============================================================

// Dispatch helpers — each mode runs with the correct HapT
template <typename HapT>
static void runIndex(const CLIOptions &opts)
{
    PBWTColumn<HapT> pbwt;
    pbwt.setVerbose(opts.verbose);
    pbwt.buildFromPanel(opts.panel_file);
    pbwt.saveIndex(opts.index_file);
}

template <typename HapT>
static void runQuery(const CLIOptions &opts)
{
    PBWTColumn<HapT> pbwt;
    pbwt.setVerbose(opts.verbose);
    pbwt.loadIndex(opts.index_file);
    pbwt.processQueryFile(opts.query_file);
    pbwt.processBML(static_cast<unsigned>(opts.L), static_cast<unsigned>(opts.k), opts.output_file);
}

template <typename HapT>
static void runBuildAndQuery(const CLIOptions &opts)
{
    PBWTColumn<HapT> pbwt(opts.panel_file, opts.query_file, opts.L,
                           static_cast<unsigned>(opts.k), opts.output_file, opts.verbose);
}

// ============================================================
// SECTION: Main
// ============================================================
int main(int argc, char* argv[]) {
    try {
        CLIOptions opts = parseArgs(argc, argv);
        if (opts.help) {
            printUsage(argv[0]);
            return 0;
        }
        validateOptions(opts);
        if (opts.threads > 0)
            omp_set_num_threads(opts.threads);

        // Determine haplotype count to choose HapT
        size_t n_haplotypes = 0;
        if (opts.mode == "index" || opts.mode == "run") {
            n_haplotypes = peekHaplotypeCount(opts.panel_file);
        } else if (opts.mode == "query") {
            n_haplotypes = peekIndexHaplotypeCount(opts.index_file);
        }

        const bool use_wide = (n_haplotypes > std::numeric_limits<uint16_t>::max());
        if (use_wide)
            std::cout << "Panel has " << n_haplotypes << " haplotypes → using 32-bit indices\n";
        else
            std::cout << "Panel has " << n_haplotypes << " haplotypes → using 16-bit indices\n";

        if (opts.mode == "index") {
            if (use_wide) runIndex<uint32_t>(opts);
            else          runIndex<uint16_t>(opts);
        }
        else if (opts.mode == "query") {
            if (use_wide) runQuery<uint32_t>(opts);
            else          runQuery<uint16_t>(opts);
        }
        else if (opts.mode == "run") {
            if (use_wide) runBuildAndQuery<uint32_t>(opts);
            else          runBuildAndQuery<uint16_t>(opts);
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
