/*
===============
Intro
===============

This is the implementation of PBML (Positional Boyer-Moore-Li): A method to enumerate SMEMs from external 
query haplotypes against a reference haplotype panel represented as PBWT using Li's forward-backward 
algorithm combined with Boyer-Moore skip logic, originally proposed by T. Gagie as BML (Boyer-Moore-Li) for BWT,
here adapted to PBWT. For forward-backward; PBML uses longest common prefix (LCP) for forward and longest common 
suffix (LCS) queries for backward matching. 

To materialize the LCP and LCS, PBML needs two PBWTs: a forward PBWT built left-to-right, and a reverse PBWT 
built right-to-left. Both are stored using RLE (run-length encoded) data structures. LCP queries 
on the forward PBWT (extending right) and LCS queries on the reverse PBWT (extending left) to do forward-backward 
using rank queries for interval tracking. Boyer-Moore skipping logic is applied to ensure a site is not
traversed again and again.

PBML uses a two-phase query strategy:
  Phase 1: Run BML for all queries in parallel, collecting SMEM descriptors (intervals + positions)
           without resolving haplotype IDs. Only the RLE structures (~80MB) are needed.
  Phase 2: Single left-to-right reconstruction pass over the forward PBWT, replaying the prefix
           array column-by-column. At each column where an SMEM ends, read off the haplotype IDs
           directly from the current permutation. Cost: O(N × M) total, amortized across all queries.

This eliminates the need for φ (phi) operations, successor arrays, and prefix array samples entirely.
The index contains only the RLE-compressed forward and reverse PBWTs.

PBML Features:
  - Finds all SMEMs of variable minimum length L against a reference panel
  - Supports k-SMEMs: matches occurring at least k times in the panel
  - Index serialization for fast repeated queries

All persistent data structures are RLE-compressed. The temporary all_columns array (used during construction
to build the PBWTs from memory) is freed after construction and not required for querying.

================================================================================
TIME COMPLEXITY ANALYSIS
================================================================================

Let:
  N = number of variant sites
  M = number of haplotypes in panel
  R = total number of runs in forward PBWT (R_fwd ≈ R_rev in practice)
  r = R / N = average runs per column
  m = query haplotype length (typically m = N)
  Q = number of query haplotypes
  P = number of parallel threads

CONSTRUCTION:
-------------
  Panel I/O:           O(N × M)       - Reading and storing all columns
  Forward PBWT:        O(N × M)       - Single pass: build RLE structures
  Reverse PBWT:        O(N × M)       - Single pass right-to-left
  ─────────────────────────────────────
  Total Construction:  O(N × M)

QUERY:
------
  Phase 1 - BML (per haplotype):
    LCS (backward extension):  O(L × r) per call, where L = match length found
    LCP (forward extension):   O(L × r) per call
    ─────────────────────────────────────
    Total per query:           O(m × r) expected
    For Q queries with P threads:
    Phase 1 Total:             O(Q × m × r / P)

  Phase 2 - Reconstruction:
    Single pass over N columns: O(N × M) total
    ─────────────────────────────────────
    Phase 2 Total:             O(N × M)

  Total Query Time:            O(Q × m × r / P + N × M)

SPACE COMPLEXITY:
-----------------
  Forward RLE:         O(R) for runLens
  Reverse RLE:         O(R_rev) for runLens_rev
  Column metadata:     O(N) for colPtrs, colCs, startBits
  ─────────────────────────────────────
  Total Index Size:    O(R + N)

INDEX FILE SIZE:
----------------
  Approximately: R × sizeof(uint16_t)       [runLens_fwd]
               + R_rev × sizeof(uint16_t)   [runLens_rev]
               + 2(N+1) × sizeof(uint32_t)  [colPtrs_fwd, colPtrs_rev]
               + 2N × sizeof(uint16_t)      [colCs]
               + 2N × sizeof(char)          [startBits]
  ≈ 2(R + R_rev) + 12N bytes

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
    std::cerr << "PBML - Positional Boyer-Moore-Li for SMEM finding\n\n"
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
              << "  -k, --min-occ <int>    Minimum occurrences [default: 1]\n"
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
    
    // First argument is mode
    opts.mode = first_arg;
    if (opts.mode != "index" && opts.mode != "query" && opts.mode != "run") {
        throw std::runtime_error("Unknown mode: " + opts.mode + ". Use 'index', 'query', or 'run'.");
    }
    
    // Parse remaining arguments
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
// SECTION: PBWT Column Class
// ============================================================

class PBWTColumn
{
public:
    // Index file constants
    static constexpr uint32_t PBML_MAGIC = 0x4C4D4250;  // "PBML" in little-endian
    static constexpr uint32_t PBML_VERSION = 2;          // v2: no phi/successor

    // Return type for LCP function
    struct BMLreturn
    {
        unsigned int lce;
        std::pair<uint16_t, uint16_t> interval;
    };

    // SMEM descriptor collected during Phase 1 (BML), resolved during Phase 2 (reconstruction)
    struct SMEMDescriptor
    {
        uint32_t query_id;
        int32_t smem_start;
        int32_t smem_end;
        uint16_t interval_start;
        uint16_t interval_end;
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
    std::unique_ptr<uint16_t[]> runLens_fwd;
    std::unique_ptr<char[]> startBits_fwd;
    std::unique_ptr<uint32_t[]> colPtrs_fwd;
    std::unique_ptr<uint16_t[]> colCs_fwd;
    size_t total_runs_fwd;

    // RLE structures for reverse PBWT
    std::unique_ptr<uint16_t[]> runLens_rev;
    std::unique_ptr<char[]> startBits_rev;
    std::unique_ptr<uint32_t[]> colPtrs_rev;
    std::unique_ptr<uint16_t[]> colCs_rev;
    size_t total_runs_rev;

    // Essential structures
    size_t n_sites;
    size_t n_haplotypes;
    std::vector<sdsl::bit_vector> queries;
    std::vector<sdsl::bit_vector> all_columns;  // temporary, freed after construction

    bool verbose;

public:
    // Default constructor
    PBWTColumn()
        : total_runs_fwd(0), total_runs_rev(0),
          n_sites(0), n_haplotypes(0),
          verbose(false)
    {
    }

    // Constructor for "run" mode (build + query)
    PBWTColumn(const std::string &panel_file, const std::string &query_file,
               const int L, const unsigned int k, const std::string &output_file,
               bool verbose_flag = false)
        : total_runs_fwd(0), total_runs_rev(0),
          n_sites(0), n_haplotypes(0),
          verbose(verbose_flag)
    {
        if (L <= 0)
            throw std::invalid_argument("L must be positive");

        auto start_total = std::chrono::high_resolution_clock::now();
        buildFromPanel(panel_file);
        auto end_construction = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> construction_time = end_construction - start_total;
        std::cout << "Built PBWTs of " << n_haplotypes << " haplotypes and "
                  << n_sites << " variable sites in " << construction_time.count() / 1000.0 << "s\n";

        processQueryFile(query_file);

        auto start_query = std::chrono::high_resolution_clock::now();
        processBML(static_cast<unsigned>(L), k, output_file);
        auto end_query = std::chrono::high_resolution_clock::now();
        auto end_total = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> query_time = end_query - start_query;
        std::chrono::duration<double, std::milli> total_time = end_total - start_total;

        std::cout << "Queried in " << query_time.count() / 1000.0 << "s\n";
        std::cout << "Total time: " << total_time.count() / 1000.0 << "s\n";
    }

    ~PBWTColumn() = default;

    // ============================================================
    // SECTION: Public Interface
    // ============================================================

    void setVerbose(bool v) { verbose = v; }

    void buildFromPanel(const std::string &filename)
    {
        auto start = std::chrono::high_resolution_clock::now();
        
        readAndStoreAllColumns(filename);

        colPtrs_fwd = std::make_unique<uint32_t[]>(n_sites + 1);
        colCs_fwd = std::make_unique<uint16_t[]>(n_sites);
        startBits_fwd = std::make_unique<char[]>(n_sites);

        colPtrs_rev = std::make_unique<uint32_t[]>(n_sites + 1);
        colCs_rev = std::make_unique<uint16_t[]>(n_sites);
        startBits_rev = std::make_unique<char[]>(n_sites);

        buildForwardPBWT();
        buildReversePBWT();

        // Free raw columns - reconstruction uses forward PBWT RLE only
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

    /*
        Two-phase SMEM finding:
        Phase 1: Parallel BML to collect SMEM descriptors (intervals without haplotype IDs)
        Phase 2: Single-pass prefix array reconstruction to resolve haplotype IDs
    */
    void processBML(unsigned int L, unsigned int k, const std::string &output_file)
    {
        if (queries.empty())
        {
            std::cerr << "ERROR: No queries loaded!\n";
            return;
        }

        auto start = std::chrono::high_resolution_clock::now();

        // ---- Phase 1: Collect SMEM descriptors in parallel ----
        const size_t num_threads = static_cast<size_t>(std::max(1, omp_get_max_threads()));
        std::vector<std::vector<SMEMDescriptor>> thread_smems(num_threads);
        int n = static_cast<int>(queries.size());
        int B = 2;
        int chunk_size = std::max(1, n / static_cast<int>(num_threads + B));

        for (auto &smems : thread_smems)
            smems.reserve(1024);

#pragma omp parallel for schedule(dynamic, chunk_size)
        for (int i = 0; i < n; ++i)
        {
            int thread_id = omp_get_thread_num();
            BML(queries[i], L, k, static_cast<size_t>(i), thread_smems[thread_id]);
        }

        // Sort each thread's SMEMs locally by smem_end
        size_t total_smems = 0;
        for (auto &ts : thread_smems)
        {
            total_smems += ts.size();
            std::sort(ts.begin(), ts.end(),
                      [](const SMEMDescriptor &a, const SMEMDescriptor &b)
                      { return a.smem_end < b.smem_end; });
        }

        auto phase1_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> phase1_time = phase1_end - start;
        std::cout << "Phase 1 (BML): " << total_smems << " SMEMs from "
                  << queries.size() << " queries in " << phase1_time.count() / 1000.0 << "s\n";

        // ---- Phase 2: Prefix array reconstruction ----
        if (total_smems == 0)
        {
            FILE *fp = std::fopen(output_file.c_str(), "wb");
            if (fp) std::fclose(fp);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> elapsed = end - start;
            std::cout << "Queried " << queries.size() << " haplotypes in " << elapsed.count() / 1000.0 << "s\n";
            return;
        }

        FILE *fp = std::fopen(output_file.c_str(), "wb");
        if (!fp)
            throw std::runtime_error("Failed to open output file: " + output_file);

        std::unique_ptr<uint16_t[]> pref = std::make_unique<uint16_t[]>(n_haplotypes);
        std::iota(pref.get(), pref.get() + n_haplotypes, static_cast<uint16_t>(0));
        std::unique_ptr<char[]> pbwt_column = std::make_unique<char[]>(n_haplotypes);
        std::unique_ptr<uint16_t[]> temp = std::make_unique<uint16_t[]>(n_haplotypes);

        // One cursor per thread
        std::vector<size_t> cursors(num_threads, 0);
        size_t resolved = 0;

        for (size_t site_idx = 0; site_idx < n_sites && resolved < total_smems; ++site_idx)
        {
            decompressForwardColumn(site_idx, pbwt_column.get());
            updatePBWTState(pref.get(), pbwt_column.get(), temp.get());

            for (size_t t = 0; t < num_threads; ++t)
            {
                while (cursors[t] < thread_smems[t].size() &&
                       thread_smems[t][cursors[t]].smem_end == static_cast<int32_t>(site_idx))
                {
                    const SMEMDescriptor &smem = thread_smems[t][cursors[t]];
                    int smem_length = smem.smem_end - smem.smem_start + 1;

                    for (uint16_t pos = smem.interval_start; pos <= smem.interval_end; ++pos)
                    {
                        std::fprintf(fp, "%u\t%u\t%d\t%d\t%d\n",
                                     smem.query_id, pref[pos],
                                     smem.smem_start, smem.smem_end, smem_length);
                    }
                    ++cursors[t];
                    ++resolved;
                }
            }
        }

        std::fclose(fp);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> phase2_time = end - phase1_end;
        std::chrono::duration<double, std::milli> elapsed = end - start;
        std::cout << "Phase 2 (reconstruction): " << phase2_time.count() / 1000.0 << "s\n";
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
        
        // Forward PBWT
        out.write(reinterpret_cast<const char*>(runLens_fwd.get()), total_runs_fwd * sizeof(uint16_t));
        out.write(reinterpret_cast<const char*>(startBits_fwd.get()), n_sites * sizeof(char));
        out.write(reinterpret_cast<const char*>(colPtrs_fwd.get()), (n_sites + 1) * sizeof(uint32_t));
        out.write(reinterpret_cast<const char*>(colCs_fwd.get()), n_sites * sizeof(uint16_t));
        
        // Reverse PBWT
        out.write(reinterpret_cast<const char*>(runLens_rev.get()), total_runs_rev * sizeof(uint16_t));
        out.write(reinterpret_cast<const char*>(startBits_rev.get()), n_sites * sizeof(char));
        out.write(reinterpret_cast<const char*>(colPtrs_rev.get()), (n_sites + 1) * sizeof(uint32_t));
        out.write(reinterpret_cast<const char*>(colCs_rev.get()), n_sites * sizeof(uint16_t));
        
        size_t file_size = out.tellp();
        out.close();
        std::cout << "Saved index to " << filename << " (" << file_size / (1024.0 * 1024.0) << " MB)\n";
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
            throw std::runtime_error("Unsupported index version: " + std::to_string(version) +
                                     " (expected " + std::to_string(PBML_VERSION) + ")");
        
        // Metadata
        in.read(reinterpret_cast<char*>(&n_sites), sizeof(n_sites));
        in.read(reinterpret_cast<char*>(&n_haplotypes), sizeof(n_haplotypes));
        in.read(reinterpret_cast<char*>(&total_runs_fwd), sizeof(total_runs_fwd));
        in.read(reinterpret_cast<char*>(&total_runs_rev), sizeof(total_runs_rev));
        
        // Allocate arrays
        runLens_fwd = std::make_unique<uint16_t[]>(total_runs_fwd);
        startBits_fwd = std::make_unique<char[]>(n_sites);
        colPtrs_fwd = std::make_unique<uint32_t[]>(n_sites + 1);
        colCs_fwd = std::make_unique<uint16_t[]>(n_sites);
        
        runLens_rev = std::make_unique<uint16_t[]>(total_runs_rev);
        startBits_rev = std::make_unique<char[]>(n_sites);
        colPtrs_rev = std::make_unique<uint32_t[]>(n_sites + 1);
        colCs_rev = std::make_unique<uint16_t[]>(n_sites);
        
        // Forward PBWT
        in.read(reinterpret_cast<char*>(runLens_fwd.get()), total_runs_fwd * sizeof(uint16_t));
        in.read(reinterpret_cast<char*>(startBits_fwd.get()), n_sites * sizeof(char));
        in.read(reinterpret_cast<char*>(colPtrs_fwd.get()), (n_sites + 1) * sizeof(uint32_t));
        in.read(reinterpret_cast<char*>(colCs_fwd.get()), n_sites * sizeof(uint16_t));
        
        // Reverse PBWT
        in.read(reinterpret_cast<char*>(runLens_rev.get()), total_runs_rev * sizeof(uint16_t));
        in.read(reinterpret_cast<char*>(startBits_rev.get()), n_sites * sizeof(char));
        in.read(reinterpret_cast<char*>(colPtrs_rev.get()), (n_sites + 1) * sizeof(uint32_t));
        in.read(reinterpret_cast<char*>(colCs_rev.get()), n_sites * sizeof(uint16_t));
        
        in.close();
        
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;
        
        std::cout << "Loaded index: " << n_haplotypes << " haplotypes, " 
                  << n_sites << " sites in " << elapsed.count() / 1000.0 << "s\n";
        std::cout << "Forward PBWT: " << total_runs_fwd << " runs\n";
        std::cout << "Reverse PBWT: " << total_runs_rev << " runs\n";
    }

private:
    // ============================================================
    // SECTION: PBWT Construction
    // ============================================================

    /*
        Build forward PBWT with RLE compression.
        Single pass: iterate sites left-to-right, RLE-encode directly into growing arrays.
        No intermediate structures.
        
        Time: O(N × M)
        Space: O(N × M) temporary for all_columns, O(R + N) for final structures
    */
    void buildForwardPBWT()
    {
        std::unique_ptr<uint16_t[]> pref = std::make_unique<uint16_t[]>(n_haplotypes);
        std::iota(pref.get(), pref.get() + n_haplotypes, 0u);
        std::unique_ptr<char[]> pbwt_column = std::make_unique<char[]>(n_haplotypes);
        std::unique_ptr<uint16_t[]> temp = std::make_unique<uint16_t[]>(n_haplotypes);

        // Grow flat RLE array directly — ~160 runs/col × 1M sites ≈ 160M entries
        std::vector<uint16_t> runLens_vec;
        runLens_vec.reserve(n_sites * 16);  // conservative initial estimate

        colPtrs_fwd[0] = 0;

        for (size_t site_idx = 0; site_idx < n_sites; ++site_idx)
        {
            const sdsl::bit_vector &original_column = all_columns[site_idx];

            // Build permuted column + inline RLE extraction
            uint16_t zero_count = 0;
            for (size_t i = 0; i < n_haplotypes; ++i)
            {
                char val = original_column[pref[i]] ? 1 : 0;
                pbwt_column[i] = val;
                zero_count += (val == 0);
            }

            // RLE-encode directly from pbwt_column
            startBits_fwd[site_idx] = pbwt_column[0];
            colCs_fwd[site_idx] = zero_count;

            char current_bit = pbwt_column[0];
            uint16_t run_len = 1;
            for (size_t i = 1; i < n_haplotypes; ++i)
            {
                if (pbwt_column[i] == current_bit)
                {
                    run_len++;
                }
                else
                {
                    runLens_vec.push_back(run_len);
                    current_bit = pbwt_column[i];
                    run_len = 1;
                }
            }
            runLens_vec.push_back(run_len);

            colPtrs_fwd[site_idx + 1] = static_cast<uint32_t>(runLens_vec.size());

            updatePBWTState(pref.get(), pbwt_column.get(), temp.get());
        }

        // Move to final flat array
        total_runs_fwd = runLens_vec.size();
        runLens_fwd = std::make_unique<uint16_t[]>(total_runs_fwd);
        std::memcpy(runLens_fwd.get(), runLens_vec.data(), total_runs_fwd * sizeof(uint16_t));
    }

    /*
        Build reverse PBWT with RLE compression.
        Processes columns right-to-left to enable backward extension (LCS queries).
        Single pass: RLE-encode directly into growing arrays.
        
        Time: O(N × M) - single pass over all columns
        Space: O(R_rev) for RLE structures
    */
    void buildReversePBWT()
    {
        std::unique_ptr<uint16_t[]> pref = std::make_unique<uint16_t[]>(n_haplotypes);
        std::iota(pref.get(), pref.get() + n_haplotypes, 0u);
        std::unique_ptr<char[]> pbwt_column = std::make_unique<char[]>(n_haplotypes);
        std::unique_ptr<uint16_t[]> temp = std::make_unique<uint16_t[]>(n_haplotypes);

        // Process right-to-left, accumulate runs in processing order
        std::vector<uint16_t> runLens_vec;
        runLens_vec.reserve(n_sites * 16);
        std::vector<uint32_t> num_runs_per_site(n_sites, 0);

        if (!all_columns.empty())
        {
            for (int i = static_cast<int>(n_sites) - 1; i >= 0; --i)
            {
                size_t site_idx = static_cast<size_t>(i);
                const sdsl::bit_vector &original_column = all_columns[site_idx];

                uint16_t zero_count = 0;
                for (size_t j = 0; j < n_haplotypes; ++j)
                {
                    char val = original_column[pref[j]] ? 1 : 0;
                    pbwt_column[j] = val;
                    zero_count += (val == 0);
                }

                startBits_rev[site_idx] = pbwt_column[0];
                colCs_rev[site_idx] = zero_count;

                // RLE-encode directly
                uint32_t col_runs = 0;
                char current_bit = pbwt_column[0];
                uint16_t run_len = 1;
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

        // Rearrange: runLens_vec is in reverse site order (N-1 first, 0 last)
        // Build flat array in site order (0 first, N-1 last)
        total_runs_rev = runLens_vec.size();
        runLens_rev = std::make_unique<uint16_t[]>(total_runs_rev);

        // Compute colPtrs_rev in site order
        colPtrs_rev[0] = 0;
        for (size_t site_idx = 0; site_idx < n_sites; ++site_idx)
        {
            colPtrs_rev[site_idx + 1] = colPtrs_rev[site_idx] + num_runs_per_site[site_idx];
        }

        // Copy runs from processing order (reverse) to site order
        // runLens_vec has: [site N-1 runs][site N-2 runs]...[site 0 runs]
        // We need:         [site 0 runs][site 1 runs]...[site N-1 runs]
        size_t src_offset = runLens_vec.size();
        for (size_t site_idx = 0; site_idx < n_sites; ++site_idx)
        {
            uint32_t nr = num_runs_per_site[site_idx];
            src_offset -= nr;
            std::memcpy(runLens_rev.get() + colPtrs_rev[site_idx],
                        runLens_vec.data() + src_offset,
                        nr * sizeof(uint16_t));
        }
    }

    /*
        Read VCF/BCF panel and store all columns in memory.
        Time: O(N × M) for reading + storing
        Space: O(N × M) bits for all_columns (freed after construction)
    */
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

    /*
        Update PBWT prefix array after processing a column.
        Stable sort: zeros before ones, preserving relative order within each group.
        Time: O(M) where M = n_haplotypes
        Space: O(1) — uses caller-provided temp buffer
    */
    void updatePBWTState(uint16_t *pref, const char *column, uint16_t *temp) const
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

        std::memcpy(pref, temp, n_haplotypes * sizeof(uint16_t));
    }

    // ============================================================
    // SECTION: RLE Decompression & Rank Queries
    // ============================================================

    /*
        Decompress a forward PBWT column from RLE to a char array.
        Time: O(M) where M = n_haplotypes
        Space: O(1) (output written to caller-provided buffer)
    */
    void decompressForwardColumn(size_t site_idx, char *pbwt_column) const
    {
        uint32_t p_start = colPtrs_fwd[site_idx];
        uint32_t p_end = colPtrs_fwd[site_idx + 1];
        char bit = startBits_fwd[site_idx];
        size_t pos = 0;

        for (uint32_t r = p_start; r < p_end; ++r)
        {
            uint16_t len = runLens_fwd[r];
            std::memset(pbwt_column + pos, bit, len);
            pos += len;
            bit = 1 - bit;
        }
    }

    /*
        Combined rank and bit query on forward PBWT column j.
        Returns rank_1(i_start), rank_1(i_end), and bit value at i_end in a single scan.
        
        Time: O(r) where r = number of runs in column j
        Space: O(1)
    */
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

    /*
        Combined rank query on reverse PBWT column j.
        Returns rank_1(i_start) and rank_1(i_end) in a single scan.
        
        Time: O(r) where r = number of runs in column j
        Space: O(1)
    */
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
    // SECTION: LCS/LCP Queries
    // ============================================================

    /*
        Longest Common Suffix query on reverse PBWT.
        Extends backward from col_idx to find the longest match with at least k occurrences.
        
        Time: O(L × r) where L = match length found, r = average runs per column
        Space: O(1)
    */
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

    /*
        Longest Common Prefix query on forward PBWT.
        Extends forward from col_idx to find the longest match with at least k occurrences.
        Returns match length and the interval in the forward PBWT.
        
        Time: O(L × r) where L = match length found, r = average runs per column
        Space: O(1)
    */
    BMLreturn LCP(const sdsl::bit_vector &query, int col_idx, unsigned int k = 1)
    {
        if (n_sites == 0 || col_idx < 0 || col_idx >= static_cast<int>(n_sites))
            return {0, {0, 0}};

        unsigned int start = 0;
        unsigned int end = static_cast<unsigned int>(n_haplotypes - 1);
        unsigned int lcp = 0;
        std::pair<unsigned int, unsigned int> interval = {start, end};

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

            start = new_start;
            end = new_end - 1u;
            interval = {start, end};
        }

        return {lcp, interval};
    }

    // ============================================================
    // SECTION: BML Algorithm
    // ============================================================

    /*
        BML (Boyer-Moore-Li) algorithm adapted for PBWT.
        Finds all k-SMEMs of length >= L for a single query haplotype.
        
        Phase 1 only: collects SMEM descriptors without resolving haplotype IDs.
        Uses LCS on reverse PBWT to find potential match start, then LCP on forward 
        PBWT to extend and verify. Skips forward using Boyer-Moore-style jumps.
        
        Time: O(m × r) expected, where m = query length, r = average runs per column
        Space: O(S) where S = number of SMEMs found
    */
    void BML(const sdsl::bit_vector &query, unsigned int L, unsigned int k, size_t query_index,
             std::vector<SMEMDescriptor> &smems)
    {
        const int m = static_cast<int>(query.size());
        const int L_int = static_cast<int>(L);
        int left = 0;

        while (left + L_int - 1 < m)
        {
            const int col = left + L_int - 1;
            const int lcs_result = LCS(query, col, k);

            if (lcs_result < L_int)
            {
                // Boyer-Moore skip: jump forward based on partial match length
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
                smems.push_back({
                    static_cast<uint32_t>(query_index),
                    smem_start,
                    smem_end,
                    lcp_result.interval.first,
                    lcp_result.interval.second
                });
            }

            // Skip past current SMEM, maintaining L-1 overlap for next potential SMEM
            left = smem_end - L_int + 2;
        }
    }
};

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
        
        if (opts.mode == "index") {
            PBWTColumn pbwt;
            pbwt.setVerbose(opts.verbose);
            pbwt.buildFromPanel(opts.panel_file);
            pbwt.saveIndex(opts.index_file);
        }
        else if (opts.mode == "query") {
            PBWTColumn pbwt;
            pbwt.setVerbose(opts.verbose);
            pbwt.loadIndex(opts.index_file);
            pbwt.processQueryFile(opts.query_file);
            pbwt.processBML(static_cast<unsigned>(opts.L), static_cast<unsigned>(opts.k), opts.output_file);
        }
        else if (opts.mode == "run") {
            PBWTColumn pbwt(opts.panel_file, opts.query_file, opts.L, static_cast<unsigned>(opts.k), opts.output_file, opts.verbose);
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
