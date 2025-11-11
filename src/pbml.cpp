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
#include <omp.h>

class PBWTColumn
{
public:
    struct CompressedColumn
    {
        sdsl::sd_vector<> bv;

        CompressedColumn() = default;
        explicit CompressedColumn(const sdsl::bit_vector &input) : bv(input) {}
        CompressedColumn(CompressedColumn &&other) noexcept : bv(std::move(other.bv)) {}

        CompressedColumn &operator=(CompressedColumn &&other) noexcept
        {
            if (this != &other)
            {
                bv = std::move(other.bv);
            }
            return *this;
        }

        CompressedColumn(const CompressedColumn &other) = delete;
        CompressedColumn &operator=(const CompressedColumn &other) = delete;
        ~CompressedColumn() = default;
    };

    struct BMLreturn
    {
        unsigned int lce;
        std::pair<uint16_t, uint16_t> interval;
        int interval_col;
        uint16_t interval_end_hap;
    };
    
    struct PhiInfo
    {
        uint16_t hap;
        uint32_t idx;
    };
    
    struct SuccessorInfo
    {
        uint32_t* sites;
        uint16_t* predecessors;
        uint32_t* pred_hints;
        size_t count;
        
        SuccessorInfo() : sites(nullptr), predecessors(nullptr), pred_hints(nullptr), count(0) {}
    };

private:
    std::unique_ptr<CompressedColumn[]> forward_pbwt;
    std::unique_ptr<CompressedColumn[]> reverse_pbwt;
    std::unique_ptr<uint16_t[]> zero_counts;
    size_t n_sites;
    size_t n_haplotypes;
    std::vector<sdsl::bit_vector> queries;
    std::vector<sdsl::bit_vector> all_columns;
    std::unique_ptr<uint16_t[]> run_begin_positions;
    std::unique_ptr<uint16_t[]> run_beg_prefs;
    std::unique_ptr<PhiInfo[]> run_phi_info;
    std::unique_ptr<uint16_t[]> end_prefs;
    std::unique_ptr<uint32_t[]> col_run_index;
    SuccessorInfo* hap_successor;
    uint32_t* hap_successor_sites;
    uint16_t* hap_successor_preds;
    uint32_t* hap_successor_hints;
    size_t run_capacity;
    size_t run_count;
    size_t max_successor_count;
    
    static unsigned compute_bits(size_t n)
    {
        unsigned b = 1;
        while ((1ULL << b) < n)
            ++b;
        return b;
    }

public:
    PBWTColumn(const std::string &panel_file, const std::string &query_file,
               const int L, const std::string &output_file)
        : n_sites(0), n_haplotypes(0), hap_successor(nullptr), 
          hap_successor_sites(nullptr), hap_successor_preds(nullptr), 
          hap_successor_hints(nullptr), run_capacity(0), run_count(0),
          max_successor_count(0)
    {
        if (L <= 0)
            throw std::invalid_argument("L must be positive");

        auto start_total = std::chrono::high_resolution_clock::now();
        buildPBWTs(panel_file);
        auto end_construction = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> construction_time = end_construction - start_total;
        std::cout << "Built PBWTs of " << n_haplotypes << " haplotypes and "
                  << n_sites << " variable sites in " << construction_time.count() / 1000.0 << "s\n";

        processQueryFile(query_file);

        auto start_query = std::chrono::high_resolution_clock::now();
        processBML(static_cast<unsigned>(L), output_file);
        auto end_query = std::chrono::high_resolution_clock::now();
        auto end_total = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> query_time = end_query - start_query;
        std::chrono::duration<double, std::milli> total_time = end_total - start_total;

        std::cout << "Queried in " << query_time.count() / 1000.0 << "s\n";
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

private:
    void buildPBWTs(const std::string &filename)
    {
        readAndStoreAllColumns(filename);
        forward_pbwt = std::make_unique<CompressedColumn[]>(n_sites);
        reverse_pbwt = std::make_unique<CompressedColumn[]>(n_sites);
        zero_counts = std::make_unique<uint16_t[]>(n_sites);
        end_prefs = std::make_unique<uint16_t[]>(n_sites);
        
        col_run_index = std::make_unique<uint32_t[]>(n_sites + 2);
        
        run_capacity = n_sites * 4 ;
        run_begin_positions = std::make_unique<uint16_t[]>(run_capacity);
        run_beg_prefs = std::make_unique<uint16_t[]>(run_capacity);
        run_phi_info = std::make_unique<PhiInfo[]>(run_capacity);
        
        buildForwardPBWT();
        buildReversePBWT();
        all_columns.clear();
        all_columns.shrink_to_fit();
    }
    
    void resize_runs_if_needed()
    {
        if (run_count >= run_capacity)
        {
            size_t new_capacity = run_capacity * 2;
            auto new_begin = std::make_unique<uint16_t[]>(new_capacity);
            auto new_beg_prefs = std::make_unique<uint16_t[]>(new_capacity);
            auto new_phi_info = std::make_unique<PhiInfo[]>(new_capacity);
            
            std::copy(run_begin_positions.get(), run_begin_positions.get() + run_count, new_begin.get());
            std::copy(run_beg_prefs.get(), run_beg_prefs.get() + run_count, new_beg_prefs.get());
            std::copy(run_phi_info.get(), run_phi_info.get() + run_count, new_phi_info.get());
            
            run_begin_positions = std::move(new_begin);
            run_beg_prefs = std::move(new_beg_prefs);
            run_phi_info = std::move(new_phi_info);
            run_capacity = new_capacity;
        }
    }

    void buildForwardPBWT()
    {
        hap_successor = new SuccessorInfo[n_haplotypes];
        hap_successor_sites = new uint32_t[n_haplotypes * n_sites];
        hap_successor_preds = new uint16_t[n_haplotypes * n_sites];
        hap_successor_hints = new uint32_t[n_haplotypes * n_sites];
        
        for (size_t i = 0; i < n_haplotypes; ++i)
        {
            hap_successor[i].sites = hap_successor_sites + (i * n_sites);
            hap_successor[i].predecessors = hap_successor_preds + (i * n_sites);
            hap_successor[i].pred_hints = hap_successor_hints + (i * n_sites);
            hap_successor[i].count = 0;
        }
        
        std::unique_ptr<uint16_t[]> pref = std::make_unique<uint16_t[]>(n_haplotypes);
        std::iota(pref.get(), pref.get() + n_haplotypes, 0u);
        std::unique_ptr<char[]> pbwt_column = std::make_unique<char[]>(n_haplotypes);
        
        size_t total_runs = 0;
        run_count = 0;
        col_run_index[0] = 0;
        
        std::unique_ptr<uint16_t[]> last_pref;
        sdsl::bit_vector last_bv;

        for (size_t site_idx = 0; site_idx < n_sites; ++site_idx)
        {
            const sdsl::bit_vector &original_column = all_columns[site_idx];
            sdsl::bit_vector bv(n_haplotypes, 0);
            uint16_t zero_count = 0;

            for (size_t i = 0; i < n_haplotypes; ++i)
            {
                bool val = original_column[pref[i]];
                pbwt_column[i] = val;
                bv[i] = val;
                zero_count += !val;
            }

            size_t column_runs = 1;
            resize_runs_if_needed();
            run_begin_positions[run_count] = 0;
            run_beg_prefs[run_count] = pref[0];
            run_phi_info[run_count] = {pref[0], static_cast<uint32_t>(hap_successor[pref[0]].count)};
            run_count++;

            for (size_t idx = 1; idx < n_haplotypes; ++idx)
            {
                if (bv[idx] != bv[idx - 1])
                {
                    uint16_t current_hap = pref[idx];
                    uint16_t pred_hap = pref[idx - 1];
                    
                    auto &succ = hap_successor[current_hap];
                    
                    resize_runs_if_needed();
                    run_begin_positions[run_count] = static_cast<uint16_t>(idx);
                    run_beg_prefs[run_count] = current_hap;
                    run_phi_info[run_count] = {current_hap, static_cast<uint32_t>(succ.count)};
                    run_count++;

                    succ.sites[succ.count] = static_cast<uint32_t>(site_idx);
                    succ.predecessors[succ.count] = pred_hap;
                    succ.pred_hints[succ.count] = 0;
                    succ.count++;
                    
                    column_runs++;
                }
            }

            end_prefs[site_idx] = pref[n_haplotypes - 1];
            
            total_runs += column_runs;
            col_run_index[site_idx + 1] = static_cast<uint32_t>(total_runs);
            forward_pbwt[site_idx] = CompressedColumn(bv);
            zero_counts[site_idx] = zero_count;

            updatePBWTState(pref.get(), pbwt_column.get());

            if (site_idx == n_sites - 1)
            {
                last_pref = std::make_unique<uint16_t[]>(n_haplotypes);
                std::copy(pref.get(), pref.get() + n_haplotypes, last_pref.get());
                last_bv = bv;
            }
        }

        size_t final_column_runs = 1;
        resize_runs_if_needed();
        run_begin_positions[run_count] = 0;
        run_beg_prefs[run_count] = last_pref[0];
        run_phi_info[run_count] = {last_pref[0], static_cast<uint32_t>(hap_successor[last_pref[0]].count)};
        run_count++;

        for (size_t idx = 1; idx < n_haplotypes; ++idx)
        {
            uint16_t current_hap = last_pref[idx];
            uint16_t pred_hap = last_pref[idx - 1];
            
            auto &succ = hap_successor[current_hap];
            
            if (last_bv[idx] != last_bv[idx - 1])
            {
                resize_runs_if_needed();
                run_begin_positions[run_count] = static_cast<uint16_t>(idx);
                run_beg_prefs[run_count] = current_hap;
                run_phi_info[run_count] = {current_hap, static_cast<uint32_t>(succ.count)};
                run_count++;
                final_column_runs++;
            }
            
            succ.sites[succ.count] = static_cast<uint32_t>(n_sites - 1);
            succ.predecessors[succ.count] = pred_hap;
            succ.pred_hints[succ.count] = 0;
            succ.count++;
        }
        
        total_runs += final_column_runs;
        col_run_index[n_sites + 1] = static_cast<uint32_t>(total_runs);
        // Now compute hints with the resized arrays
        for (uint16_t hap = 0; hap < n_haplotypes; ++hap)
        {
            auto &succ = hap_successor[hap];
            
            for (size_t i = 0; i < succ.count; ++i)
            {
                uint32_t site = succ.sites[i];
                uint16_t pred = succ.predecessors[i];

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
        std::unique_ptr<uint16_t[]> pref = std::make_unique<uint16_t[]>(n_haplotypes);
        std::iota(pref.get(), pref.get() + n_haplotypes, 0u);
        std::unique_ptr<char[]> pbwt_column = std::make_unique<char[]>(n_haplotypes);

        if (!all_columns.empty())
        {
            int start_i = static_cast<int>(all_columns.size()) - 1;
            for (int i = start_i; i >= 0; --i)
            {
                const sdsl::bit_vector &original_column = all_columns[static_cast<size_t>(i)];

                for (size_t j = 0; j < n_haplotypes; ++j)
                {
                    pbwt_column[j] = original_column[pref[j]] ? 1 : 0;
                }

                sdsl::bit_vector bv(n_haplotypes, 0);
                for (size_t j = 0; j < n_haplotypes; ++j)
                    bv[j] = static_cast<bool>(pbwt_column[j]);

                reverse_pbwt[i] = CompressedColumn(bv);
                updatePBWTState(pref.get(), pbwt_column.get());
            }
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

    void updatePBWTState(uint16_t* pref, const char* column) const
    {
        size_t count0 = 0;
        for (size_t i = 0; i < n_haplotypes; ++i)
        {
            count0 += (column[i] == 0);
        }

        std::unique_ptr<uint16_t[]> new_pref = std::make_unique<uint16_t[]>(n_haplotypes);
        size_t idx0 = 0, idx1 = count0;

        for (size_t i = 0; i < n_haplotypes; ++i)
        {
            if (column[i] == 0)
            {
                new_pref[idx0++] = pref[i];
            }
            else
            {
                new_pref[idx1++] = pref[i];
            }
        }

        std::copy(new_pref.get(), new_pref.get() + n_haplotypes, pref);
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
    }

    int LCS(const sdsl::bit_vector &query, int col_idx)
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
            const CompressedColumn &column = reverse_pbwt[query_pos];
            sdsl::sd_vector<>::rank_1_type rank1(&column.bv);
            unsigned int new_start, new_end;

            if (query_bit == 0)
            {
                new_start = (start == 0) ? 0u : start - rank1(start);
                new_end = (end + 1u) - rank1(end + 1u);
            }
            else
            {
                new_start = zero_counts[query_pos] + ((start == 0) ? 0u : rank1(start));
                new_end = zero_counts[query_pos] + rank1(end + 1u);
            }

            if (new_start >= new_end)
            {
                break;
            }

            start = new_start;
            end = new_end - 1u;
        }

        return static_cast<int>(lcs);
    }

    BMLreturn LCP(const sdsl::bit_vector &query, int col_idx)
    {
        if (n_sites == 0 || col_idx < 0 || col_idx >= static_cast<int>(n_sites))
            return {0, {0, 0}, col_idx, 0};

        unsigned int start = 0;
        unsigned int end = static_cast<unsigned int>(n_haplotypes - 1);
        unsigned int lcp = 0;
        std::pair<unsigned int, unsigned int> interval = {start, end};
        int last_valid_col = col_idx - 1;
        
        uint16_t last_pref = (col_idx > 0) ? end_prefs[col_idx - 1] : static_cast<uint16_t>(n_haplotypes - 1);

        for (lcp = 0; col_idx + static_cast<int>(lcp) < static_cast<int>(n_sites); ++lcp)
        {
            int current_col = col_idx + static_cast<int>(lcp);
            const CompressedColumn &column = forward_pbwt[current_col];
            sdsl::sd_vector<>::rank_1_type rank1(&column.bv);
            unsigned int new_start, new_end;
            bool query_bit = query[current_col];

            if (query_bit == 0)
            {
                new_start = (start == 0) ? 0u : start - rank1(start);
                new_end = (end + 1u) - rank1(end + 1u);
            }
            else
            {
                new_start = zero_counts[current_col] + ((start == 0) ? 0u : rank1(start));
                new_end = zero_counts[current_col] + rank1(end + 1u);
            }

            if (new_start >= new_end)
                break;

            bool last_bit = column.bv[end];
            
            if (last_bit != query_bit)
            {
                uint32_t run_idx = col_run_index[current_col];
                uint32_t num_runs_in_col = col_run_index[current_col + 1] - run_idx;
                
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

        return {lcp, interval, last_valid_col, last_pref};
    }

    std::pair<uint16_t, uint32_t> phi_with_hint(uint16_t hap, uint32_t col, uint32_t hint) const
    {
        const auto &succ = hap_successor[hap];
        if (succ.count == 0)
            return {hap, 0};
        
        uint32_t idx = hint;
        
        while (idx > 0 && succ.sites[idx - 1] >= col)
        {
            --idx;
        }
        
        if (succ.sites[idx] >= col)
        {
            return {succ.predecessors[idx], idx};
        }
        
        return {succ.predecessors[succ.count - 1], static_cast<uint32_t>(succ.count - 1)};
    }

    std::vector<uint16_t> getOriginalIndicesForInterval(size_t col, uint16_t start, uint16_t end, uint16_t interval_end_hap)
    {
        std::vector<uint16_t> result(end - start + 1);
        
        uint16_t current_hap = interval_end_hap;
        uint16_t prev_calling_hap = interval_end_hap;
        uint32_t prev_result_idx = 0;
        
        uint32_t run_idx = col_run_index[col];
        uint32_t num_runs_in_col = col_run_index[col + 1] - run_idx;
        
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
        
        PhiInfo initial_phi = run_phi_info[global_run_idx];
        prev_result_idx = initial_phi.idx;
        
        for (uint16_t pos = end;; --pos)
        {
            result[pos - start] = current_hap;
            
            if (pos == start)
                break;
            
            uint16_t calling_hap = current_hap;
            uint32_t hint = hap_successor[prev_calling_hap].pred_hints[prev_result_idx];
            std::pair<uint16_t, uint32_t> result_pair = phi_with_hint(current_hap, static_cast<uint32_t>(col), hint);
            
            current_hap = result_pair.first;
            prev_calling_hap = calling_hap;
            prev_result_idx = result_pair.second;
        }
        
        return result;
    }
    
    void processBML(unsigned int L, const std::string &output_file)
    {
        const size_t num_threads = static_cast<size_t>(std::max(1, omp_get_max_threads()));
        std::vector<std::string> thread_outputs(num_threads);
        int n = static_cast<int>(queries.size());
        int B = 2;
        int chunk_size = (1 > n / (num_threads + B) ? 1 : n / (num_threads + B));

        for (auto &output : thread_outputs)
            output.reserve(512 * 1024);

#pragma omp parallel for schedule(dynamic, chunk_size)
        for (int i = 0; i < n; ++i)
        {
            int thread_id = omp_get_thread_num();
            BML(queries[i], L, static_cast<size_t>(i), thread_outputs[thread_id]);
        }

        std::ofstream out(output_file, std::ios::binary);
        out.rdbuf()->pubsetbuf(nullptr, 65536);
        for (const auto &thread_output : thread_outputs)
        {
            if (!thread_output.empty())
            {
                out.write(thread_output.data(), static_cast<std::streamsize>(thread_output.size()));
            }
        }
        out.close();
    }
    
    void BML(const sdsl::bit_vector &query, unsigned int L, size_t query_index, std::string &output)
    {
        const int m = static_cast<int>(query.size());
        const int L_int = static_cast<int>(L);
        int left = 0;
        char buffer[128];

        while (left + L_int - 1 < m)
        {
            const int col = left + L_int - 1;
            const int lcs_result = LCS(query, col);

            if (lcs_result < L_int)
            {
                const int partial_match_len = lcs_result;
                const int jump = std::max(1, (partial_match_len > 0) ? (L_int - partial_match_len) : L_int);
                left += jump;
                continue;
            }

            const int smem_start = col - lcs_result + 1;
            BMLreturn lcp_result = LCP(query, smem_start);
            const int smem_length = static_cast<int>(lcp_result.lce);
            const int smem_end = smem_start + smem_length - 1;

            if (smem_length >= L_int)
            {
                std::vector<uint16_t> original_indices = getOriginalIndicesForInterval(
                    smem_end + 1,
                    lcp_result.interval.first,
                    lcp_result.interval.second,
                    lcp_result.interval_end_hap);

                for (size_t idx = 0; idx < original_indices.size(); ++idx)
                {
                    uint16_t original_row = original_indices[idx];
                    int len = std::snprintf(buffer, sizeof(buffer), "%zu\t%u\t%d\t%d\t%d\n",
                                            query_index, original_row, smem_start, smem_end, smem_length);
                    if (len > 0)
                        output.append(buffer, static_cast<size_t>(len));
                }
            }

            left = smem_end - L_int + 2;
        }
    }
};

int main(int argc, char *argv[])
{
    if (argc < 4 || argc > 5)
    {
        std::cerr << "Usage: " << argv[0] << " <panel.bcf> <query.bcf> <L> [output_file]\n";
        return 1;
    }

    const std::string panel_file = argv[1];
    const std::string query_file = argv[2];

    int L = 0;
    try
    {
        L = std::stoi(argv[3]);
    }
    catch (...)
    {
        std::cerr << "Invalid L value: " << argv[3] << "\n";
        return 1;
    }
    if (L <= 0)
    {
        std::cerr << "L must be positive.\n";
        return 1;
    }

    std::string output_file = "smems.tsv";
    if (argc > 4)
        output_file = argv[4];

    try
    {
        PBWTColumn pbwt(panel_file, query_file, L, output_file);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}