# PBML

**Fast, memory-efficient enumeration of set-maximal exact matches (SMEMs) from query haplotypes against a reference panel.**

PBML finds all haplotype segments in a query that match at least *k* panel haplotypes and cannot be extended without losing matches (SMEMs) — with a configurable minimum length *L*. No other PBWT-based tool supports minimum-length SMEM constraints (*kL*-SMEMs).

## Highlights

- **kL-SMEMs** — Find SMEMs of length ≥ *L* with at least *k* occurrences. Unique to PBML; competing tools are limited to L=1.
- **4.6× faster** than μ-PBWT, **2.4× faster** than Durbin's PBWT, using **1.3–25.6× less memory** (1KGP, single-threaded).
- **O(r) index size** — Run-length compressed PBWTs scale with panel compressibility, not raw size.
- **Parallel queries** — 8.2× speedup at 16 threads with near-constant memory overhead.

## Quick Start

```bash
git clone https://github.com/uwaiseibna/PBML.git
cd PBML && mkdir build && cd build
cmake .. && make -j
./pbml run -p panel.bcf -q query.bcf -L 5 -k 1 -o smems.tsv
```

Requires [htslib](https://github.com/samtools/htslib) and [SDSL](https://github.com/simongog/sdsl-lite). Input files must be `.bcf` or `.vcf`.

## Performance

Median across chromosomes 1–22 on 1000 Genomes Project Phase 3 (4,008 panel haplotypes, 1,000 queries, L=1, single-threaded). Build and query times normalized per million variant sites.

| Method | Build (s/M sites) | Query (s/M sites) | Peak Memory (GB) |
|--------|:-:|:-:|:-:|
| **PBML** | 45.19 (43.86–46.12) | **49.08** (44.86–67.42) | **6.4** (2.4–12.4) |
| [μ-PBWT](https://github.com/dlcgold/muPBWT) | 63.76 (62.72–65.88) | 233.90 (200.48–276.75) | 9.1 (2.8–17.4) |
| [Dynamic μ-PBWT](https://github.com/ucfcbb/Dynamic-mu-PBWT) | 212.45 (203.46–351.44) | 582.28 (485.38–756.82) | 24.2 (9.4–46.0) |
| [PBWT](https://github.com/richarddurbin/pbwt) | **20.20** (14.95–20.60) | 121.75 (113.36–135.03) | 176.2 (52.3–336.2) |

Values are median (min–max). All methods report identical SMEM sets. With 16 threads, PBML achieves an 8.2× query speedup with <1% memory increase, while μ-PBWT sees 15.9% memory growth for a 2.4× speedup.

## Usage

```
./pbml <mode> [options]
```

### Modes

| Mode | Description |
|------|-------------|
| `index` | Build and save index from panel |
| `query` | Query from a saved index |
| `run` | Build and query in one step (no index saved) |

### Options

| Flag | Description | Default |
|------|-------------|---------|
| `-p, --panel` | Panel BCF/VCF file (index/run) | — |
| `-q, --query` | Query BCF/VCF file (query/run) | — |
| `-i, --index` | Index file path | — |
| `-o, --output` | Output file | `smems.tsv` |
| `-L, --length` | Minimum SMEM length | `1` |
| `-k, --min-occ` | Minimum occurrences in panel | `1` |
| `-t, --threads` | Number of threads | all available |

### Examples

```bash
# Build index once, query many times
./pbml index -p panel.bcf -i panel.pbml
./pbml query -i panel.pbml -q query.bcf -L 10 -k 2 -o smems.tsv

# One-shot build + query
./pbml run -p panel.bcf -q query.bcf -L 5 -k 1

# Single-threaded
OMP_NUM_THREADS=1 ./pbml run -p panel.bcf -q query.bcf -L 5
```

## Output Format

Tab-separated, one line per (query, panel haplotype, SMEM) triple:

```
query_id    panel_haplotype    start_site    end_site    length
```

## Implementations

Two variants are provided in `src/`, both producing identical output:

| Binary | Haplotype recovery | Best for |
|--------|-------------------|----------|
| `pbml` | φ-based lookups during querying | Multi-threaded workloads, large panels (≥5K haplotypes) |
| `pbmlRecon` | Sequential prefix array replay after querying | Single-threaded on small panels, lowest memory |

`pbml` resolves haplotype IDs independently per query via constant-time φ operations, so all queries run fully in parallel. `pbmlRecon` eliminates the φ/successor structures for a smaller index and faster construction, but its sequential O(w×h) reconstruction pass limits parallel scalability and becomes a bottleneck as panel size grows.

## How It Works

PBML adapts the BML (Boyer-Moore-Li) algorithm — originally proposed by T. Gagie for BWT — to the PBWT framework. It builds forward and reverse run-length encoded PBWTs and uses LCP/LCS queries with Boyer-Moore skip logic to enumerate SMEMs efficiently. See the [paper](https://www.biorxiv.org/content/10.64898/2025.12.01.691644v1)for details.

## Citation

```bibtex
@article{islam2025scalable,
  title   = {Scalable PBWT Queries with Minimum-Length SMEM Constraints},
  author  = {Islam, Uwaise Ibna and Cozzi, Davide and Gagie, Travis and
             Varki, Rahul and Colonna, Vincenza and Garrison, Erik and
             Bonizzoni, Paola and Boucher, Christina},
  journal = {bioRxiv},
  pages   = {2025--12},
  year    = {2025},
  publisher = {Cold Spring Harbor Laboratory}
}
```
