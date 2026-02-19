# PBML

PBML (Positional Boyer-Moore-Li) adapts the BML (Boyer-Moore-Li) algorithm — combining Boyer-Moore string search with Heng Li's forward-backward algorithm, originally proposed by T. Gagie for BWT — into the PBWT (Positional Burrows-Wheeler Transform) framework. It builds forward and reverse RLE-compressed PBWTs to perform LCP/LCS queries and reports k-SMEMs (Set Maximal Exact Matches) of configurable minimum length from query haplotypes against a reference panel.

PBML uses SDSL and htslib libraries. Input files must be in `.bcf/.vcf` format.

## Implementations

Two implementations are provided in `src/`:

- **`pbml.cpp`** — Uses φ (phi) operations with successor arrays and hint-guided lookups (inspired by MOVI) to recover haplotype IDs directly during parallel querying. Larger index (~20R + 12N bytes) but no reconstruction pass needed.

- **`pbmlRecon.cpp`** — Two-phase approach: Phase 1 runs BML in parallel collecting SMEM descriptors (intervals only), Phase 2 replays the prefix array left-to-right to resolve haplotype IDs. Smaller index (~2(R+R_rev) + 12N bytes) with no phi/successor structures.

Both produce identical output.

## Build

```bash
git clone https://github.com/uwaiseibna/PBML.git
cd PBML
mkdir build && cd build
cmake .. && make -j
```

This builds two binaries in `build/`: `pbml` and `pbmlRecon`.

## Usage

```
./pbml <mode> [options]

Modes:
  index    Build and save index from panel
  query    Query from saved index
  run      Build and query (no index saved)

Options:
  -p, --panel <file>     Panel BCF/VCF file (index/run modes)
  -q, --query <file>     Query BCF/VCF file (query/run modes)
  -i, --index <file>     Index file (query: input, index: output)
  -o, --output <file>    Output file [default: smems.tsv]
  -L, --length <int>     Minimum SMEM length [default: 1]
  -k, --min-occ <int>    Minimum occurrences [default: 1]
  -t, --threads <int>    Number of threads [default: all]
```

### Examples

```bash
# Build index
./pbml index -p panel.bcf -i panel.pbml

# Query from index
./pbml query -i panel.pbml -q query.bcf -L 5 -o smems.tsv

# Build and query in one step
./pbml run -p panel.bcf -q query.bcf -L 5 -k 1

# Single-threaded
OMP_NUM_THREADS=1 ./pbml run -p panel.bcf -q query.bcf -L 5
```

Same usage applies to `pbmlRecon`.

## Output Format

Tab-separated, one line per SMEM match:

```
query_id    panel_haplotype    start_site    end_site    length
```

## Citation

```
@article{islam2025scalable,
  title={Scalable PBWT Queries with Minimum-Length SMEM Constraints},
  author={Islam, Uwaise Ibna and Cozzi, Davide and Gagie, Travis and Varki, Rahul and Colonna, Vincenza and Garrison, Erik and Bonizzoni, Paola and Boucher, Christina},
  journal={bioRxiv},
  pages={2025--12},
  year={2025},
  publisher={Cold Spring Harbor Laboratory}
}
```
