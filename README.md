# PBML
PBML is adaptation of BML (Boyer-Moore-Li) trick into indexing and querying of PBWT (Positional Burrows--Wheeler Transform). Boyer-Moore string search algorithm and Heng Li's forward backward was combined into one algorithm named BML by T. Gagie. PBML adapts BML into the PBWT framework and performs and contiguous memory allocation trick. PBML builds two PBWTs one in the forward direction another in the reverse direction to perform LCP and LCS and a run-based haplotype-retrieval and can report SMEMs (Set Maximal Exact Match) of configurable minimum length given a panel and a set of queries (.bcf) format, it uses SDSL libraries. It can outperform available run-sampled PBWT tools in querying and indexing time with competitive memory usage. PBML is highly scalable and multi-thread querying performance achieves signficant querying speedup with near constant memory footprint.

## Build
clone the repository and inside the repository directory run the following commands.
```
#Requires CMake installed, once CMake is installed run:

mkdir build
cd build
cmake ..
make -j
```
The binary will be stored in the build directory named `pbml`.


## Usage
```
#mono-thread
OMP_NUM_THREADS=1 ./build/pbml <dir/.../panel.bcf> <dir/.../query.bcf> <minimum_SMEM_length> <output_name>

#multi-thread (n threads)
OMP_NUM_THREADS=n ./build/pbml <dir/.../panel.bcf> <dir/.../query.bcf> <minimum_SMEM_length> <output_name>

```

## Output format
```
#Each line reports a SMEM with following tab seperated information:
<haplotype> <starting column> <ending column> <length>
```
