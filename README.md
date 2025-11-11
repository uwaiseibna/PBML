# PBML

## Build
clone the repository and inside the repository directory run the following commands.
```
#requires CMake installed
#Once CMake is installed run:

mkdir build
cd build
cmake ..
make -j
```


## Usage
```
#mono-thread
OMP_NUM_THREADS=1 ./build/pbml <dir/.../panel.bcf> <dir/.../query.bcf> <minimum_SMEM_length> <output_name>

#multi-thread
OMP_NUM_THREADS=n ./build/pbml <dir/.../panel.bcf> <dir/.../query.bcf> <minimum_SMEM_length> <output_name>

```

## Output format
```
#Each line reports a SMEM with following tab seperated information:
<haplotype> <starting column> <ending column> <length>
```
