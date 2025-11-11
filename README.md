# PBML

## Build

```
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
