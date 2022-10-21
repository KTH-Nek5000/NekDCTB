# Installation and Dependencies

## Dependencies

1. The following dependencies are required for Nek5000
  * FORTRAN 77 compilers
  * C and C++ compilers
  * MPI 

2. The following dependencies are required for ADIOS2
  * CMake

In a super computer such as [Dardel](https://www.pdc.kth.se/hpc-services/computing-systems/about-dardel-1.1053338) at PDC in Stockholm, loading the following modules is enough for compilation and execution:
```sh
module swap PrgEnv-cray PrgEnv-gnu
module load PDC
module load CMake/3.21.2
module swap craype-network-ofi craype-network-ucx
module swap cray-mpich cray-ucx
module load cray-mpich-ucx
```
