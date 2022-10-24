# Installation and Dependencies

## Note
This toolbox is self-contained, therefore the solver Nek5000 and needed libraries as ADIOS2 are shipped in this repository. We note that the procedure to use the toolbox is more involved that simply activating a runtime parameter. This is done to take into consideration the Nek5000 philosophy of modifing the source as little as possible and running user specified functions from a so called ```user file```. In the course of these instructions we will specify what commands need to be added to this file for everything to run smoothly. We expect that any Nek5000 user will be able to follow these instructions without much consideration. If dificulties arise, do not doubt to contact us.

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

## Installation and Usage

## Compiling Nek5000
Nek5000 is a FORTRAN77 solver that must be compiled for each case to be run, therefore the tools must be compiled as well. The third party libraries, like ADIOS2, do not need to be compiled everytime changes happen in the code, therefore we recomend that this is done only once.

Just as for any other used-defined function in Nek5000, 3 main files must be modified to be able to use the toolbox:
1. [```compile_script```](./compile_script.md)
2. [```makefile_usr.inc```](./makefile_usr_inc.md)
3. [```<casename>.usr```](./casename_usr.md)

Templates of the modified files can be found at the [DCTB](https://github.com/KTH-Nek5000/NekDCTB/tree/main/DCTB) and [examples](https://github.com/KTH-Nek5000/NekDCTB/blob/main/examples/turbPipe/compile/) folders, however, a list of modifications to the files can be found by using the hyperlinks.

Having modified all the files accordingly, compilation follows by ensuring that the ```SIZE``` file is present in the folder and executing the command:
```sh
./compile_script --all
```

If ```ADIOS2``` has not been compiled, the following message will prompt and ADIOS2 will be compiled:
```sh
Message from adios2
```

In this case, after the compilation of ```ADIOS2```, Nek5000 will not automatically restart the compilation of the case, therefore after all the 3rd party libraries have been compiled, the user must run ```./compile_script --all``` once again to compile the case. If the process is succesful, the executable  ```Nek5000``` and the folder  ```config``` will be added to the  ```compile``` sub-directory.


### Linking Nek5000 with ADIOS2
Nek5000 is a FORTRAN77 solver that must be compiled for each case to be run, therefore the tools must be compiled as well. The third party libraries, like ADIOS2, do not need to be compiled everytime changes happen in the code, therefore we recomend that this is done only once. In this repository we include a [Nek5000](https://github.com/KTH-Nek5000/NekDCTB/tree/main/Nek5000) version that has been modified for linkage with the ADIOS2 library. For the sake of convenience, we have included a tested version of ADIOS2 under the [/Nek5000/3rd_party/adios2](https://github.com/KTH-Nek5000/NekDCTB/tree/main/Nek5000/3rd_party/adios2) that can be easily installed by executing the install.sh script in the folder. We will list the changes done to Nek5000 in later sections to facilitate the integration with users that want to use their own versions of the solver. 
