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

3 main files must be modified to be able to use the toolbox.
1. ```compile_script```
2. ```makefile_usr.inc```
3. ```<casename>.usr```


### compile_script
We recomend to use a compile script as the one shown in the [examples](https://github.com/KTH-Nek5000/NekDCTB/blob/main/examples/turbPipe/compile/compile_script) folder. The structure is the same as would be needed to compile Nek5000 but has the following aditions:
  1. The path to the data compression tool box must be specified.
  ```sh
  export COMPRESS_SRC="$PWD/../../../DCTB"
  export COMPRESS_INC=" -I"${COMPRESS_SRC}
  COMPRESS_INC+=" -I"${COMPRESS_SRC}"/framework_files"
  ```
  2. The folders with the subroutines must be added to the shearch space of the fortran compiler.
  ```sh
  export FC="mpifort"${COMPRESS_INC}
  ```
  3. The following user defined functions must included in the compilation of Nek5000. Users of the [KTH-Framework](https://github.com/KTH-Nek5000/KTH_Framework) will notice that the first row correspond to supporting files. The data compression toolbox itself is contained in ```io_trunc.o```
  ```sh
  export USR="frame.o mntrlog_block.o mntrlog.o mntrtmr_block.o mntrtmr.o rprm_block.o rprm.o io_tools_block.o io_tools.o"
  USR+=" io_trunc.o sperri_trunc.o"
  ```
  4. Finally, some files from the toolbox folder must be copied into the compile script. These are files that will later be used by ADIOS2 for runtime confugiration of  compression.
  ```sh
  cp -r ${COMPRESS_SRC}/config .
  cp -r config/synch_config.xml config/config.xml
  cp ${NEK_SOURCE_ROOT}/core/3rd_party/opt/single_synch_nek_adios2.cpp ${NEK_SOURCE_ROOT}/core/3rd_party/nek_adios2.cpp
  ```
### makefile_usr.inc
The ```makefile_usr.inc``` must be modified to include the user functions. This for the toolbox to work, we must include the following in the file:
```sh
$(OBJDIR)/frame.o           :${COMPRESS_SRC}/framework_files/frame.f;              $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/mntrlog_block.o   :${COMPRESS_SRC}/framework_files/mntrlog_block.f;      $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/mntrlog.o         :${COMPRESS_SRC}/framework_files/mntrlog.f;            $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/mntrtmr_block.o   :${COMPRESS_SRC}/framework_files/mntrtmr_block.f;      $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/mntrtmr.o         :${COMPRESS_SRC}/framework_files/mntrtmr.f;            $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/rprm_block.o      :${COMPRESS_SRC}/framework_files/rprm_block.f;         $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/rprm.o            :${COMPRESS_SRC}/framework_files/rprm.f;               $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/io_tools_block.o  :${COMPRESS_SRC}/framework_files/io_tools_block.f;     $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/io_tools.o        :${COMPRESS_SRC}/framework_files/io_tools.f;           $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/sperri_trunc.o        :${COMPRESS_SRC}/framework_files/sperri_trunc.f;           $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/io_trunc.o        :${COMPRESS_SRC}/io_trunc.f;                           $(F77) -c $(FL2) -I./ $< -o $@
```


### Linking Nek5000 with ADIOS2
Nek5000 is a FORTRAN77 solver that must be compiled for each case to be run, therefore the tools must be compiled as well. The third party libraries, like ADIOS2, do not need to be compiled everytime changes happen in the code, therefore we recomend that this is done only once. In this repository we include a [Nek5000](https://github.com/KTH-Nek5000/NekDCTB/tree/main/Nek5000) version that has been modified for linkage with the ADIOS2 library. For the sake of convenience, we have included a tested version of ADIOS2 under the [/Nek5000/3rd_party/adios2](https://github.com/KTH-Nek5000/NekDCTB/tree/main/Nek5000/3rd_party/adios2) that can be easily installed by executing the install.sh script in the folder. We will list the changes done to Nek5000 in later sections to facilitate the integration with users that want to use their own versions of the solver. 
