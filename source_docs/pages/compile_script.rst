.. _compile_script:

compile_script
==============

We recomend to use a compile script as the one shown in the
`examples <https://github.com/KTH-Nek5000/NekDCTB/blob/main/examples/turbPipe/compile/compile_script>`__
folder. The structure is the same as would be needed to compile Nek5000
but has the following aditions: 

1. The path to the data compression tool box must be specified.

.. code:: sh

   export COMPRESS_SRC="$PWD/../../../DCTB"
   export COMPRESS_INC=" -I"${COMPRESS_SRC}   
   COMPRESS_INC+=" -I"${COMPRESS_SRC}"/framework_files"

2. The folders with the subroutines must be added to the shearch space
of the fortran compiler. 

.. code:: sh

   export FC="mpifort"${COMPRESS_INC} 

3. ADIOS2 instrumentation must be activated. This is done by including the
following flag in the preprocesor list. 

.. code:: sh

    export PPLIST="ADIOS2"

4. The following user defined functions must included in the compilation
   of Nek5000. Users of the
   `KTH-Framework <https://github.com/KTH-Nek5000/KTH_Framework>`__ will
   notice that the first row correspond to supporting files. The data
   compression toolbox itself is contained in ``io_trunc.o``

.. code:: sh

   export USR="frame.o mntrlog_block.o mntrlog.o mntrtmr_block.o mntrtmr.o rprm_block.o rprm.o io_tools_block.o io_tools.o"
   USR+=" io_trunc.o sperri_trunc.o"

5. Finally, some files from the toolbox folder must be copied into the
   compile script. These are files that will later be used by ADIOS2 for
   runtime configuration of the lossless compression.

.. code:: sh

   cp -r ${COMPRESS_SRC}/config .
   cp -r config/synch_config.xml config/config.xml
   cp ${NEK_SOURCE_ROOT}/core/3rd_party/opt/single_synch_nek_adios2.cpp ${NEK_SOURCE_ROOT}/core/3rd_party/nek_adios2.cpp
