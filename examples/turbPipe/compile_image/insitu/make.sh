#!/bin/bash

export OSMESA=/path/to/OSMESA

export LIBDIR=$OSMESA/lib:$LIBDIR
export LD_LIBRARY_PATH=$OSMESA/lib:$LD_LIBRARY_PATH
export LD_RUN_PATH=$OSMESA/lib:$LD_RUN_PATH

export OSMESA_INCLUDE_DIR=$OSMESA/include
export OSMESA_LIBRARY=$OSMESA/lib

export PARAVIEW=/path/to/paraview
export PATH=$PARAVIEW/bin:$PATH
export LD_LIBRARY_PATH=$PARAVIEW/lib:$LD_LIBRARY_PATH

export PYTHONPATH=$PARAVIEW/lib/python3.8/site-packages:$PYTHONPATH
export PYTHONPATH=$PARAVIEW/lib/python3.8/site-packages/paraview/:$PYTHONPATH

export ADIOS2_DIR=/path/to/adios
export PATH=$ADIOS2_DIR/bin:$PATH

export ADIOS2=1
export ADIOS2_INCS=`adios2-config --cxx-flags`
export ADIOS2_LIBS=`adios2-config --cxx-libs`

export CATALYST=1
export CATALYST_LIBS+=`paraview-config -l -c Catalyst`
export CATALYST_INCS=`paraview-config -f -c Catalyst`

export CXXFLAG=" -g -O2 "

mpicxx -c $CXXFLAG $CATALYST_INCS myCPPythonAdaptorAPI.cpp -o myCPPythonAdaptorAPI.o
mpicxx -c $CXXFLAG $CATALYST_INCS nek_catalyst.cpp -o nek_catalyst.o
mpicxx $CXXFLAG $CATALYST_INCS $ADIOS2_INCS adios_catalyst.cpp nek_catalyst.o myCPPythonAdaptorAPI.o $ADIOS2_LIBS $CATALYST_LIBS -o catalystSpace
