#!/bin/bash
rm d* e* k* f* hat* crt* rct* *.nek5000
#rm -r *.bp
rm nek5000
rm compressor
rm -r config
cp ../compile/nek5000 .
cp ../compile/compressor .
cp -r ../compile/config .
export PATH=/home/adalberto/Nek5000_adios2/bin:$PATH
echo mixlay > SESSION.NAME
echo `pwd`'/' >> SESSION.NAME
mpirun -np 12 ./nek5000 : -np 2 ./compressor

