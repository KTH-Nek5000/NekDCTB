#!/bin/bash
rm d* e* k* f* hat* crt* rct* *.nek5000
#rm -r *.bp
rm nek5000
rm -r config
cp ../compile_compression/nek5000 .
cp -r ../compile_compression/config .
cp ../compile_compression/compressor
echo turbPipe > SESSION.NAME
echo `pwd`'/' >> SESSION.NAME
cp data/* .
mpirun -np 4 ./nek5000 : -np 2 ./compressor

