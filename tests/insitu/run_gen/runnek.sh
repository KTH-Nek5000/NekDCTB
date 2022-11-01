#!/bin/bash
rm d* e* k* f* hat* crt* rct* *.nek5000
#rm -r *.bp
rm nek5000
rm -r config
cp ../compile/nek5000 .
cp -r ../compile/config .
echo turbPipe > SESSION.NAME
echo `pwd`'/' >> SESSION.NAME
mpirun -np 4 ./nek5000

