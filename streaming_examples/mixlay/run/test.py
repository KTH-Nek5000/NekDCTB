#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# helloBPReaderHeatMap2D.py
#
#
#  Created on: Dec 5th, 2017
#      Author: William F Godoy godoywf@ornl.gov
#
import sys
sys.path.append('/home/adalberto/NekDCTB/Nek5000/3rd_party/adios2/lib/python3.8/site-packages/adios2')


import sys
sys.path.append('./modules/')
from nek_snaps import dbCreator,dbReader,massMatrixReader
from plotters import contour2d,scatter

import numpy as np

numelx=80
numely=20
lx1=8
n=numelx*numely*lx1*lx1
m=500
k=30
iSnap=120
iMode=15
    
print('Reading the POD modes')
U=np.loadtxt("U.txt").reshape(n,k)
Usnap=np.loadtxt("VX.txt").reshape(n,m)
#print('Reading the mass matrix')
#bm1sqrt=np.loadtxt("bm1sqrt.txt")
#print('Reading D matrix')
#D=np.loadtxt("D.txt")
print('Reading Vt matrix')
Vt=np.loadtxt("Vt.txt").reshape(k,m)


print((U.T)@U)


#print('Calculate time coefficients')
#T=np.diag(D)@Vt

#print('Reconstructing one snapshot from POD')

