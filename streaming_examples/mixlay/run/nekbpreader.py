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

from mpi4py import MPI
import numpy as np
import matplotlib.pyplot as plt
import adios2

import sys
sys.path.append('./modules/')
from nek_snaps import dbCreator,dbReader,massMatrixReader
from plotters import contour2d,scatter

# MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# ADIOS portion
adios = adios2.ADIOS(comm)

if rank == 0:
    # ADIOS IO
    ioRead = adios.DeclareIO("ioReader")
    ioRead.SetEngine('BPFile')
    #ioRead.SetEngine('InSituMPI')

    ibpStream = ioRead.Open('geo.bp', adios2.Mode.Read, MPI.COMM_SELF)

    stepStatus = ibpStream.BeginStep()

    if stepStatus == adios2.StepStatus.OK:

        currentStep = ibpStream.CurrentStep() 
        var_inVX = ioRead.InquireVariable("X")
        var_inVY = ioRead.InquireVariable("Y")

        if var_inVX is not None:
            print('THE SHAPE OF X IS')
            print(var_inVX.Shape()) 
            #Set up array sizes in the first step
            if currentStep == 0:
                total_size=var_inVX.Shape()[0]
                my_count= int(total_size)
                my_start= int(0)
                #if total_size % size != 0:
                #    if rank < (total_size % size):
                #        my_count += 1
                #        my_start += rank
                #    else:
                #        my_start += (total_size%size)
                        
                #Allocate the vector
                inVX = np.zeros((my_count,1), dtype=np.double)
                inVY = np.zeros((my_count,1), dtype=np.double)

            #Select the data in the global array that belongs to me
            var_inVX.SetSelection([[my_start], [my_count]])
            var_inVY.SetSelection([[my_start], [my_count]])

            #Read the variable into array
            ibpStream.Get(var_inVX, inVX)
            ibpStream.Get(var_inVY, inVY)

            ibpStream.EndStep()

    ibpStream.Close()

    numelx=80
    numely=20
    lx1=8
    n=numelx*numely*lx1*lx1
    m=500
    k=30
    iSnap=120
    iMode=1

    x=inVX
    y=inVY

    print('Reading the global element number')
    lglel=np.loadtxt("LGLEL.txt")
    print('Reading the snapshot matrix')
    Ux=np.loadtxt("VX.txt").reshape(n,m)
    print('Reading the scaled snapshot matrix')
    utot_w=np.loadtxt("scaledVX.txt").reshape(n,m)
    print('Reading the POD modes')
    U=np.loadtxt("U.txt").reshape(n,k)
    print('Reading the mass matrix')
    bm1sqrt=np.loadtxt("bm1sqrt.txt")
    print('Reading D matrix')
    D=np.loadtxt("D.txt")
    print('Reading Vt matrix')
    Vt=np.loadtxt("Vt.txt").reshape(k,m)
    snapshot=Ux[:,iSnap]
    mode=U[:,iMode]
    
    print('Calculate time coefficients')
    T=np.diag(D)@Vt

    print('For streaming:')
    print('The shape of U and T are:')
    print(U.shape)
    print(T.shape)

    print('Reconstructing one snapshot from POD')
    ux=U@T

    snaprec=ux[:,iSnap]

    #ind=np.where(lglel==1)[0][0]
    #ind2=np.where(lglel==81)[0][0]
    #print(x[(ind)*lx1*lx1:(ind+1)*lx1*lx1])
    #print(y[(ind)*lx1*lx1:(ind+1)*lx1*lx1])
    #print(x[(ind2)*lx1*lx1:(ind2+1)*lx1*lx1])
    #print(y[(ind2)*lx1*lx1:(ind2+1)*lx1*lx1])


    print('Calculating the SVD of the snapshot matrix')
    Ur2,Dr2,Vtr2=np.linalg.svd(utot_w, full_matrices=False)

    mode_standar=Ur2[:,iMode]

    print('Rconstructing from stardan SVD')
    uxs=Ur2[:,0:k]@np.diag(Dr2[0:k])@Vtr2[0:k,:]

    print('For standar:')
    print('The shape of U and T are:')
    print(Ur2.shape)
    print((np.diag(Dr2)@Vtr2).shape)

    snaprecs=uxs[:,iSnap]

    for j in range(0,n):
        mode_standar[j]= mode_standar[j]/bm1sqrt[j]
        snaprecs[j]= snaprecs[j]/bm1sqrt[j]




    print('Reorganizing the data for plotting')
    X=np.zeros((numely*lx1,numelx*lx1))
    Y=np.zeros((numely*lx1,numelx*lx1))
    UU=np.zeros((numely*lx1,numelx*lx1))
    Urec=np.zeros((numely*lx1,numelx*lx1))
    Urecs=np.zeros((numely*lx1,numelx*lx1))
    Umod=np.zeros((numely*lx1,numelx*lx1))
    Umods=np.zeros((numely*lx1,numelx*lx1))
    kk=1
    for i in range(0,numely):
        for j in range(0,numelx):
            iwantelem=kk
            kk=kk+1
            elind=np.where(lglel==iwantelem)[0][0]
            elemdatax=(x[(elind)*lx1*lx1:(elind+1)*lx1*lx1])
            elemdatay=(y[(elind)*lx1*lx1:(elind+1)*lx1*lx1])
            elemdatasnap=(snapshot[(elind)*lx1*lx1:(elind+1)*lx1*lx1])
            elemdatamode=(mode[(elind)*lx1*lx1:(elind+1)*lx1*lx1])
            elemdatamodes=(mode_standar[(elind)*lx1*lx1:(elind+1)*lx1*lx1])
            elemdatasnaprec=(snaprec[(elind)*lx1*lx1:(elind+1)*lx1*lx1])
            elemdatasnaprecs=(snaprecs[(elind)*lx1*lx1:(elind+1)*lx1*lx1])

            for i2 in range(0,lx1):
                for j2 in range(0,lx1):
                    X[i*lx1+i2,j*lx1+j2]=elemdatax[i2*lx1+j2]
                    Y[i*lx1+i2,j*lx1+j2]=elemdatay[i2*lx1+j2]
                    UU[i*lx1+i2,j*lx1+j2]=elemdatasnap[i2*lx1+j2]
                    Umod[i*lx1+i2,j*lx1+j2]=elemdatamode[i2*lx1+j2]
                    Umods[i*lx1+i2,j*lx1+j2]=elemdatamodes[i2*lx1+j2]
                    Urec[i*lx1+i2,j*lx1+j2]=elemdatasnaprec[i2*lx1+j2]
                    Urecs[i*lx1+i2,j*lx1+j2]=elemdatasnaprecs[i2*lx1+j2]
                    

    print('Plotting results')

    contour2d(UU.T,X,Y,title='snapshot %d' %iSnap)
    contour2d(Urec.T,X,Y,title='Streaming - Reconstructed snapshot' +repr(iSnap)+', from truncating to '+repr(k)+' Modes ')
    contour2d(Urecs.T,X,Y,title='Standard - Reconstructed snapshot' +repr(iSnap)+', from truncating to '+repr(k)+' Modes ')

    #results from streaming
    A=(Urec-Urecs)
    plt.figure(figsize=(7,5))
    ax = plt.subplot(1,1,1)
    p = ax.pcolormesh(A)
    plt.colorbar(p)
    plt.show()

    contour2d(Umod.T,X,Y,title='Mode %d from streaming algorithm' %iMode)
    contour2d(Umods.T,X,Y,title='Mode %d from standard algorithm' %iMode)

    #results from streaming
    A=np.abs((np.abs(Umod)-np.abs(Umods)))
    plt.figure(figsize=(7,5))
    ax = plt.subplot(1,1,1)
    p = ax.pcolormesh(A)
    plt.colorbar(p)
    plt.show()


    #results from streamingi
    print(D)
    print(Dr2[0:k])
    #print(np.diag(Dr2[0:k]))
    A=np.abs(np.diag(D)-np.diag(Dr2[0:k]))
    print(A)
    plt.figure(figsize=(7,5))
    ax = plt.subplot(1,1,1)
    p = ax.pcolormesh(A)
    plt.colorbar(p)
    plt.show()


    #results from streaming
    #print(Vt)
    #print(Vtr2[0:k,:])
    A=np.abs(Vt-Vtr2[0:k,:])
    plt.figure(figsize=(7,5))
    ax = plt.subplot(1,1,1)
    p = ax.pcolormesh(A)
    plt.colorbar(p)
    plt.show()
