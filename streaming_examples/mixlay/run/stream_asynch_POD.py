import sys
import os
import copy
os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "1" # export OPENBLAS_NUM_THREADS=4 
os.environ["MKL_NUM_THREADS"] = "1" # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "1" # export NUMEXPR_NUM_THREADS=6

from mpi4py import MPI #equivalent to the use of MPI_init() in C
import numpy as np
import scipy.optimize
import time
import matplotlib.pyplot as plt
import pymech as pm
from pymech.neksuite import readnek
from pymech.neksuite import writenek
from tqdm import tqdm
import mpi_spSVD as spSVD
import pandas as pd

# Import adios2
sys.path.append('/home/adalberto/NekDCTB/Nek5000/3rd_party/adios2/lib/python3.8/site-packages/adios2')
import adios2

#####################################
# Case parameters
## Number of snapshots
nsnap=1000
## Batch size
p=5 # Number of snapshots in a batch
## Pod of which quantity
qoi=0 #0=vx, 1=vy, 2=vz

# Data path
casename='mixlay'
path='./snapshots/'

# number of modes
k=10
maxk = 1000
# Options if autodecide number of modes
## Autoupdate?
ifautoupdate=False
## If autoupdate has been done, gather ALL modes?
ifget_all_modes=False
## Minimun percentage of new snapshots that need to be parallel
min_ratio=0.01

# Update gbl or lcl?
ifgbl_update=True
num_recs = 10

#### Automatic, do not touch #####
if ifgbl_update==True: 
    outprefix='gbl'
else:
    outprefix='lcl'


if __name__ == "__main__":

    # MPI - MPMD
    worldcomm = MPI.COMM_WORLD
    worldrank = worldcomm.Get_rank()
    worldsize = worldcomm.Get_size()
    col = 1
    comm = worldcomm.Split(col,worldrank)
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    # ADIOS portion
    adios = adios2.ADIOS(comm)

    # ADIOS IO - Engine
    ioRead = adios.DeclareIO("ioReader")
    ioRead.SetEngine('InSituMPI')

    if rank ==0: print(outprefix)
    
    #define the class
    spPOD = spSVD.spSVD(k=k, p=p, nsnap=nsnap, qoi= qoi, path =path, casename=casename, ifautoupdate=ifautoupdate, min_ratio=min_ratio, ifget_all_modes=ifget_all_modes,ifgbl_update=ifgbl_update)
    #retreive case information
    spPOD.get_case_info(comm)
    nel = spPOD.nel
    lx1 = spPOD.lx1
    ly1 = spPOD.ly1
    lz1 = spPOD.lz1
    nxyz= spPOD.nxyz
    nxyze= spPOD.nxyze

    print('#############')
    print(nel)
    print(lx1)
    print(ly1)
    print(lz1)
    print(nxyz)
    print(nxyze)

    print('#############')

    #Create buffer for reading data and redistributing it
    if (rank==0): 
        v=np.empty((nxyze,1))
        bm1sqrt=np.empty((nxyze,1))
        lglel=np.empty((nel,1),dtype=np.int32)
    else:
        v=None
        bm1sqrt=None
        lglel=None

    # Create dummies for the streaming POD
    U_1t = None
    D_1t = None
    Vt_1t = None

    # Change k to a low value if autoupdate is required
    if ifautoupdate==True: spPOD.k=p

    running_ra=[]

    #This is the streaming process
    # Open the insitu array... It is open until stream ends
    ibpStream = ioRead.Open("globalArray", adios2.Mode.Read, comm)
    
    isnap = -1
    while True:
 
        stepStatus = ibpStream.BeginStep()

        if stepStatus == adios2.StepStatus.OK:

            isnap = int(isnap + 1)
            currentStep = ibpStream.CurrentStep()

            # Get the data
            if (currentStep == 0) :
                Xi, Xtemp, buff, bm1i, bm1sqrti, lgleli = spPOD.io_allocate_adios(ioRead, p, comm)
 
            Xi, bm1sqrti, lgleli = spPOD.io_read_adios(ioRead, ibpStream, Xi, bm1i, bm1sqrti, lgleli, comm)
            ibpStream.EndStep()


            print("scale data")
            #Scale the data with the mass matrix
            Xi=spPOD.scale_data(Xi,bm1sqrti,int(nxyze/size),1,'mult')

            # Update the number of snapshots processed up to now
            m=isnap+1

            print("fill buffer")
            # fill in the buffer
            buff[:,np.mod(isnap,spPOD.p)] = Xi[:,0]
            cbf = np.mod(isnap,spPOD.p)


            # Calculate the residual and check if basis needs to be expanded 
            if isnap>=spPOD.p: 
                if ifautoupdate==True:
                    if ifgbl_update == False:
                        ra=spPOD.get_perp_ratio(U_1t,Xi.reshape((-1,1)))
                        running_ra.append(ra)
                    else:
                        ra=spPOD.mpi_get_perp_ratio(U_1t,Xi.reshape((-1,1)),comm)
                        running_ra.append(ra)
                else:
                    ra=0
                    running_ra.append(ra)
                if spPOD.ifautoupdate==True and ra>=spPOD.min_ratio and spPOD.k<maxk: spPOD.k+=1

            # Update
            if np.mod(isnap+1,spPOD.p) == 0 or (isnap+1)==nsnap:
                # Print the local state
                print("rank "+repr(rank)+" has k="+repr(spPOD.k))
                # Update the svd with each new snapshot 
                if ifgbl_update==True:
                    U_1t,D_1t,Vt_1t = spPOD.gblSVD_update_fromBatch(U_1t,D_1t,Vt_1t,buff[:,:(cbf+1)],isnap,comm)
                else:
                    U_1t,D_1t,Vt_1t = spPOD.lclSVD_update_fromBatch(U_1t,D_1t,Vt_1t,buff[:,:(cbf+1)],isnap)
               
                ##Temporal test
                #print("temporal test")
                #temp=U_1t@np.diag(D_1t)@Vt_1t
                #print(np.max(Xi[:,0]-temp[:,-1]))
    
            #pbar.update(1)
        #pbar.close()
        elif stepStatus == adios2.StepStatus.EndOfStream:
            print("End of stream")
            break

    ibpStream.Close()
    print("closed stream")

    #Do this in case you want to keep more data than you actually collected
    if spPOD.setk>=m: spPOD.setk=m

    print("scaling next data")
    ##Scale the modes back before gathering them
    U_1t=spPOD.scale_data(U_1t,bm1sqrti,int(nxyze/size),spPOD.k,'div')

    ## If local update, reconstruct for later comparison
    if ifgbl_update == False: Si = U_1t@np.diag(D_1t)@Vt_1t[:,0:num_recs]

    # Obtain a global svd from a local one only if the svd has been updated locally
    if ifgbl_update == False: U_1t,D_1t,Vt_1t = spPOD.lclSVD_to_gblSVD(U_1t,D_1t,Vt_1t,comm)
    
    if ifgbl_update == False: 
        ortho = comm.gather(running_ra,root=0)
    else:
        ortho = running_ra

    print(U_1t.shape)

    #Once the streaming is done. Gather the modes at one rank
    if ifget_all_modes==True: kk=m
    if ifget_all_modes==False:
        if ifgbl_update == True:
            if spPOD.setk <= spPOD.k: # This line here is to protect the code from truncating when self updating, other wise they are the same
                kk=spPOD.setk
            else:
                kk=spPOD.k
        else:
            kk=spPOD.setk

        U_1t=np.copy(U_1t[:,0:kk])
        D_1t=np.copy(D_1t[0:kk])
        Vt_1t=np.copy(Vt_1t[0:kk,:])

    U, bm1sqrt2 = spPOD.gatherModesandMass_atRank0(U_1t,bm1sqrti,spPOD.nxyze,kk,comm)

    comm.Gather([lgleli,MPI.INT], [lglel, MPI.INT], root=0)


    if ifgbl_update == False: S, _ = spPOD.gatherModesandMass_atRank0(Si,bm1sqrti,spPOD.nxyze,num_recs,comm)

    #Write the modes
    if rank==0:
        
        # Save the running orthogonality and others
        np.save(spPOD.path+'PODrra'+'_'+outprefix+'_'+spPOD.casename,np.array(ortho))
    
        print(U.shape)
        print(D_1t.shape)
        print(Vt_1t.shape)

        ## Do this because when globaly autoupdating, you are constrained
        #if ifgbl_update== True and ifautoupdate== True:
        #    kkk=spPOD.k
        #else:
        #    kkk=spPOD.setk

        spPOD.write_modes(U,lglel,0,kk,outprefix)
        spPOD.write_timecoeff(D_1t,Vt_1t,outprefix)
        
        if ifgbl_update == False:
            spPOD.write_rec(S,0,num_recs,outprefix)
