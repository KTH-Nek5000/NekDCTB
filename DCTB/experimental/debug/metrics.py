import matplotlib as mpl
import matplotlib.pyplot as plt
import pymech as pm
from pymech.neksuite import readnek
from pymech.neksuite import writenek
import copy
import numpy as np
import sys


#=================================================================

#Define some functions

def get_metrics(rctr_data,orgl_data,bm1_data,qoi):

    #Get some information from the data
    nel=orgl_data.nel
    lxr=orgl_data.lr1
    nxyz=np.prod(lxr) 

    #Allocate vectors for the metrics
    elem=np.zeros(nel)
    maxabserr_v=np.zeros(nel)
    maxrelerr_v=np.zeros(nel)
    rmse_v=np.zeros(nel)
    nrmse_v=np.zeros(nel)
    rho_v=np.zeros(nel)
    psnr_v=np.zeros(nel)
    qoi=qoi

    for e in range(0,nel):
        #Assign the element number
        elem[e]=e+1

        #Read the data in each element (ed=element data)
        rctr_ed=rctr_data.elem[e]  #for the reconstructed data
        bm1_ed=orgl_data.elem[e]  #for the original data

        #Rearange the velocities in an easier vector to read.
        orgl_v=orgl_data.elem[e].vel[qoi,:,:,:].reshape((nxyz,1),order='F')
        rctr_v=rctr_data.elem[e].vel[qoi,:,:,:].reshape((nxyz,1),order='F')
        bm1_v=bm1_data.elem[e].vel[qoi,:,:,:].reshape((nxyz,1),order='F')
        R_v=np.max(orgl_v)-np.min(orgl_v)

        #1. Calculate the first metric (Absolute error)
        abserr_v=orgl_v-rctr_v
        #get the maximun to print only one value
        maxabserr_v[e]=np.max(np.abs(abserr_v))

        #2. Calculate the second metric (Relative error)
        relerr_v=abserr_v/R_v
        #get the maximun to print only one value
        maxrelerr_v[e]=np.max(np.abs(relerr_v))

        #3. Calculate the third metric (RMSE)
        ##3.1 - weighted
        temp=0
        tempvol=0
        for i in range(0,nxyz):
            temp=temp+abserr_v[i,0]*abserr_v[i,0]*bm1_v[i,0]
            tempvol=tempvol+bm1_v[i,0]
        rmse_v[e]=np.sqrt(temp)/np.sqrt(tempvol)

        ##3.2 - not weighted
        #temp=np.matmul(abserr_v.T,abserr_v)/nxyz
        #rmse_v[e]=np.sqrt(temp[0,0])

        # Also calculate the normalized one
        nrmse_v[e]=rmse_v[e]/R_v

        #4. Calculate the pearson correlation coefficient.
        rho=np.corrcoef(rctr_v[:,0],orgl_v[:,0])
        rho_v[e]=rho[0,1]

        #Calculate the psnr
        maxI=np.max(orgl_v)
        minI=np.min(orgl_v)
        psnr_v[e]=20*np.log10((maxI-minI)/rmse_v[e])



    return elem,maxabserr_v,maxrelerr_v,rmse_v,nrmse_v,rho_v,psnr_v



def get_3Dmetrics(rctr_data,orgl_data,bm1_data,filename,isnap):

    #Get some information from the data
    nel=orgl_data.nel
    lxr=orgl_data.lr1
    nxyz=np.prod(lxr) 

    #Allocate vectors for the metrics
    elem=np.zeros(nel)
    maxabserr_v=np.zeros(nel)
    maxrelerr_v=np.zeros(nel)
    rmse_v=np.zeros(nel)
    nrmse_v=np.zeros(nel)
    rho_v=np.zeros(nel)
    psnr_v=np.zeros(nel)

    #Copy the nek data structure
    data=orgl_data
    abserr_data=copy.deepcopy(data)
    rms_data=copy.deepcopy(data)
    nrms_data=copy.deepcopy(data)
    rho_data=copy.deepcopy(data)
    psnr_data=copy.deepcopy(data)

    for e in range(0,nel):
        #Assign the element number
        elem[e]=e+1

        #Read the data in each element (ed=element data)
        rctr_ed=rctr_data.elem[e]  #for the reconstructed data
        bm1_ed=orgl_data.elem[e]  #for the original data

        #Do the same operations for the three velocity components
        for qoi in range(0,3):
            #Rearange the velocities in an easier vector to read.
            orgl_v=orgl_data.elem[e].vel[qoi,:,:,:].reshape((nxyz,1),order='F')
            rctr_v=rctr_data.elem[e].vel[qoi,:,:,:].reshape((nxyz,1),order='F')
            bm1_v=bm1_data.elem[e].vel[qoi,:,:,:].reshape((nxyz,1),order='F')
            R_v=np.max(orgl_v)-np.min(orgl_v)

            #1. Calculate the first metric (Absolute error)
            abserr_v=orgl_v-rctr_v
            #get the maximun to print only one value
            maxabserr_v[e]=np.max(np.abs(abserr_v))

            #2. Calculate the second metric (Relative error)
            relerr_v=abserr_v/R_v
            #get the maximun to print only one value
            maxrelerr_v[e]=np.max(np.abs(relerr_v))

            #3. Calculate the third metric (RMSE)
            #3.1 - weighted
            temp=0
            tempvol=0
            for i in range(0,nxyz):
                temp=temp+abserr_v[i,0]*abserr_v[i,0]*bm1_v[i,0]
                tempvol=tempvol+bm1_v[i,0]
            rmse_v[e]=np.sqrt(temp)/np.sqrt(tempvol)

            ##3.2 - not weighted
            #temp=np.matmul(abserr_v.T,abserr_v)/nxyz
            #rmse_v[e]=np.sqrt(temp[0,0])

            # Also calculate the normalized one
            nrmse_v[e]=rmse_v[e]/R_v

            #4. Calculate the pearson correlation coefficient.
            rho=np.corrcoef(rctr_v[:,0],orgl_v[:,0])
            rho_v[e]=rho[0,1]

            #5. Calculate the psnr
            maxI=np.max(orgl_v)
            minI=np.min(orgl_v)
            psnr_v[e]=20*np.log10((maxI-minI)/rmse_v[e])

            #copy into data structure to be read by nek
            abserr_data.elem[e].vel[qoi,:,:,:]=np.reshape(np.abs(abserr_v),(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')
            rms_data.elem[e].vel[qoi,:,:,:]=np.reshape(np.full((data.lr1[2]*data.lr1[1]*data.lr1[0],1),rmse_v[e]),(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')
            nrms_data.elem[e].vel[qoi,:,:,:]=np.reshape(np.full((data.lr1[2]*data.lr1[1]*data.lr1[0],1),nrmse_v[e]),(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')
            rho_data.elem[e].vel[qoi,:,:,:]=np.reshape(np.full((data.lr1[2]*data.lr1[1]*data.lr1[0],1),rho_v[e]),(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')
            psnr_data.elem[e].vel[qoi,:,:,:]=np.reshape(np.full((data.lr1[2]*data.lr1[1]*data.lr1[0],1),psnr_v[e]),(data.lr1[2],data.lr1[1],data.lr1[0]),order='F')


    writenek('abserr'+filename+'0.f'+str(isnap).zfill(5),abserr_data)
    writenek('rms'+filename+'0.f'+str(isnap).zfill(5),rms_data)
    writenek('nrms'+filename+'0.f'+str(isnap).zfill(5),nrms_data)
    writenek('rho'+filename+'0.f'+str(isnap).zfill(5),rho_data)
    writenek('psnr'+filename+'0.f'+str(isnap).zfill(5),psnr_data)

    return 

def plot_metrics(elem,maxabserr_v,maxrelerr_v,rmse_v,nrmse_v,rho_v,qoi):
    mpl.rcParams["text.usetex"]= 'true'
    mpl.rcParams["figure.figsize"]= [10,7]
    if qoi==0:
        tag='$v_x$'
    if qoi==1:
        tag='$v_y$'
    if qoi==2:
        tag='$v_z$'

    fig, axs = plt.subplots(2, 2)
    fig.suptitle(r'Metrics for '+tag)
    axs[0, 0].plot(elem, maxabserr_v)
    axs[0,0].set_ylabel(r'$e=X-\tilde{X}$')
    axs[0, 1].plot(elem, rmse_v,)
    axs[0,1].set_ylabel(r'$rmse$')
    axs[1, 0].plot(elem, nrmse_v, )
    axs[1,0].set_ylabel(r'$psnr$')
    axs[1, 1].plot(elem, rho_v, )
    axs[1,1].set_ylabel(r'$\rho$')
    #axs[0,0].set_ylim(-0.1,0.5)
    #axs[0,1].set_ylim(-0.1,0.5)
    #axs[1,0].set_ylim(-0.1,0.5)
    #axs[1,1].set_ylim(0,1.1)
    plt.show()
    
    return

#====================================================================

#Define the names of the files to compare
#path
orgl_path='../run_10e-5/alldata/'
rctr_path='../run_10e-3/olddata/'
#rctr_path='data/'
bm1_path='../run_10e-5/data/'
#case name
casename='turbPipe'
isnap=1
nsnap=50

#default names
orgl_filename=orgl_path+casename+'0.f'+str(isnap).zfill(5)
rctr_filename=rctr_path+'rct'+casename+'0.f'+str(isnap).zfill(5)
bm1_filename=bm1_path+'bm1'+casename+'0.f'+str(1).zfill(5)
out_filename='rct'+casename

#Read the whole files information
orgl_data=readnek(orgl_filename)
rctr_data=readnek(rctr_filename)
bm1_data=readnek(bm1_filename)

#get the 3D metrics
get_3Dmetrics(rctr_data,orgl_data,bm1_data,out_filename,isnap)

#get the metrics
qoi=2  #0-vx 1-vy 2-vz
[elemz,maxabserr_vz,maxrelerr_vz,rmse_vz,nrmse_vz,rho_vz,psnr_vz]=get_metrics(rctr_data,orgl_data,bm1_data,qoi)

#Plot the metrics
plot_metrics(elemz,maxabserr_vz,maxrelerr_vz,rmse_vz,psnr_vz,rho_vz,qoi)
