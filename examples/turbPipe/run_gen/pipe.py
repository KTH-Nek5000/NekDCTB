
#--------------------------------------------------------------

# Global timestep output options
timeStepToStartOutputAt=0
forceOutputAtFirstCall=False

# Global screenshot output options
imageFileNamePadding=5
rescale_lookuptable=False

# Whether or not to request specific arrays from the adaptor.
requestSpecificArrays=False

# a root directory under which all Catalyst output goes
rootDirectory=''

# makes a cinema D index table
make_cinema_table=False

#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# paraview version 5.6.2
#--------------------------------------------------------------

from paraview.simple import *
from paraview import coprocessing
from mpi4py import MPI

#import vmslice, vmisosurface, vmvolrender
#import vmslice, vmisosurface,
#import vmslice, savedata
import vmslice

# ---------------------- Pipeline Mode Setup ----------------------

global_freq = 2
samples_per_mode = 101
break_between_modes = 9
# Do one break before the first mode
# First timestep in which visualization is called is global_freq 
#(simulation steps start with 1, Catalyst seems to expect 0)
startTime = global_freq + break_between_modes * global_freq + 100

# Samples are inclusive with respect to start and end (->-1), break are exclusive (->=+1)

vmslice.setTimeSettings(samples_per_mode, break_between_modes, startTime)

#startTime += vmslice.numModes * global_freq * (samples_per_mode-1 + break_between_modes+1)
'''
vmisosurface.setTimeSettings(samples_per_mode, break_between_modes, startTime)

#startTime += vmisosurface.numModes * global_freq * (samples_per_mode-1 + break_between_modes+1)
#vmvolrender.setTimeSettings(samples_per_mode, break_between_modes, startTime)

startTime += vmisosurface.numModes * global_freq * (samples_per_mode-1 + break_between_modes+1)
'''
#savedata.setTimeSettings(samples_per_mode, break_between_modes, startTime)

if (MPI.COMM_WORLD.Get_rank() == 0):
    print("############## Visualization Modes ##############")
    vmslice.printModes()
    #vmisosurface.printModes()
    #vmvolrender.printModes()
    #savedata.printModes()
    print("#################################################")


# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global vmslice
    #global vmisosurface
    #global vmvolrender
    #global savedata

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    vmslice.RequestDataDescription(datadescription)
    #vmisosurface.RequestDataDescription(datadescription)
    #vmvolrender.RequestDataDescription(datadescription)
    #savedata.RequestDataDescription(datadescription)


# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global vmslice
    #global vmisosurface
    #global vmvolrender
    #global savedata

    timestep = datadescription.GetTimeStep()
    '''    
    if (timestep >= vmslice.startTime and timestep <= vmslice.endTime):
        vmslice.DoCoProcessing(datadescription)

    if (timestep >= vmisosurface.startTime and timestep <= vmisosurface.endTime):
        vmisosurface.DoCoProcessing(datadescription)
    
    #if (timestep >= vmvolrender.startTime and timestep <= vmvolrender.endTime):
    #    vmvolrender.DoCoProcessing(datadescription)
    if (timestep >= savedata.startTime and timestep <= savedata.endTime):
        savedata.DoCoProcessing(datadescription)
    '''
    vmslice.DoCoProcessing(datadescription)
    #savedata.DoCoProcessing(datadescription)
