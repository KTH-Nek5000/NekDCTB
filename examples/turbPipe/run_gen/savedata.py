
#--------------------------------------------------------------

# Global timestep output options
timeStepToStartOutputAt=0
forceOutputAtFirstCall=False

# Global screenshot output options
imageFileNamePadding=0
rescale_lookuptable=False

# Whether or not to request specific arrays from the adaptor.
requestSpecificArrays=False

# a root directory under which all Catalyst output goes
rootDirectory=''

# makes a cinema D index table
make_cinema_table=False

# TODO: Change this number for the real case (to e.g. 10)
global_freq = 2

#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# paraview version 5.6.0
#--------------------------------------------------------------

from paraview.simple import *
from paraview import coprocessing

# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.6.0

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # trace generated using paraview version 5.6.0
      #
      # To ensure correct image size when batch processing, please search
      # for and uncomment the line `# renderView*.ViewSize = [*,*]`

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # create a new 'PVTrivialProducer'
      # create a producer from a simulation input
      input = coprocessor.CreateProducer(datadescription, 'input')

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(input)
      # ----------------------------------------------------------------

      # Now any catalyst writers
      xMLPUnstructuredGridWriter1 = servermanager.writers.XMLPUnstructuredGridWriter(Input=input)
      coprocessor.RegisterWriter(xMLPUnstructuredGridWriter1, filename='data/input_%t.pvtu', freq=global_freq)

    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'input': [global_freq]}
  coprocessor.SetUpdateFrequencies(freqs)
  if requestSpecificArrays:
    arrays = [['pressure', 0], ['temperature', 0], ['velocity', 0]]
    coprocessor.SetRequestedArrays('input', arrays)
  coprocessor.SetInitialOutputOptions(timeStepToStartOutputAt,forceOutputAtFirstCall)

  if rootDirectory:
      coprocessor.SetRootDirectory(rootDirectory)

  if make_cinema_table:
      coprocessor.EnableCinemaDTable()

  return coprocessor


#--------------------------------------------------------------
# Global variable that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView and the update frequency
coprocessor.EnableLiveVisualization(False, 1)

# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

# TODO: Change this number for the real case (to e.g. 10 and 1)
samples_per_mode = 2
break_between_modes = 2

# Mode: Save data
startTime = 1 + 27 * global_freq * (samples_per_mode-1 + break_between_modes+1)
modeStartTimes = [startTime]
modeEndTimes = [startTime + global_freq*(samples_per_mode-1)]
endTime = modeEndTimes[-1]
numModes = len(modeEndTimes)

def setTimeSettings(set_samples_per_mode, set_break_between_modes, set_startTime):
    global samples_per_mode, break_between_modes
    global modeStartTimes, modeEndTimes
    global startTime, endTime
    samples_per_mode = set_samples_per_mode
    break_between_modes = set_break_between_modes
    startTime = set_startTime
    modeStartTimes = [startTime]
    modeEndTimes = [startTime + global_freq*(samples_per_mode-1)]
    endTime = modeEndTimes[-1]

def printModes():
    print("Time: [%i, %i],\tData: All, Mode: Save."%(modeStartTimes[0],modeEndTimes[0]))


def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=rescale_lookuptable,
        image_quality=0, padding_amount=imageFileNamePadding)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
