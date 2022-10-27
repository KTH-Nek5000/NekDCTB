

#--------------------------------------------------------------



# Global timestep output options

timeStepToStartOutputAt=0

forceOutputAtFirstCall=True



# Global screenshot output options

imageFileNamePadding=5

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

# paraview version 5.6.2

#--------------------------------------------------------------



from paraview.simple import *

from paraview import coprocessing



# ----------------------- CoProcessor definition -----------------------



def CreateCoProcessor():

  def _CreatePipeline(coprocessor, datadescription):

    class Pipeline:

      # state file generated using paraview version 5.6.2



      # ----------------------------------------------------------------

      # setup views used in the visualization

      # ----------------------------------------------------------------



      # trace generated using paraview version 5.6.2

      #

      # To ensure correct image size when batch processing, please search 

      # for and uncomment the line `# renderView*.ViewSize = [*,*]`



      #### disable automatic camera reset on 'Show'

      paraview.simple._DisableFirstRenderCameraReset()



      # get the material library

      #materialLibrary1 = GetMaterialLibrary()



      # Create a new 'Render View'

      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [960, 540]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.StereoType = 'Crystal Eyes'
      renderView1.CameraFocalDisk = 1.0
      renderView1.CameraPosition = [16.69621680161044, 0.0, 4.083334922790527]
      renderView1.CameraFocalPoint = [0.0, 0.0, 4.083334922790527]
      renderView1.CameraViewUp = [0.0, -1.0, 2.220446049250313e-16]
      renderView1.CameraParallelScale = 4.321298889417477


      # register the view with coprocessor

      # and provide it with information such as the filename to use,

      # how frequently to write the images, etc.

      coprocessor.RegisterView(renderView1,

          filename='fig/Pipe_VelocityMagnitude_Slice_%t.png', freq=global_freq, fittoscreen=0, magnification=1, width=1920, height=1280, cinema={})

      renderView1.ViewTime = datadescription.GetTime()



      # ----------------------------------------------------------------

      # restore active view

      SetActiveView(renderView1)

      # ----------------------------------------------------------------



      # ----------------------------------------------------------------

      # setup the data processing pipelines

      # ----------------------------------------------------------------



      # create a new 'XML Partitioned Unstructured Grid Reader'

      # create a producer from a simulation input

      input = coprocessor.CreateProducer(datadescription, 'input')


      # show data in view

      inputpvtuDisplay = Show(input, renderView1, 'UnstructuredGridRepresentation')



      # trace defaults for the display properties.

      inputpvtuDisplay.Representation = 'Surface'
      inputpvtuDisplay.ColorArrayName = [None, '']
      inputpvtuDisplay.SelectTCoordArray = 'None'
      inputpvtuDisplay.SelectNormalArray = 'None'
      inputpvtuDisplay.SelectTangentArray = 'None'
      inputpvtuDisplay.OSPRayScaleArray = 'pressure'
      inputpvtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      inputpvtuDisplay.SelectOrientationVectors = 'None'
      inputpvtuDisplay.ScaleFactor = 0.8166669845581055
      inputpvtuDisplay.SelectScaleArray = 'None'
      inputpvtuDisplay.GlyphType = 'Arrow'
      inputpvtuDisplay.GlyphTableIndexArray = 'None'
      inputpvtuDisplay.GaussianRadius = 0.040833349227905276
      inputpvtuDisplay.SetScaleArray = ['POINTS', 'pressure']
      inputpvtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      inputpvtuDisplay.OpacityArray = ['POINTS', 'pressure']
      inputpvtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
      inputpvtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
      inputpvtuDisplay.PolarAxes = 'PolarAxesRepresentation'
      inputpvtuDisplay.ScalarOpacityUnitDistance = 0.03518722456811719
      inputpvtuDisplay.OpacityArrayName = ['POINTS', 'pressure']


# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'

      inputpvtuDisplay.ScaleTransferFunction.Points = [-0.5385369962813863, 0.0, 0.5, 0.0, 1.69633290316386, 1.0, 0.5, 0.0]


# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'

      inputpvtuDisplay.OpacityTransferFunction.Points = [-0.5385369962813863, 0.0, 0.5, 0.0, 1.69633290316386, 1.0, 0.5, 0.0]



      # show data from zSlice

      #zSliceDisplay = Show(zSlice, renderView1)

      slice1 = Slice(registrationName='Slice1', Input=input)
      slice1.SliceType = 'Plane'
      slice1.HyperTreeGridSlicer = 'Plane'
      slice1.SliceOffsetValues = [0.0]



      # init the 'Plane' selected for 'SliceType'

      slice1.SliceType.Origin = [0.0, 0.0, 4.083334922790527]


      # init the 'Plane' selected for 'HyperTreeGridSlicer'

      slice1.HyperTreeGridSlicer.Origin = [0.0, 0.0, 4.083334922790527]

      # get color transfer function/color map for 'velocity'

      velocityLUT = GetColorTransferFunction('velocity')

      # show data in view

      slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')



      # trace defaults for the display properties.

      slice1Display.Representation = 'Surface'
      slice1Display.ColorArrayName = ['POINTS', 'velocity']
      slice1Display.LookupTable = velocityLUT
      slice1Display.SelectTCoordArray = 'None'
      slice1Display.SelectNormalArray = 'None'
      slice1Display.SelectTangentArray = 'None'
      slice1Display.OSPRayScaleArray = 'pressure'
      slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      slice1Display.SelectOrientationVectors = 'None'
      slice1Display.ScaleFactor = 0.8166669845581055
      slice1Display.SelectScaleArray = 'None'
      slice1Display.GlyphType = 'Arrow'
      slice1Display.GlyphTableIndexArray = 'None'
      slice1Display.GaussianRadius = 0.040833349227905276
      slice1Display.SetScaleArray = ['POINTS', 'pressure']
      slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
      slice1Display.OpacityArray = ['POINTS', 'pressure']
      slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
      slice1Display.DataAxesGrid = 'GridAxesRepresentation'
      slice1Display.PolarAxes = 'PolarAxesRepresentation'



      # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'

      slice1Display.ScaleTransferFunction.Points = [-0.5202897808174186, 0.0, 0.5, 0.0, 1.5198318261104378, 1.0, 0.5, 0.0]



      # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'

      slice1Display.OpacityTransferFunction.Points = [-0.5202897808174186, 0.0, 0.5, 0.0, 1.5198318261104378, 1.0, 0.5, 0.0]

      # ----------------------------------------------------------------

      # setup color maps and opacity mapes used in the visualization

      # note: the Get..() functions create a new object, if needed

      # ----------------------------------------------------------------
      Hide(input, renderView1)

      slice1Display.SetScalarBarVisibility(renderView1, True)
      


      # get opacity transfer function/opacity map for 'velocity'

      velocityPWF = GetOpacityTransferFunction('velocity')

      '''

      velocityPWF.Points = [0.0, 0.0, 0.5, 0.0, 0.425, 1.0, 0.5, 0.0, 0.85, 0.0, 0.5, 0.0, 1.5, 0.0, 0.5, 0.0, 1.5, 0.0, 0.5, 0.0]

      velocityPWF.ScalarRangeInitialized = 1

      '''

      # ----------------------------------------------------------------

      # finally, restore active source

      SetActiveSource(slice1)

      # ----------------------------------------------------------------

    return Pipeline()



  class CoProcessor(coprocessing.CoProcessor):

    def CreatePipeline(self, datadescription):

      self.Pipeline = _CreatePipeline(self, datadescription)



  coprocessor = CoProcessor()

  # these are the frequencies at which the coprocessor updates.

  freqs = {'input': [global_freq, global_freq]}

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



# Resolution settings:

# - low res 4:3 

lowres = [720, 480]

# - medium res 16:9

midres = [1280, 720]

# - high res 16:9

highres = [1920, 1080]



# Modes:

# For each resolution (lowres, midres, highres)

# - Compute slice without rendering

# - Render slice without saving

# - Save slice

startTime = 1

modeStartTimes = [startTime + x * global_freq * (samples_per_mode-1 + break_between_modes+1) for x in range(9)]

modeEndTimes = [x + global_freq*(samples_per_mode-1) for x in modeStartTimes]

endTime = modeEndTimes[-1]

numModes = len(modeEndTimes)



def setTimeSettings(set_samples_per_mode, set_break_between_modes, set_startTime):

    global samples_per_mode, break_between_modes

    global modeStartTimes, modeEndTimes

    global startTime, endTime

    samples_per_mode = set_samples_per_mode

    break_between_modes = set_break_between_modes

    startTime = set_startTime

    modeStartTimes = [startTime + x * global_freq * (samples_per_mode-1 + break_between_modes+1) for x in range(9)]

    modeEndTimes = [x + global_freq*(samples_per_mode-1) for x in modeStartTimes]

    endTime = modeEndTimes[-1]



def printModes():

    print("Time: [%i, %i], Resolution: [%i, %i],\tData: Velocity Magnitude, Mode: Compute slice without rendering."%(modeStartTimes[0],modeEndTimes[0],lowres[0],lowres[1]))

    print("Time: [%i, %i], Resolution: [%i, %i],\tData: Velocity Magnitude, Mode: Compute slice without rendering."%(modeStartTimes[1],modeEndTimes[1],midres[0],midres[1]))

    print("Time: [%i, %i], Resolution: [%i, %i],\tData: Velocity Magnitude, Mode: Compute slice without rendering."%(modeStartTimes[2],modeEndTimes[2],highres[0],highres[1]))

    

    print("Time: [%i, %i], Resolution: [%i, %i],\tData: Velocity Magnitude, Mode: Render slice without saving."%(modeStartTimes[3],modeEndTimes[3],lowres[0],lowres[1]))

    print("Time: [%i, %i], Resolution: [%i, %i],\tData: Velocity Magnitude, Mode: Render slice without saving."%(modeStartTimes[4],modeEndTimes[4],midres[0],midres[1]))

    print("Time: [%i, %i], Resolution: [%i, %i],\tData: Velocity Magnitude, Mode: Render slice without saving."%(modeStartTimes[5],modeEndTimes[5],highres[0],highres[1]))

    

    print("Time: [%i, %i], Resolution: [%i, %i],\tData: Velocity Magnitude, Mode: Save slice."%(modeStartTimes[6],modeEndTimes[6],lowres[0],lowres[1]))

    print("Time: [%i, %i], Resolution: [%i, %i],\tData: Velocity Magnitude, Mode: Save slice."%(modeStartTimes[7],modeEndTimes[7],midres[0],midres[1]))

    print("Time: [%i, %i], Resolution: [%i, %i],\tData: Velocity Magnitude, Mode: Save slice."%(modeStartTimes[8],modeEndTimes[8],highres[0],highres[1]))





def DoCoProcessing(datadescription):

    "Callback to do co-processing for current timestep"

    global coprocessor



    # Update the coprocessor by providing it the newly generated simulation data.

    # If the pipeline hasn't been setup yet, this will setup the pipeline.

    coprocessor.UpdateProducers(datadescription)



    timestep = datadescription.GetTimeStep()


    '''
    if ((timestep >= modeStartTimes[0] and timestep <= modeEndTimes[0]) or 

        (timestep >= modeStartTimes[3] and timestep <= modeEndTimes[3]) or

        (timestep >= modeStartTimes[6] and timestep <= modeEndTimes[6])):

        coprocessor.Pipeline.renderView1.ViewSize = lowres



    elif ((timestep >= modeStartTimes[1] and timestep <= modeEndTimes[1]) or 

        (timestep >= modeStartTimes[4] and timestep <= modeEndTimes[4]) or

        (timestep >= modeStartTimes[7] and timestep <= modeEndTimes[7])):

        coprocessor.Pipeline.renderView1.ViewSize = midres

    

    elif ((timestep >= modeStartTimes[2] and timestep <= modeEndTimes[2]) or 

        (timestep >= modeStartTimes[5] and timestep <= modeEndTimes[5]) or

        (timestep >= modeStartTimes[8] and timestep <= modeEndTimes[8])):

        coprocessor.Pipeline.renderView1.ViewSize = highres

    else:

        return



    if (timestep <= modeEndTimes[2]):

        # Compute one the slice without rendering it

        coprocessor.Pipeline.slice1.UpdatePipeline()

    elif (timestep <= modeEndTimes[5]):



        # Viewtime needs to be updated here, so not cause errors

        coprocessor.Pipeline.renderView1.ViewTime = datadescription.GetTime()

        # Render slice, triggers computation as well

        Render(view=coprocessor.Pipeline.renderView1)

    else:

        # Write output data, if appropriate.

        coprocessor.WriteData(datadescription);



        # Write image capture (Last arg: rescale lookup table), if appropriate.

        coprocessor.WriteImages(datadescription, rescale_lookuptable=rescale_lookuptable,

          image_quality=0, padding_amount=imageFileNamePadding)

    '''
    #coprocessor.Pipeline.renderView1.ViewSize = midres
    coprocessor.WriteData(datadescription)
    coprocessor.WriteImages(datadescription, rescale_lookuptable=rescale_lookuptable, image_quality=0, padding_amount=imageFileNamePadding)
    # Live Visualization, if enabled.
    
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)



