# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
pv_version = GetParaViewVersion()
assert pv_version >= 5.7, "Wrong paraview version, we need 5.7.0."

import argparse
parser = argparse.ArgumentParser(description='Plotting optimized design')
parser.add_argument('--results_dir', dest="results_dir", type=str, help='Results directory', default='./', required=True)
parser.add_argument('--filename', dest="filename", type=str, help='filename for the screenshot', default="optimized_design", required=True)

args = parser.parse_args()

results_dir = args.results_dir
filename = args.filename

# create a new 'PVD Reader'
control_iterations_fpvd = PVDReader(FileName=results_dir + "/" + filename + ".pvd")
control_iterations_fpvd.PointArrays = ['rho filtered']

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1804, 871]

# show data in view
control_iterations_fpvdDisplay = Show(control_iterations_fpvd, renderView1)

# get color transfer function/color map for 'rho'
rhoLUT = GetColorTransferFunction('rho filtered')

# get opacity transfer function/opacity map for 'rho'
rhoPWF = GetOpacityTransferFunction('rho filtered')

# trace defaults for the display properties.
control_iterations_fpvdDisplay.Representation = 'Surface'
control_iterations_fpvdDisplay.ColorArrayName = ['POINTS', 'rho filtered']
control_iterations_fpvdDisplay.LookupTable = rhoLUT
control_iterations_fpvdDisplay.OSPRayScaleArray = 'rho filtered'
control_iterations_fpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
control_iterations_fpvdDisplay.SelectOrientationVectors = 'None'
control_iterations_fpvdDisplay.ScaleFactor = 0.16
control_iterations_fpvdDisplay.SelectScaleArray = 'rho filtered'
control_iterations_fpvdDisplay.GlyphType = 'Arrow'
control_iterations_fpvdDisplay.GlyphTableIndexArray = 'rho filtered'
control_iterations_fpvdDisplay.GaussianRadius = 0.008
control_iterations_fpvdDisplay.SetScaleArray = ['POINTS', 'rho filtered']
control_iterations_fpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
control_iterations_fpvdDisplay.OpacityArray = ['POINTS', 'rho filtered']
control_iterations_fpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
control_iterations_fpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
control_iterations_fpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
control_iterations_fpvdDisplay.ScalarOpacityFunction = rhoPWF
control_iterations_fpvdDisplay.ScalarOpacityUnitDistance = 0.09367875016563294

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
control_iterations_fpvdDisplay.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
control_iterations_fpvdDisplay.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.6, 0.75, 10000.0]
renderView1.CameraFocalPoint = [0.6, 0.75, 0.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# show color bar/color legend
control_iterations_fpvdDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

animationScene1.GoToLast()

# rescale color and/or opacity maps used to exactly fit the current data range
control_iterations_fpvdDisplay.RescaleTransferFunctionToDataRange(False, True)

# hide color bar/color legend
control_iterations_fpvdDisplay.SetScalarBarVisibility(renderView1, False)

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# create a new 'Text'
text1 = Text()

# Properties modified on text1
text1.Text = filename

# show data in view
text1Display = Show(text1, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on text1Display
text1Display.WindowLocation = 'AnyLocation'
text1Display.Position = [0.570573, 0.221482]
text1Display.Color = [0.0, 0.0, 0.0]
text1Display.FontSize = 9


# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(control_iterations_fpvd)

#### saving camera placements for all active views

# current camera placement for renderView1

renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [50.0, 50.0, 273.20508075688775]
renderView1.CameraFocalPoint = [50.0, 50.0, 0.0]
renderView1.CameraParallelScale = 70.71067811865476

# save screenshot
SaveScreenshot(results_dir + "/" + filename + ".png", renderView1, ImageResolution=[1804, 871])

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
