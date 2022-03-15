# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
import glob
import re
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
pv_version = GetParaViewVersion()
assert pv_version >= 5.7, "Wrong paraview version, we need 5.7.0."

import argparse
parser = argparse.ArgumentParser(description='Plotting optimized design')
parser.add_argument('--results_dir', dest="results_dir", type=str, help='Results directory', default='./', required=True)
parser.add_argument('--filename', dest="filename", type=str, help='filename for the screenshot', default="optimized_design", required=True)
parser.add_argument('--output_name', dest="output_name", type=str, help='output_name for the screenshot', default="optimized_design", required=True)

args = parser.parse_args()

results_dir = args.results_dir
print(results_dir)
filename = args.filename

# create a new 'PVD Reader'
all_files = glob.glob(results_dir + "/" + filename + "*")
print(all_files)
all_files = sorted(all_files, key=lambda L: list(map(int, re.findall('(\d+).pvd', L))))

control_iterations_fpvd = PVDReader(FileName=all_files[-1])
control_iterations_fpvd.PointArrays = ['rho']

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
rhoLUT = GetColorTransferFunction('rho')

# get opacity transfer function/opacity map for 'rho'
rhoPWF = GetOpacityTransferFunction('rho')

# trace defaults for the display properties.
control_iterations_fpvdDisplay.Representation = 'Surface'
control_iterations_fpvdDisplay.ColorArrayName = ['POINTS', 'rho']
control_iterations_fpvdDisplay.LookupTable = rhoLUT
control_iterations_fpvdDisplay.OSPRayScaleArray = 'rho'
control_iterations_fpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
control_iterations_fpvdDisplay.SelectOrientationVectors = 'None'
control_iterations_fpvdDisplay.ScaleFactor = 0.16
control_iterations_fpvdDisplay.SelectScaleArray = 'rho'
control_iterations_fpvdDisplay.GlyphType = 'Arrow'
control_iterations_fpvdDisplay.GlyphTableIndexArray = 'rho'
control_iterations_fpvdDisplay.GaussianRadius = 0.008
control_iterations_fpvdDisplay.SetScaleArray = ['POINTS', 'rho']
control_iterations_fpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
control_iterations_fpvdDisplay.OpacityArray = ['POINTS', 'rho']
control_iterations_fpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
control_iterations_fpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
control_iterations_fpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
control_iterations_fpvdDisplay.ScalarOpacityFunction = rhoPWF
control_iterations_fpvdDisplay.ScalarOpacityUnitDistance = 0.09367875016563294

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
rhoLUT.ApplyPreset('X Ray', True)

# Properties modified on levelSetLUT
rhoLUT.EnableOpacityMapping = 1

LoadPalette(paletteName='WhiteBackground')


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

# reset view to fit data
renderView1.ResetCamera()

# save screenshot
SaveScreenshot(results_dir + "/" + args.output_name + ".png", renderView1, ImageResolution=[1804, 871])


# create a new 'Calculator'
calculator1 = Calculator(Input=control_iterations_fpvd)
calculator1.Function = ''

# Properties modified on calculator1
calculator1.Function = '0'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1627, 880]

# get layout
layout1 = GetLayout()

# show data in view
calculator1Display = Show(calculator1, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'Result'
resultLUT = GetColorTransferFunction('Result')

# get opacity transfer function/opacity map for 'Result'
resultPWF = GetOpacityTransferFunction('Result')

# trace defaults for the display properties.
calculator1Display.Representation = 'Surface'
calculator1Display.ColorArrayName = ['POINTS', 'Result']
calculator1Display.LookupTable = resultLUT
calculator1Display.OSPRayScaleArray = 'Result'
calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator1Display.SelectOrientationVectors = 'None'
calculator1Display.ScaleFactor = 10.0
calculator1Display.SelectScaleArray = 'Result'
calculator1Display.GlyphType = 'Arrow'
calculator1Display.GlyphTableIndexArray = 'Result'
calculator1Display.GaussianRadius = 0.5
calculator1Display.SetScaleArray = ['POINTS', 'Result']
calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator1Display.OpacityArray = ['POINTS', 'Result']
calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
calculator1Display.PolarAxes = 'PolarAxesRepresentation'
calculator1Display.ScalarOpacityFunction = resultPWF
calculator1Display.ScalarOpacityUnitDistance = 1.6065154272065734

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
calculator1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
calculator1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# hide data in view
Hide(control_iterations_fpvd, renderView1)

# show color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, False)

# update the view to ensure updated data information
renderView1.Update()

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
resultLUT.ApplyPreset('X Ray', True)

# Properties modified on resultLUT
resultLUT.EnableOpacityMapping = 1

# change representation type
calculator1Display.SetRepresentationType('Surface With Edges')

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# save screenshot
SaveScreenshot(results_dir + "/" + args.output_name + "_grid.png", renderView1, ImageResolution=[7216, 3484])

#### saving camera placements for all active views
#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
