### Modified script for extracting data along a 
### given line using lineout operator.
### Extracts data upto 14 decimal places.

#Command to run from terminal (without GUI) : visit -nowin -cli -s lineout-extract-modified.py
#Author  : Anant Diwakar
#Date    : 14 April 2015

var_list = ["Density","Pressure","Energy"] # List of variable to be extracted

no_of_points = 100 # No. of points on the curve

#Lineout parameters
x_start = 0.0
x_end = 1.0
dx = (x_end - x_start)/no_of_points
y_start = 5*dx
y_end = y_start

OpenDatabase("./solution-*.vtu database", 0)
SetTimeSliderState(TimeSliderGetNStates() - 1) # Move the slider to last time

for var in var_list:

    var_path = "operators/Lineout/%s" %(var)
    AddPlot("Curve",var_path, 1, 1)

    LineoutAtts = LineoutAttributes()
    LineoutAtts.point1 = (x_start, y_start, 0)
    LineoutAtts.point2 = (x_end, y_end, 0)
    LineoutAtts.interactive = 0
    LineoutAtts.ignoreGlobal = 1
    LineoutAtts.samplingOn = 1
    LineoutAtts.numberOfSamplePoints = no_of_points
    LineoutAtts.reflineLabels = 0
    SetOperatorOptions(LineoutAtts, 1)
    DrawPlots()
    SaveWindowAtts = SaveWindowAttributes()
    SaveWindowAtts.outputToCurrentDirectory = 1
    SaveWindowAtts.outputDirectory = "."
    SaveWindowAtts.fileName = var
    SaveWindowAtts.family = 0
    SaveWindowAtts.format = SaveWindowAtts.CURVE
    SetSaveWindowAttributes(SaveWindowAtts)
    SaveWindow()
    HideActivePlots()

exit()
