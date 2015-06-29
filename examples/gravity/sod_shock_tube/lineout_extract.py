#Command : visit -cli -nowin -s lineout_extract.py
#Author  : Anant Diwakar
#Date    : 28 June 2014

######################################################

axis_font_size = 1.1

no_of_points = 100 # No. of points on the curve

#Lineout parameters
x_start = 0.0
x_end = 1.0
y_start = 0.05
y_end = 0.05

var_list = ["Density","Pressure","Energy"] # List of variable to be extracted

#######################################################

OpenDatabase("./solution-*.vtu database", 0)
SetTimeSliderState(TimeSliderGetNStates() - 1) # Move the slider to last time 

for var in var_list:

    var_path = """operators/Lineout/%s""" % (var)
    AddPlot("Curve", var_path, 1, 1)
    
    LineoutAtts = LineoutAttributes()
    LineoutAtts.point1 = (x_start, y_start, 0)
    LineoutAtts.point2 = (x_end, y_end, 0)
    LineoutAtts.interactive = 1
    LineoutAtts.ignoreGlobal = 1
    LineoutAtts.samplingOn = 1
    LineoutAtts.numberOfSamplePoints = no_of_points
    LineoutAtts.reflineLabels = 0
    SetOperatorOptions(LineoutAtts, 1)
    CurveAtts = CurveAttributes()
    CurveAtts.showLines = 1
    CurveAtts.lineStyle = CurveAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
    CurveAtts.lineWidth = 1
    CurveAtts.showPoints = 1
    CurveAtts.symbol = CurveAtts.Circle  # Point, TriangleUp, TriangleDown, Square, Circle, Plus, X
    CurveAtts.pointSize = 2
    CurveAtts.pointFillMode = CurveAtts.Static  # Static, Dynamic
    CurveAtts.pointStride = 1
    CurveAtts.symbolDensity = 50
    CurveAtts.curveColorSource = CurveAtts.Cycle  # Cycle, Custom
    CurveAtts.curveColor = (255, 0, 0, 255)
    CurveAtts.showLegend = 1
    CurveAtts.showLabels = 0
    CurveAtts.designator = ""
    CurveAtts.doBallTimeCue = 0
    CurveAtts.ballTimeCueColor = (0, 0, 0, 255)
    CurveAtts.timeCueBallSize = 0.01
    CurveAtts.doLineTimeCue = 0
    CurveAtts.lineTimeCueColor = (0, 0, 0, 255)
    CurveAtts.lineTimeCueWidth = 0
    CurveAtts.doCropTimeCue = 0
    CurveAtts.timeForTimeCue = 0
    CurveAtts.fillMode = CurveAtts.NoFill  # NoFill, Solid, HorizontalGradient, VerticalGradient
    CurveAtts.fillColor1 = (255, 0, 0, 255)
    CurveAtts.fillColor2 = (255, 100, 100, 255)
    CurveAtts.polarToCartesian = 0
    CurveAtts.polarCoordinateOrder = CurveAtts.R_Theta  # R_Theta, Theta_R
    CurveAtts.angleUnits = CurveAtts.Radians  # Radians, Degrees
    SetPlotOptions(CurveAtts)

    # End spontaneous state

    ViewCurveAtts = ViewCurveAttributes()
    ViewCurveAtts.domainCoords = (x_start, x_end)
    ViewCurveAtts.rangeCoords = (-1, 3)
    ViewCurveAtts.viewportCoords = (0.2, 0.95, 0.15, 0.9)
    ViewCurveAtts.domainScale = ViewCurveAtts.LINEAR  # LINEAR, LOG
    ViewCurveAtts.rangeScale = ViewCurveAtts.LINEAR  # LINEAR, LOG
    SetViewCurve(ViewCurveAtts)
    # Logging for SetAnnotationObjectOptions is not implemented yet.
    AnnotationAtts = AnnotationAttributes()
    AnnotationAtts.axes2D.visible = 1
    AnnotationAtts.axes2D.xAxis.title.font.scale = axis_font_size
    AnnotationAtts.axes2D.xAxis.label.font.scale = axis_font_size
    AnnotationAtts.axes2D.yAxis.title.font.scale = axis_font_size
    AnnotationAtts.axes2D.yAxis.label.font.scale = axis_font_size
    AnnotationAtts.userInfoFlag = 0
    AnnotationAtts.databaseInfoFlag = 1
    AnnotationAtts.timeInfoFlag = 1
    AnnotationAtts.legendInfoFlag = 1
    SetAnnotationAttributes(AnnotationAtts)
    #ResetView()
    DrawPlots()

    ####################################################
    # Extracting the variables to respective data files

    ExportDBAtts = ExportDBAttributes()
    ExportDBAtts.db_type = "Curve2D"
    ExportDBAtts.filename = var    # variable file name
    ExportDBAtts.dirname = "."
    ExportDBAtts.variables = ()
    ExportDBAtts.opts.types = ()
    ExportDatabase(ExportDBAtts)

exit()
