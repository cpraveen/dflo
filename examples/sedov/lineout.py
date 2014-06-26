#Command : visit -s lineout.py
#Author  : Anant Diwakar
#Date    : 26 June 2014

######################################################

axis_font_size = 1.1

no_of_points = 500 # No. of points on the curve

#Lineout parameters
x_start = -2.0
x_end = 2.0
y_start = 0.0
y_end = 0.0

#######################################################
OpenDatabase("./solution-*.vtu database", 0)
AddPlot("Curve", "operators/Lineout/Density", 1, 1)
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
ViewCurveAtts.domainCoords = (7.45058e-10, 1)
ViewCurveAtts.rangeCoords = (0.125, 1)
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

DrawPlots()
