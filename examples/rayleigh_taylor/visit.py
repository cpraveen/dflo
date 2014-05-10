#Window size
windowWidth = 350
windowHeight = 800

#Plot bounds
plot_bound_x1 = 0.35
plot_bound_x2 = 0.9
plot_bound_y1 = 0.15
plot_bound_y2 = 0.9

#Legend parameters
legend_fontHeight = 0.02
legend_xScale = 1
legend_yScale = 1
legend_position = (0.01,0.5)

MoveAndResizeWindow(1, 0, 0, windowWidth, windowHeight)

OpenDatabase("./solution-*.vtu database", 0)
AddPlot("Pseudocolor", "Density", 1, 1)

# Begin spontaneous state
View2DAtts = View2DAttributes()
View2DAtts.windowCoords = (-0.25, 0.25, -0.75, 0.75)
View2DAtts.viewportCoords = (plot_bound_x1, plot_bound_x2, plot_bound_y1, plot_bound_y2)
View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
View2DAtts.fullFrameAutoThreshold = 100
View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
View2DAtts.windowValid = 1
SetView2D(View2DAtts)

legend=GetAnnotationObject("Plot0000")
legend.position = legend_position
legend.fontHeight = legend_fontHeight
legend.xScale = legend_xScale
legend.yScale = legend_yScale

DrawPlots()
