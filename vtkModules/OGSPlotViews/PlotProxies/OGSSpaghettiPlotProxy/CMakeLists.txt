set(classes OGSSpaghettiPlotProxy)

vtk_module_add_module(OGS::OGSSpaghettiPlotProxy CLASSES ${classes})
vtk_module_set_property(OGS::OGSSpaghettiPlotProxy PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

vtk_module_link(OGS::OGSSpaghettiPlotProxy PUBLIC OGS::OGSSpaghettiPlot)
