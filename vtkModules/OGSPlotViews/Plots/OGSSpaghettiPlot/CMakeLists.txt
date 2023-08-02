set(classes OGSSpaghettiPlot)

vtk_module_add_module(OGS::OGSSpaghettiPlot CLASSES ${classes})
vtk_module_set_property(OGS::OGSSpaghettiPlot PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

