set(classes OGSHovmoellerPlot)

vtk_module_add_module(OGS::OGSHovmoellerPlot CLASSES ${classes})
vtk_module_set_property(OGS::OGSHovmoellerPlot PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

