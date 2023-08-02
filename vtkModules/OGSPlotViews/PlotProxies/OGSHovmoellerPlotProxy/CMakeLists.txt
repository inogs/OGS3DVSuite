set(classes OGSHovmoellerPlotProxy)

vtk_module_add_module(OGS::OGSHovmoellerPlotProxy CLASSES ${classes})
vtk_module_set_property(OGS::OGSHovmoellerPlotProxy PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

vtk_module_link(OGS::OGSHovmoellerPlotProxy PUBLIC OGS::OGSHovmoellerPlot)
