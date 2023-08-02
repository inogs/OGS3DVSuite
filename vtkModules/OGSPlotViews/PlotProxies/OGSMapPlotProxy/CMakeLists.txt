set(classes OGSMapPlotProxy)

vtk_module_add_module(OGS::OGSMapPlotProxy CLASSES ${classes})
vtk_module_set_property(OGS::OGSMapPlotProxy PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

vtk_module_link(OGS::OGSMapPlotProxy PUBLIC OGS::OGSMapPlot)
