set(classes OGSMapPlot)

vtk_module_add_module(OGS::OGSMapPlot CLASSES ${classes})
vtk_module_set_property(OGS::OGSMapPlot PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

