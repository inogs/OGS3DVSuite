set(classes OGSVerticalProfilePlotProxy)

vtk_module_add_module(OGS::OGSVerticalProfilePlotProxy CLASSES ${classes})
vtk_module_set_property(OGS::OGSVerticalProfilePlotProxy PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

vtk_module_link(OGS::OGSVerticalProfilePlotProxy PUBLIC OGS::OGSVerticalProfilePlot)
