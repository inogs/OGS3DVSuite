set(classes OGSTimeStatsAverage)

vtk_module_add_module(OGS::OGSTimeStatsAverage CLASSES ${classes})
vtk_module_set_property(OGS::OGSTimeStatsAverage PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

vtk_module_link(OGS::OGSTimeStatsAverage PUBLIC OGS::OGS3DVSuite)
