set(classes OGSTimeStatsAggregator)

vtk_module_add_module(OGS::OGSTimeStatsAggregator CLASSES ${classes})
vtk_module_set_property(OGS::OGSTimeStatsAggregator PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

vtk_module_link(OGS::OGSTimeStatsAggregator PUBLIC OGS::OGS3DVSuite)
