set(classes OGSTimeStatsClimatology)

vtk_module_add_module(OGS::OGSTimeStatsClimatology CLASSES ${classes})
vtk_module_set_property(OGS::OGSTimeStatsClimatology PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

vtk_module_link(OGS::OGSTimeStatsClimatology PUBLIC OGS::OGS3DVSuite)
