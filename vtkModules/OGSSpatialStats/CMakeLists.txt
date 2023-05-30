set(classes OGSSpatialStats OGSSpatialStatsFromFile)

vtk_module_add_module(OGS::OGSSpatialStats CLASSES ${classes})
vtk_module_set_property(OGS::OGSSpatialStats PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

vtk_module_link(OGS::OGSSpatialStats PUBLIC OGS::OGS3DVSuite)
