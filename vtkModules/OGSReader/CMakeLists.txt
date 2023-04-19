set(classes OGSReader)

vtk_module_add_module(OGS::OGSReader CLASSES ${classes})
vtk_module_set_property(OGS::OGSReader PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

vtk_module_link(OGS::OGSReader PUBLIC OGS::OGS3DVSuite)
