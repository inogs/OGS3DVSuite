set(classes OGSSelectCoast)

vtk_module_add_module(OGS::OGSSelectCoast CLASSES ${classes})
vtk_module_set_property(OGS::OGSSelectCoast PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

vtk_module_link(OGS::OGSSelectCoast PUBLIC OGS::OGS3DVSuite)
