set(classes OGSSelectLand)

vtk_module_add_module(OGS::OGSSelectLand CLASSES ${classes})
vtk_module_set_property(OGS::OGSSelectLand PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

vtk_module_link(OGS::OGSSelectLand PUBLIC OGS::OGS3DVSuite)
