set(classes OGSSelectOkuboWeiss)

vtk_module_add_module(OGS::OGSSelectOkuboWeiss CLASSES ${classes})
vtk_module_set_property(OGS::OGSSelectOkuboWeiss PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

vtk_module_link(OGS::OGSSelectOkuboWeiss PUBLIC OGS::OGS3DVSuite)
