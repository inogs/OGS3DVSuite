set(classes OGSSpaghetti)

vtk_module_add_module(OGS::OGSSpaghetti CLASSES ${classes})
vtk_module_set_property(OGS::OGSSpaghetti PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

vtk_module_link(OGS::OGSSpaghetti PUBLIC OGS::OGS3DVSuite)
