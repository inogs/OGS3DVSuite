set(classes OGSDepthProfile)

vtk_module_add_module(OGS::OGSDepthProfile CLASSES ${classes})
vtk_module_set_property(OGS::OGSDepthProfile PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

vtk_module_link(OGS::OGSDepthProfile PUBLIC OGS::OGS3DVSuite)
