set(classes OGSSelectBasin)

vtk_module_add_module(OGS::OGSSelectBasin CLASSES ${classes})
vtk_module_set_property(OGS::OGSSelectBasin PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

vtk_module_link(OGS::OGSSelectBasin PUBLIC OGS::OGS3DVSuite)
