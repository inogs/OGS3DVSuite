find_package(PROJ REQUIRED CONFIG)

set(classes OGSPointSource OGSDepthLineSource)

vtk_module_add_module(OGS::OGSUtils CLASSES ${classes})
vtk_module_set_property(OGS::OGSUtils PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

vtk_module_link(OGS::OGSUtils PUBLIC OGS::OGS3DVSuite)
vtk_module_link(OGS::OGSUtils PRIVATE PROJ::proj)
