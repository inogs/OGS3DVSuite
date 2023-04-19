find_package(PROJ REQUIRED CONFIG)

set(classes OGSSelectPolygon)

vtk_module_add_module(OGS::OGSSelectPolygon CLASSES ${classes})
vtk_module_set_property(OGS::OGSSelectPolygon PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

vtk_module_link(OGS::OGSSelectPolygon PUBLIC OGS::OGS3DVSuite)
vtk_module_link(OGS::OGSSelectPolygon PRIVATE PROJ::proj)
