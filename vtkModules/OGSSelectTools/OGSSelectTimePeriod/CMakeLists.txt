set(classes OGSSelectTimePeriod)

vtk_module_add_module(OGS::OGSSelectTimePeriod CLASSES ${classes})
vtk_module_set_property(OGS::OGSSelectTimePeriod PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

vtk_module_link(OGS::OGSSelectTimePeriod PUBLIC OGS::OGS3DVSuite)
