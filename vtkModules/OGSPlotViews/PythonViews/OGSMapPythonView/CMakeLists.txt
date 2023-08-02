set(classes OGSMapPythonView)

vtk_module_add_module(OGS::OGSMapPythonView CLASSES ${classes})
vtk_module_set_property(OGS::OGSMapPythonView PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

