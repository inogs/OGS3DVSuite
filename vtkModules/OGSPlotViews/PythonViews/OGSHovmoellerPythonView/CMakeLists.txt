set(classes OGSHovmoellerPythonView)

vtk_module_add_module(OGS::OGSHovmoellerPythonView CLASSES ${classes})
vtk_module_set_property(OGS::OGSHovmoellerPythonView PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

