set(classes OGSVerticalProfilePythonView)

vtk_module_add_module(OGS::OGSVerticalProfilePythonView CLASSES ${classes})
vtk_module_set_property(OGS::OGSVerticalProfilePythonView PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

vtk_module_link(OGS::OGSVerticalProfilePythonView PUBLIC OGS::OGSVerticalProfilePlotProxy)
