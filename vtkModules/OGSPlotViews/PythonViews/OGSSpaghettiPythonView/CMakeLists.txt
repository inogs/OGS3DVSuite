set(classes OGSSpaghettiPythonView)

vtk_module_add_module(OGS::OGSSpaghettiPythonView CLASSES ${classes})
vtk_module_set_property(OGS::OGSSpaghettiPythonView PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)

vtk_module_link(OGS::OGSSpaghettiPythonView PUBLIC OGS::OGSSpaghettiPlotProxy)
