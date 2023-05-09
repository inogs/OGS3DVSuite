set(classes OGSAnnotateDateTime)

if (WITH_PARAVIEW_PLUGINS)
    vtk_module_add_module(OGS::OGSAnnotateDateTime CLASSES ${classes})
    vtk_module_set_property(OGS::OGSAnnotateDateTime PROPERTY INSTALL_RPATH_USE_LINK_PATH VALUE TRUE)
endif ()
