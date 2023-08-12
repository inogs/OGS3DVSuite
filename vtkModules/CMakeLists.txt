set(_rejected_modules "")

SET(PythonViewModules
        OGS::OGSHovmoellerPythonView
        OGS::OGSMapPythonView
        OGS::OGSSpaghettiPythonView
        OGS::OGSVerticalProfilePythonView)

SET(PythonPlotModules
        OGS::OGSHovmoellerPlot
        OGS::OGSMapPlot
        OGS::OGSSpaghettiPlot
        OGS::OGSVerticalProfilePlot)

SET(PlotProxyModules
        OGS::OGSHovmoellerPlotProxy
        OGS::OGSMapPlotProxy
        OGS::OGSSpaghettiPlotProxy
        OGS::OGSVerticalProfilePlotProxy)


if (NOT WITH_PARAVIEW_PLUGINS)
    list(APPEND _reject_modules "OGS::OGSAnnotateDateTime")
    foreach(module IN ITEMS ${PythonViewModules})
        list(APPEND _reject_modules ${module})
    endforeach()
    foreach(module IN ITEMS ${PythonPlotModules})
        list(APPEND _reject_modules ${module})
    endforeach()
    foreach(module IN ITEMS ${PlotProxyModules})
        list(APPEND _reject_modules ${module})
    endforeach()
endif ()

vtk_module_find_modules(module_files ${CMAKE_CURRENT_SOURCE_DIR})

vtk_module_scan(
        MODULE_FILES ${module_files}
        PROVIDES_MODULES OGSModules
        WANT_BY_DEFAULT ON
        REJECT_MODULES ${_reject_modules}
)

set(clean_modules "")
foreach(module IN ITEMS ${OGSModules})
    string(REPLACE "OGS::" "" module_name ${module})
    list(APPEND clean_modules ${module_name})
endforeach()
set(OGS3DVSuiteModules ${clean_modules} PARENT_SCOPE)


# I do not know why, but without this OGSReader is compiled as a static library
set(BUILD_SHARED_LIBS ON)

vtk_module_build(MODULES ${OGSModules} BUILD_WITH_KITS OFF ENABLE_WRAPPING ON)

if (WITH_VTK_PYTHON_MODULES)
    vtk_module_python_default_destination(python_destination)
    vtk_module_wrap_python(MODULES ${OGSModules}
            TARGET OGS::OGS3DVSuitePythonModules
            WRAPPED_MODULES wrapped
            PYTHON_PACKAGE "vtkOGS"
            MODULE_DESTINATION "${python_destination}"
            CMAKE_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/"
            LIBRARY_DESTINATION "${CMAKE_INSTALL_LIBDIR}"
            BUILD_STATIC OFF
            INSTALL_HEADERS OFF)
endif ()
