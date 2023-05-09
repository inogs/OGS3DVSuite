set(REQUIRED_MODULES
        OGS::OGSSelectBasin
        OGS::OGSSelectCoast
        OGS::OGSSelectLand
        OGS::OGSSelectOkuboWeiss
        OGS::OGSSelectPolygon
        OGS::OGSSelectTimePeriod)

paraview_add_plugin(OGSSelectToolsPlugin
        VERSION "1.0"
        REQUIRED_ON_SERVER
        SERVER_MANAGER_XML OGSSelectTools.xml
        MODULES ${REQUIRED_MODULES}
        )


target_link_libraries(OGSSelectToolsPlugin PUBLIC OGS::OGS3DVSuite ${REQUIRED_MODULES})
set_target_properties(OGSSelectToolsPlugin PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
