paraview_add_plugin(OGSUtilsPlugin
        VERSION "1.0"
        REQUIRED_ON_SERVER
        SERVER_MANAGER_XML OGSUtils.xml
        MODULES OGS::OGSUtils
        )


target_link_libraries(OGSUtilsPlugin PRIVATE OGS::OGS3DVSuite OGS::OGSUtils)
set_target_properties(OGSUtilsPlugin PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
