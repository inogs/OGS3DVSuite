paraview_add_plugin(OGSSpaghettiPlugin
        VERSION "1.0"
        REQUIRED_ON_SERVER
        SERVER_MANAGER_XML OGSSpaghetti.xml
        MODULES OGS::OGSSpaghetti
        )


target_link_libraries(OGSSpaghettiPlugin PRIVATE OGS::OGS3DVSuite OGS::OGSSpaghetti)
set_target_properties(OGSSpaghettiPlugin PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
