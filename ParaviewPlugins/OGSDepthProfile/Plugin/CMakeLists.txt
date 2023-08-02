paraview_add_plugin(OGSDepthProfilePlugin
        VERSION "1.0"
        REQUIRED_ON_SERVER
        SERVER_MANAGER_XML OGSDepthProfile.xml
        MODULES OGS::OGSDepthProfile
)


target_link_libraries(OGSDepthProfilePlugin PRIVATE OGS::OGS3DVSuite OGS::OGSDepthProfile)
set_target_properties(OGSDepthProfilePlugin PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
