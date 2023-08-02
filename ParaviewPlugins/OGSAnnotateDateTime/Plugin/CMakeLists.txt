paraview_add_plugin(OGSAnnotateDateTimePlugin
        VERSION "1.0"
        REQUIRED_ON_SERVER
        SERVER_MANAGER_XML OGSAnnotateDateTime.xml
        MODULES OGS::OGSAnnotateDateTime
)


target_link_libraries(OGSAnnotateDateTimePlugin PRIVATE OGS::OGSAnnotateDateTime)
set_target_properties(OGSAnnotateDateTimePlugin PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
