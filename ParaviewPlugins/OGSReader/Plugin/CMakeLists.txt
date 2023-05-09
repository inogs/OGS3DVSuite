paraview_add_plugin(OGSReaderPlugin
        VERSION "1.0"
        REQUIRED_ON_SERVER
        SERVER_MANAGER_XML OGSReader.xml
        MODULES OGS::OGSReader
        )


target_link_libraries(OGSReaderPlugin PRIVATE OGS::OGS3DVSuite OGS::OGSReader)
set_target_properties(OGSReaderPlugin PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
