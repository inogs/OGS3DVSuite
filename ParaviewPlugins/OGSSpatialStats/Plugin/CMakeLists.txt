paraview_add_plugin(OGSSpatialStatsPlugin
        VERSION "1.0"
        REQUIRED_ON_SERVER
        SERVER_MANAGER_XML OGSSpatialStats.xml
        MODULES OGS::OGSSpatialStats
)


target_link_libraries(OGSSpatialStatsPlugin PRIVATE OGS::OGS3DVSuite OGS::OGSSpatialStats)
set_target_properties(OGSSpatialStatsPlugin PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
