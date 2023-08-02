set(REQUIRED_MODULES
        OGS::OGSTimeStatsAggregator
        OGS::OGSTimeStatsAverage
        OGS::OGSTimeStatsClimatology)

paraview_add_plugin(OGSTimeStatsPlugin
        VERSION "1.0"
        REQUIRED_ON_SERVER
        SERVER_MANAGER_XML OGSTimeStats.xml
        MODULES ${REQUIRED_MODULES}
)


target_link_libraries(OGSTimeStatsPlugin PRIVATE OGS::OGS3DVSuite ${REQUIRED_MODULES})
set_target_properties(OGSTimeStatsPlugin PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
