include(GoogleTest)

set(_target tests)

set(_test_files
        simulation.cpp mesh_vertex_coords.cpp)

configure_file(globals.h.in ${CMAKE_BINARY_DIR}/tests/globals.h)

include_directories(${CMAKE_BINARY_DIR}/tests)

add_executable(${_target} main.cpp ${_test_files})
target_link_libraries(${_target} OGS::OGS3DVSuite gtest gtest_main)

gtest_discover_tests(${_target})
