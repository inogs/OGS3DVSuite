cmake_minimum_required(VERSION 3.24)

project(OGS3DVSuite
        VERSION 0.9
        DESCRIPTION "Utilities for postprocessing of scientific models"
        LANGUAGES C CXX
        )

set(CMAKE_CXX_STANDARD 17)
# There is hardly any good reason for ony enabling -fPIC; moreover,
# this is useful when building static python modules for VTK (but now
# we use shared modules)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# We use the FetchContent mechanism for ensuring that GoogleTest is
# available
include(FetchContent)
# Search for local packages before downloading new ones
set(FETCHCONTENT_TRY_FIND_PACKAGE_MODE "ALWAYS")

# Find external modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Fix install paths
include(GNUInstallDirs)

# Set which files must be compiled
add_library(OGS3DVSuite SHARED
        sources/OGS3DVSuite.cpp
        sources/netcdfio.cpp
        sources/fieldOperations.cpp
        sources/vtkFields.cpp
        sources/vtkOperations.cpp)

set_target_properties(OGS3DVSuite PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
        )


option(WITH_TESTS "Compile the tests for the library" OFF)
option(WITH_VTK_MODULES "Compile the VTK modules (needed for the Paraview plugins)" ON)
option(WITH_VTK_PYTHON_MODULES "Compile python wrapper for the VTK modules" OFF)
option(WITH_PARAVIEW_PLUGINS "Transform the VTK modules into ParaView plugins" ON)

if (WITH_VTK_PYTHON_MODULES)
    if (NOT WITH_VTK_MODULES)
        message(FATAL_ERROR "Python modules require the option \"WITH_VTK_MODULE\" to be ON")
    endif ()
endif ()


if (WITH_PARAVIEW_PLUGINS)
    if (NOT WITH_VTK_MODULES)
        message(FATAL_ERROR "Paraview plugins require the option \"WITH_VTK_MODULE\" to be ON")
    endif ()
endif ()


# Here we search for the external dependencies of this module
message(STATUS "Checking for NetCDF library. If this search fails, "
        "you can give a hint by setting an environment\n   variable "
        "called NETCDF_DIR to the path of your netCDF installation "
        "that contains the netcdf.h file and,\n   if this still fails, "
        "you can set a variable named NETCDF_LIB_DIR that points to "
        "the directory that\n   contains your libnetcdf.so file")
find_package(NetCDF REQUIRED)

find_package(OpenMP REQUIRED)

if (WITH_PARAVIEW_PLUGINS)
    find_package(ParaView REQUIRED)
    add_compile_definitions(WITH_PARAVIEW)
endif ()

find_package(VTK REQUIRED)

if (WITH_TESTS)
    FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG 58d77fa8070e8cec2dc1ed015d66b454c8d78850 # release-1.12.1
    )
    FetchContent_MakeAvailable(googletest)
endif ()

# Compile tests
if (WITH_TESTS)
    add_subdirectory("tests")
endif ()

# Add Vtk Modules
if (WITH_VTK_MODULES)
    add_subdirectory("vtkModules")
endif ()

if (WITH_PARAVIEW_PLUGINS)
    add_subdirectory("ParaviewPlugins")
endif ()

# Set the compilation type, i.e. if we are in "Release" mode or "Debug"
set(default_build_type "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
            STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
            "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()


# Report all warnings and consider them as errors
if (MSVC)
    target_compile_options(OGS3DVSuite PRIVATE /W4 /WX)
else ()
    target_compile_options(OGS3DVSuite PRIVATE -Wall -Wextra -Wpedantic -Werror -Wno-unknown-pragmas)
endif ()


# Set the dependencies for the compilation
target_include_directories(OGS3DVSuite PRIVATE "${NETCDF_INCLUDES}")
target_link_libraries(OGS3DVSuite PRIVATE ${VTK_LIBRARIES} OpenMP::OpenMP_CXX)
target_link_libraries(OGS3DVSuite PRIVATE "${NETCDF_LIBRARIES}")


# Write the RPATH in the library
set_target_properties(OGS3DVSuite PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)


# Specify the headers for the compilation
target_include_directories(
        OGS3DVSuite
        PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        PRIVATE
        ${NETCDF_INCLUDES})


# Fix properties of the previous library
set_target_properties(OGS3DVSuite PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(OGS3DVSuite PROPERTIES SOVERSION 1)


# Set headers that must be exported
set(PUBLIC_HEADER_LIST
        "include/OGS/OGS3DVSuite.h"
        "include/OGS/field.h"
        "include/OGS/fieldOperations.h"
        "include/OGS/geometry.h"
        "include/OGS/macros.h"
        "include/OGS/netcdfio.h"
        "include/OGS/ParaviewImports.h"
        "include/OGS/projection.h"
        "include/OGS/Season.h"
        "include/OGS/TimeInterval.h"
        "include/OGS/TimeList.h"
        "include/OGS/TimeObject.h"
        "include/OGS/TimeRequestors.h"
        "include/OGS/vtkFields.h"
        "include/OGS/V3.h")

set_target_properties(OGS3DVSuite PROPERTIES PUBLIC_HEADER "${PUBLIC_HEADER_LIST}")


# This is useful is this project is imported as a subdir
add_library(OGS::OGS3DVSuite ALIAS OGS3DVSuite)


# Build a pkg file starting from the provided template
configure_file(ogs3dvsuite.pc.in ogs3dvsuite.pc @ONLY)


# Create also a ConfigVersion.cmake file: when other cmake projects link to
# this one, they use this file to check the version of OGS3DVSuite
include(CMakePackageConfigHelpers)
write_basic_package_version_file(OGS3DVSuiteConfigVersion.cmake
        VERSION ${PACKAGE_VERSION}
        COMPATIBILITY AnyNewerVersion)

SET(_targets OGS3DVSuite)
if (WITH_VTK_MODULES)
    SET(_targets ${_targets} OGSReader)
    SET(_targets ${_targets} OGSSelectBasin)
    SET(_targets ${_targets} OGSSelectCoast)
    SET(_targets ${_targets} OGSSelectLand)
    SET(_targets ${_targets} OGSSelectOkuboWeiss)
    SET(_targets ${_targets} OGSSelectPolygon)
    SET(_targets ${_targets} OGSSelectTimePeriod)
    if (WITH_PARAVIEW_PLUGINS)
        SET(_targets ${_targets} OGSAnnotateDateTime)
    endif ()
    SET(_targets ${_targets} OGSSpaghetti)
endif ()

install(TARGETS ${_targets}
        EXPORT OGS3DVSuiteTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/OGS")

# Now we create the *Targets.cmake file and, later, the *Config.cmake file
# that includes the former
install(EXPORT OGS3DVSuiteTargets
        FILE OGS3DVSuiteTargets.cmake
        NAMESPACE OGS::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake)


configure_file(OGS3DVSuiteConfig.cmake.in OGS3DVSuiteConfig.cmake @ONLY)


# Install the auxiliary files
install(FILES ${CMAKE_BINARY_DIR}/ogs3dvsuite.pc
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/OGS3DVSuiteConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/OGS3DVSuiteConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake)


# This is useful if, during development, it is required to link to OGS3DVSuite
# directly from the build directory without performing installation
# export(TARGETS ${_targets} NAMESPACE OGS:: FILE OGS3DVSuiteTargets.cmake)
# set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
# export(PACKAGE OGS3DVSuite)
