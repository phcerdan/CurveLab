cmake_minimum_required(VERSION 3.5)
project(CurveLab)

################# VERSION ##################
set(CurveLab_VERSION_MAJOR 2 CACHE STRING "major version" FORCE)
set(CurveLab_VERSION_MINOR 1 CACHE STRING "minor version" FORCE)
set(CurveLab_VERSION_PATCH 2 CACHE STRING "patch version" FORCE)
set(CurveLab_VERSION ${CurveLab_VERSION_MAJOR}.${CurveLab_VERSION_MINOR}.${CurveLab_VERSION_PATCH} CACHE STRING "version" FORCE)

################# INSTALL_DIR ##################
set(CurveLab_INCLUDE_INSTALL_DIR include/ CACHE PATH "INCLUDE_INSTALL_DIR")
set(CurveLab_LIB_INSTALL_DIR lib/ CACHE PATH "LIB_INSTALL_DIR")

option(CurveLab_BUILD_TESTING "BUILD_TESTING for CurveLab" OFF)

################# ExternalProject ##################
#CMake files are broken in FFTW3
# find_package(FFTW3 REQUIRED)

message(STATUS "ExternalProject: fftw...")
include(ExternalProject)
ExternalProject_Add(ep_fftw
    # URL http://www.fftw.org/fftw-2.1.5.tar.gz
    URL http://www.fftw.org/fftw-3.3.8.tar.gz
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/fftw
    CONFIGURE_COMMAND
    ${CMAKE_CURRENT_BINARY_DIR}/fftw/src/ep_fftw/configure
    --prefix=${CMAKE_CURRENT_BINARY_DIR}/fftw/install
    --enable-shared=yes
    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/fftw/install
    )

##### IGNORE
# add_library(fftw2 STATIC IMPORTED)
# set(lib_fftw_name ${CMAKE_STATIC_LIBRARY_PREFIX}fftw${CMAKE_STATIC_LIBRARY_SUFFIX})
# add_dependencies(fftw ep_fftw)

# set(fftw_install_include_dir ${CMAKE_CURRENT_BINARY_DIR}/fftw/install/include)
# file(MAKE_DIRECTORY ${fftw_install_include_dir})
# set_target_properties(fftw PROPERTIES
#     IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/fftw/install/lib/${lib_fftw_name}
#     INTERFACE_INCLUDE_DIRECTORIES $<BUILD_INTERFACE:${fftw_install_include_dir}>
#     # INTERFACE_INCLUDE_DIRECTORIES $<INSTALL_INTERFACE:${CurveLab_INCLUDE_INSTALL_DIR}>
#     )
# install(TARGETS fftw EXPORT CurveLabTargets DESTINATION ${CurveLab_LIB_INSTALL_DIR})
# install(FILES ${fftw_install_include_dir}/fftw.h DESTINATION ${CurveLab_INCLUDE_INSTALL_DIR})

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
if("${FFTW3_DIR}" STREQUAL "")
    set(FFTW3_DIR ${CMAKE_CURRENT_BINARY_DIR}/fftw/install)
endif()
# TODO Fix: This requires double run of cmake configuration. Superbuild is not handled well.
find_package(FFTW3 MODULE) # This is required, but we leave it optional to build the ExternalProject
message(WARNING "FFTW3_DIR: ${FFTW3_DIR}")
message(STATUS "FFTW3_INCLUDE_DIR: ${FFTW3_INCLUDE_DIR}")
message(STATUS "FFTW3_LIBRARIES: ${FFTW3_LIBRARIES}")

################# Subdirectories ##################

set(CurveLab_CONFIG_INSTALL_DIR ${CurveLab_LIB_INSTALL_DIR}/CurveLab/cmake)
set(CurveLab_CONFIG_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})
# Clear Targets.cmake (exported targets are appended to this file from each subdirectory)
file(WRITE ${CurveLab_CONFIG_BUILD_DIR}/CurveLabTargets.cmake "")

if(CurveLab_BUILD_TESTING)
    enable_testing()
endif()
add_subdirectory(fdct3d)
if(_disable)
add_subdirectory(fdct_usfft_cpp)
add_subdirectory(fdct_wrapping_cpp)
add_subdirectory(fdct3d_outcore)
endif()

# WIP, PETSc not working
# Set PETSC_DIR pointing to base directory
# option(CurveLab_USE_MPI "Use MPI, require PETSc" OFF)
if(CurveLab_USE_MPI)
    add_subdirectory(fdct3d_mpi)
endif()

################# Config ##################
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/CurveLabConfigVersion.cmake
    VERSION ${CurveLab_VERSION}
    COMPATIBILITY SameMajorVersion )

# Config for build time
set(CurveLab_INCLUDE_BUILD_DIR ${PROJECT_BINARY_DIR}/include)
set(CurveLab_INCLUDE_DIR ${CurveLab_INCLUDE_BUILD_DIR})
set(CurveLab_CONFIG_DIR  ${CurveLab_CONFIG_BUILD_DIR})
set(CurveLab_CMAKE_DIR   ${CMAKE_CURRENT_BINARY_DIR})
configure_package_config_file(
    cmake/CurveLabConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/CurveLabConfig.cmake
    INSTALL_DESTINATION ${PROJECT_BINARY_DIR}
    PATH_VARS
    CurveLab_INCLUDE_DIR
    CurveLab_CONFIG_DIR
    CurveLab_CMAKE_DIR)

# Configure at install time. Located out of the way from the build tree.
set(CurveLab_INCLUDE_DIR ${CurveLab_INCLUDE_INSTALL_DIR})
set(CurveLab_CONFIG_DIR ${CurveLab_CONFIG_INSTALL_DIR})
configure_package_config_file(
    cmake/CurveLabConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CurveLabConfig.cmake
    INSTALL_DESTINATION ${CurveLab_CONFIG_INSTALL_DIR}
    PATH_VARS CurveLab_INCLUDE_DIR CurveLab_CONFIG_DIR)

################# Install ##################
install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CurveLabConfig.cmake
    DESTINATION ${CurveLab_CONFIG_INSTALL_DIR})
install(
    FILES cmake/FindFFTW3.cmake
    DESTINATION ${CurveLab_CONFIG_INSTALL_DIR})
file(COPY cmake/FindFFTW3.cmake
    DESTINATION ${CurveLab_CONFIG_BUILD_DIR})

install(EXPORT CurveLabTargets
    NAMESPACE CurveLab::
    DESTINATION ${CurveLab_CONFIG_INSTALL_DIR})
