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
include(ExternalProject)

message(STATUS "ExternalProject: fftw...")

ExternalProject_Add(fftw
    URL http://www.fftw.org/fftw-2.1.5.tar.gz
    # URL http://www.fftw.org/fftw-3.3.8.tar.gz
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/fftw
    CONFIGURE_COMMAND
    ${CMAKE_CURRENT_BINARY_DIR}/fftw/src/fftw/configure
    --prefix=${CMAKE_CURRENT_BINARY_DIR}/fftw/install
    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/fftw/install
    )

add_library(FFTW STATIC IMPORTED)
set(lib_fftw_name ${CMAKE_STATIC_LIBRARY_PREFIX}fftw${CMAKE_STATIC_LIBRARY_SUFFIX})
add_dependencies(FFTW fftw)
set(fftw_install_include_dir ${CMAKE_CURRENT_BINARY_DIR}/fftw/install/include)
file(MAKE_DIRECTORY ${fftw_install_include_dir})
set_target_properties(FFTW PROPERTIES
    IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/fftw/install/lib/${lib_fftw_name}
    INTERFACE_INCLUDE_DIRECTORIES ${fftw_install_include_dir}
    )

################# Subdirectories ##################

if(CurveLab_BUILD_TESTING)
    enable_testing()
endif()
#CMake files are broken in FFTW3
#find_package(FFTW3 REQUIRED)
add_subdirectory(fdct3d)
add_subdirectory(fdct_usfft_cpp)
add_subdirectory(fdct_wrapping_cpp)
add_subdirectory(fdct3d_outcore)

# WIP, PETSc not working
# Set PETSC_DIR pointing to base directory
# option(CurveLab_USE_MPI "Use MPI, require PETSc" OFF)
if(CurveLab_USE_MPI)
    add_subdirectory(fdct3d_mpi)
endif()

################# Config ##################
include(CMakePackageConfigHelpers)
configure_package_config_file(
    CurveLabConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/CurveLabConfig.cmake
    INSTALL_DESTINATION ${LIB_INSTALL_DIR}/CurveLab/cmake
    PATH_VARS CurveLab_INCLUDE_INSTALL_DIR)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/CurveLabConfigVersion.cmake
    VERSION ${CurveLab_VERSION}
    COMPATIBILITY SameMajorVersion )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/CurveLabConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/CurveLabConfigVersion.cmake
    DESTINATION ${LIB_INSTALL_DIR}/CurveLab/cmake )
