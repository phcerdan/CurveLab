cmake_minimum_required(VERSION 3.5)
project(CurveLab)

option(CurveLab_BUILD_TESTING "BUILD_TESTING for CurveLab" OFF)

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

if(CurveLab_BUILD_TESTING)
    enable_testing()
endif()
#CMake files are broken in FFTW3
#find_package(FFTW3 REQUIRED)
add_subdirectory(fdct3d)
add_subdirectory(fdct_usfft_cpp)
add_subdirectory(fdct_wrapping_cpp)

# WIP, PETSc not working
# Set PETSC_DIR pointing to base directory
# option(CurveLab_USE_MPI "Use MPI, require PETSc" OFF)
if(CurveLab_USE_MPI)
    add_subdirectory(fdct3d_mpi)
endif()
