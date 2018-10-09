set(source_files
    cpxnumtnsblkd.cpp
    cpxcrvletprtd.cpp
    fdct3d.cpp
    fdct3d_param.cpp
    fdct3d_forward.cpp
    fdct3d_inverse.cpp
    )

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)
find_package(PETSc REQUIRED)

add_library(fdct3d_mpi ${source_files})
target_link_libraries(fdct3d_mpi FFTW)
target_link_libraries(fdct3d_mpi ${PETSC_LIBRARIES})
target_include_directories(fdct3d_mpi PUBLIC ${PETSC_INCLUDES})

if(CurveLab_BUILD_TESTING)
    add_executable(test_fdct3d_mpi test.cpp)
    target_link_libraries(test_fdct3d_mpi fdct3d)
    file(COPY options DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    add_test(
        NAME test_fdct3d
        COMMAND test_fdct3d_mpi ${CMAKE_CURRENT_BINARY_DIR}/options)
endif()
