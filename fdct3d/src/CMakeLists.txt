set(source_files
    fdct3d_forward.cpp
    fdct3d_inverse.cpp
    fdct3d_param.cpp
    )

add_library(fdct3d ${source_files})
target_link_libraries(fdct3d FFTW)

if(CurveLab_BUILD_TESTING)
    add_executable(test_fdct3d test.cpp)
    target_link_libraries(test_fdct3d fdct3d)
    file(COPY options DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    add_test(
        NAME test_fdct3d
        COMMAND test_fdct3d ${CMAKE_CURRENT_BINARY_DIR}/options)
endif()
