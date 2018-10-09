set(source_files
    fdct3d_param.cpp
    fdct3d_forward.cpp
    fdct3d_inverse.cpp
    cpxcrvletocr.cpp
   )

add_library(fdct3d_outcore ${source_files})
target_link_libraries(fdct3d_outcore FFTW)

if(CurveLab_BUILD_TESTING)
    add_executable(test_fdct3d_outcore test.cpp)
    target_link_libraries(test_fdct3d_outcore fdct3d_outcore)
    file(COPY options DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    add_test(
        NAME test_fdct3d_outcore
        COMMAND test_fdct3d_outcore ${CMAKE_CURRENT_BINARY_DIR}/options)
endif()
