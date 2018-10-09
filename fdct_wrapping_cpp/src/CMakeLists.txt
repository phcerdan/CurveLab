set(source_files
    fdct_wrapping.cpp
    ifdct_wrapping.cpp
    fdct_wrapping_param.cpp
    )

add_library(fdct_wrapping ${source_files})
target_link_libraries(fdct_wrapping FFTW)

if(CurveLab_BUILD_TESTING)
    add_executable(test_fdct_wrapping test.cpp)
    target_link_libraries(test_fdct_wrapping fdct_wrapping)
    file(COPY options DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    add_test(
        NAME test_fdct_wrapping
        COMMAND test_fdct_wrapping ${CMAKE_CURRENT_BINARY_DIR}/options)
endif()
