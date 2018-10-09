set(source_files
    fdct_usfft.cpp
    afdct_usfft.cpp
    ifdct_usfft.cpp
    fdct_usfft_param.cpp
    )

add_library(fdct_usfft ${source_files})
target_link_libraries(fdct_usfft FFTW)

if(CurveLab_BUILD_TESTING)
    add_executable(test_fdct_usfft test.cpp)
    target_link_libraries(test_fdct_usfft fdct_usfft)
    file(COPY options DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    add_test(
        NAME test_fdct_usfft
        COMMAND test_fdct_usfft ${CMAKE_CURRENT_BINARY_DIR}/options)
endif()
