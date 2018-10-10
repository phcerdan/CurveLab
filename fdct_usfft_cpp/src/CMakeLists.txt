set(header_files
    fdct_usfft.hpp
    fdct_usfft_inc.hpp
    fdct_usfft_inline.hpp
    # mexaux.hpp
    nummat.hpp
    numvec.hpp
    offmat.hpp
    offvec.hpp
    )
set(source_files
    fdct_usfft.cpp
    afdct_usfft.cpp
    ifdct_usfft.cpp
    fdct_usfft_param.cpp
    )

add_library(fdct_usfft ${source_files})
target_link_libraries(fdct_usfft FFTW)

install(TARGETS fdct_usfft EXPORT CurveLabTargets DESTINATION ${CurveLab_LIB_INSTALL_DIR})
install(FILES ${header_files} DESTINATION ${CurveLab_INCLUDE_INSTALL_DIR}/fdct_usfft)
file(COPY ${header_files} DESTINATION ${PROJECT_BINARY_DIR}/include/fdct_usfft)

if(CurveLab_BUILD_TESTING)
    add_executable(test_fdct_usfft test.cpp)
    target_link_libraries(test_fdct_usfft fdct_usfft)
    file(COPY options DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    add_test(
        NAME test_fdct_usfft
        COMMAND test_fdct_usfft ${CMAKE_CURRENT_BINARY_DIR}/options)
endif()
