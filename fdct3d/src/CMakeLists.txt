set(header_files
    commoninc.hpp
    fdct3d.hpp
    fdct3dinline.hpp
    # mexaux.hpp
    nummat.hpp
    numtns.hpp
    numvec.hpp
    offmat.hpp
    offtns.hpp
    offvec.hpp
    )
set(source_files
    fdct3d_forward.cpp
    fdct3d_inverse.cpp
    fdct3d_param.cpp
    )

# Configure
add_library(fdct3d ${source_files})
add_dependencies(fdct3d ep_fftw)
target_link_libraries(fdct3d PUBLIC fftw2)

file(COPY ${header_files} DESTINATION ${PROJECT_BINARY_DIR}/include/fdct3d)
export(TARGETS fdct3d NAMESPACE CurveLab::
    APPEND FILE ${CurveLab_CONFIG_BUILD_DIR}/CurveLabTargets.cmake)

# Install
install(TARGETS fdct3d EXPORT CurveLabTargets DESTINATION ${CurveLab_LIB_INSTALL_DIR})
install(FILES ${header_files} DESTINATION ${CurveLab_INCLUDE_INSTALL_DIR}/fdct3d)

if(CurveLab_BUILD_TESTING)
    add_executable(test_fdct3d test.cpp)
    target_link_libraries(test_fdct3d PUBLIC fdct3d)
    file(COPY options DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    add_test(
        NAME test_fdct3d
        COMMAND test_fdct3d ${CMAKE_CURRENT_BINARY_DIR}/options)
endif()
