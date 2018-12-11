# - Find FFTW
# Find the native FFTW includes and library
# This module defines
#  FFTW3_INCLUDE_DIR, where to find fftw.h, etc.
#  FFTW3_LIBRARIES, the libraries needed to use FFTW.
#  FFTW3_FOUND, If false, do not try to use FFTW.
# also defined, but not for general use are
#  FFTW3_LIBRARIES, where to find the FFTW library.

# message("FFTW3_DIR set to ${FFTW3_DIR}" )

find_path(FFTW3_INCLUDE_DIR fftw3.h
  ${FFTW3_DIR}/include
  /usr/pkgs64/include
  /usr/include
  /usr/local/include
)

find_library(FFTW3_LIBRARIES
  NAMES fftw3
  PATHS ${FFTW3_DIR}/lib
  "${FFTW3_DIR}\\win33\\lib"
  /usr/lib/x86_64-linux-gnu
  /usr/pkgs64/lib
  /usr/lib64
  /usr/lib
  /usr/local/lib
  NO_DEFAULT_PATH
)

if (FFTW3_LIBRARIES AND FFTW3_INCLUDE_DIR)
  set(FFTW3_LIBRARIES ${FFTW3_LIBRARIES})
  set(FFTW3_FOUND "YES")
  # if (NOT FFTW3_FIND_QUIETLY)
  #   message(STATUS "FindFFTW3 LIBRARIES: ${FFTW3_LIBRARIES}")
  #   message(STATUS "FindFFTW3 INCLUDE_DIR: ${FFTW3_INCLUDE_DIR}")
  # endif()
else()
  set(FFTW3_FOUND "NO")
endif()

if(FFTW3_FOUND)
  # Populate alternative names.
  set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})
  if(NOT FFTW3_LIBRARY)
    set(FFTW3_LIBRARY ${FFTW3_LIBRARIES})
  endif()

  # Add a imported library
  add_library(fftw3 UNKNOWN IMPORTED)
  set_target_properties(fftw3 PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR}")
  set_property(TARGET fftw3 APPEND PROPERTY
    IMPORTED_LOCATION "${FFTW3_LIBRARIES}")

  if (NOT FFTW3_FIND_QUIETLY)
    message(STATUS "Found FFTW: ${FFTW3_LIBRARIES}")
  endif()

else()

  if(FFTW3_FIND_REQUIRED)
    message(FATAL_ERROR "Could not find FFTW library")
  endif()

endif()

# mark_as_advanced(
#   FFTW3_LIBRARIES
#   FFTW3_INCLUDE_DIR
# )
