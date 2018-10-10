# - Find FFTW
# Find the native FFTW includes and library
# This module defines
#  FFTW2_INCLUDE_DIR, where to find fftw.h, etc.
#  FFTW2_LIBRARIES, the libraries needed to use FFTW.
#  FFTW2_FOUND, If false, do not try to use FFTW.
# also defined, but not for general use are
#  FFTW2_LIBRARIES, where to find the FFTW library.

# message("FFTW2_DIR set to ${FFTW2_DIR}" )

find_path(FFTW2_INCLUDE_DIR fftw.h
  ${FFTW2_DIR}/include
  /usr/pkgs64/include
  /usr/include
  /usr/local/include
)

find_library(FFTW2_LIBRARIES
  NAMES fftw
  PATHS ${FFTW2_DIR}/lib
  "${FFTW2_DIR}\\win32\\lib"
  /usr/lib/x86_64-linux-gnu
  /usr/pkgs64/lib
  /usr/lib64
  /usr/lib
  /usr/local/lib
  NO_DEFAULT_PATH
)

if (FFTW2_LIBRARIES AND FFTW2_INCLUDE_DIR)
  set(FFTW2_LIBRARIES ${FFTW2_LIBRARIES})
  set(FFTW2_FOUND "YES")
  # if (NOT FFTW2_FIND_QUIETLY)
  #   message(STATUS "FindFFTW2 LIBRARIES: ${FFTW2_LIBRARIES}")
  #   message(STATUS "FindFFTW2 INCLUDE_DIR: ${FFTW2_INCLUDE_DIR}")
  # endif()
else()
  set(FFTW2_FOUND "NO")
endif()

if(FFTW2_FOUND)
  # Populate alternative names.
  set(FFTW2_INCLUDE_DIRS ${FFTW2_INCLUDE_DIR})
  if(NOT FFTW2_LIBRARY)
    set(FFTW2_LIBRARY ${FFTW2_LIBRARIES})
  endif()

  # Add a imported library
  add_library(fftw2 UNKNOWN IMPORTED)
  set_target_properties(fftw2 PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${FFTW2_INCLUDE_DIR}")
  set_property(TARGET fftw2 APPEND PROPERTY
    IMPORTED_LOCATION "${FFTW2_LIBRARIES}")

  if (NOT FFTW2_FIND_QUIETLY)
    message(STATUS "Found FFTW: ${FFTW2_LIBRARIES}")
  endif()

else()

  if(FFTW2_FIND_REQUIRED)
    message(FATAL_ERROR "Could not find FFTW library")
  endif()

endif()

# mark_as_advanced(
#   FFTW2_LIBRARIES
#   FFTW2_INCLUDE_DIR
# )
