# - Try to find Thrust headers
# Once done this will define
#  THRUST_FOUND - System has thrust
#  THRUST_INCLUDE_DIRS - The thrust include directories
#  THRUST_DEFINITIONS - Compiler switches required for using thrust

find_path(THRUST_INCLUDE_DIR thrust/version.h
        HINTS /usr/include/cuda /usr/local/include /usr/local/cuda/include /opt/cuda/include ../external/thrust)

set(THRUST_INCLUDE_DIRS ${THRUST_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set THRUST to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(Thrust DEFAULT_MSG
        THRUST_INCLUDE_DIRS)

mark_as_advanced(THRUST_INCLUDE_DIRS)
