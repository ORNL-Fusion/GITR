# - Try to find TBB headers and libraries
# Once done this will define
#  TBB_FOUND - System has TBB
#  TBB_INCLUDE_DIRS - The TBB include directories
#  TBB_LIBRARIES - The libraries needed to use TBB
#  TBB_DEFINITIONS - Compiler switches required for using TBB

find_package(PkgConfig)
pkg_check_modules(PC_TBB tbb)
set(TBB_DEFINITIONS ${PC_TBB_CFLAGS_OTHER})

find_path(TBB_INCLUDE_DIR tbb/tbb.h
        HINTS ${PC_TBB_INCLUDEDIR} ${PC_TBB_INCLUDE_DIRS})

find_library(TBB_LIBRARY NAMES tbb tbbmalloc`
        HINTS ${PC_TBB_LIBDIR} ${PC_TBB_LIBRARY_DIRS})

set(TBB_LIBRARIES ${TBB_LIBRARY})
set(TBB_INCLUDE_DIRS ${TBB_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set TBB to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(TBB DEFAULT_MSG
        TBB_LIBRARY TBB_INCLUDE_DIR)

mark_as_advanced(TBB_INCLUDE_DIR TBB_LIBRARY)
