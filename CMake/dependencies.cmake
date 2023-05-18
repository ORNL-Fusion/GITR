include( CMake/filepaths_default.cmake )

include( CMake/locate_files.cmake )

# this is a global hack until headers can be found as-needed
include_directories( ${libconfig_c_headers} ${libconfig_cxx_headers} )
include_directories( ${netcdf_c_headers} ${netcdf_cxx_headers} )
