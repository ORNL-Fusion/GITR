# default filepath locations for installed dependency files

if( NOT DEFINED NETCDF_CXX_SHARED_LIB_DIR  )

  set( NETCDF_CXX_SHARED_LIB_DIR  "/usr/lib/x86_64-linux-gnu" "/usr/lib" 
       CACHE STRING "${description}" )
else()

  set( NETCDF_CXX_SHARED_LIB_DIR ${NETCDF_CXX_SHARED_LIB_DIR}
       CACHE STRING "${description}" )
endif()

if( NOT DEFINED NETCDF_CXX_HEADERS_DIR )

  set( NETCDF_CXX_HEADERS_DIR  "/usr/lib/x86_64-linux-gnu" "/usr/lib" 
       CACHE STRING "${description}" )
else()

  set( NETCDF_CXX_HEADERS_DIR  ${NETCDF_CXX_HEADERS_DIR}
       CACHE STRING "${description}" )
endif()

if( NOT DEFINED NETCDF_C_SHARED_LIB_DIR )

set( NETCDF_C_SHARED_LIB_DIR "/usr/lib/x86_64-linux-gnu" "/usr/lib" 
     CACHE STRING "${description}" )

else()

set( NETCDF_C_SHARED_LIB_DIR ${NETCDF_C_SHARED_LIB_DIR}
     CACHE STRING "${description}" )

endif()

if( NOT DEFINED NETCDF_C_HEADERS_DIR )
set( NETCDF_C_HEADERS_DIR "/usr/include" CACHE STRING "${description}" )
else()
set( NETCDF_C_HEADERS_DIR ${NETCDF_C_HEADERS_DIR} CACHE STRING "${description}" )
endif()

if( NOT DEFINED ${LIBCONFIG_CXX_LIB_DIR})
set( LIBCONFIG_CXX_LIB_DIR "/usr/lib/x86_64-linux-gnu" "/usr/lib"
     CACHE STRING "${description}" )
else()
set( LIBCONFIG_CXX_LIB_DIR ${LIBCONFIG_CXX_LIB_DIR}
     CACHE STRING "${description}" )
endif()

if( NOT DEFINED ${LIBCONFIG_C_LIB_DIR} )
set( LIBCONFIG_C_LIB_DIR "/usr/lib/x86_64-linux-gnu" "/usr/lib"
     CACHE STRING "${description}" )
else()
set( LIBCONFIG_C_LIB_DIR ${LIBCONFIG_C_LIB_DIR}
     CACHE STRING "${description}" )
endif()

if( NOT DEFINED ${LIBCONFIG_CXX_HEADERS_DIR} )
set( LIBCONFIG_CXX_HEADERS_DIR "/usr/lib/x86_64-linux-gnu" "/usr/lib"
     CACHE STRING "${description}" )
else()
set( LIBCONFIG_CXX_HEADERS_DIR ${LIBCONFIG_CXX_HEADERS_DIR}
     CACHE STRING "${description}" )
endif()
if( NOT DEFINED ${LIBCONFIG_C_HEADERS_DIR} )
set( LIBCONFIG_C_HEADERS_DIR "/usr/lib/x86_64-linux-gnu" "/usr/lib"
     CACHE STRING "${description}" )
else()
set( LIBCONFIG_C_HEADERS_DIR ${LIBCONFIG_C_HEADERS_DIR}
     CACHE STRING "${description}" )
endif()
