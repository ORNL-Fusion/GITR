add_library( libconfig INTERFACE )

target_include_directories( libconfig INTERFACE "/usr/include/" )

target_link_libraries( libconfig INTERFACE "/usr/lib/libconfig++.so" "/usr/lib/libconfig.so" )
