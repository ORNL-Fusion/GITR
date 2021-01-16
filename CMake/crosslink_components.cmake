# link together singleton dependencies
# Here are the source targets:
#     efield_interp - needs no linking?
#     interp2d - needs no linking
#     particle
#     utils
#     flags
#     setup
#     GITR

# Here are the test targets:
#     atomic_tests
#     coulomb_tests
#     field_tests
#     file_io_tests

# Here are the external interface library targets:
# thrust
# libconfig
# netcdf

# notes above, code below
# link source targets
target_link_libraries( interp2d thrust )
target_link_libraries( flags libconfig thrust )
target_link_libraries( utils libconfig thrust )

# link test targets
target_link_libraries( coulomb_tests test_utils libconfig thrust )
target_link_libraries( atomic_tests test_utils )
target_link_libraries( field_tests test_utils )
target_link_libraries( file_io_tests test_utils )

