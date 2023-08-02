# define test components - CMake "targets" - as separate compilation components
include( CTest ) 

enable_testing()

add_test(NAME sample_test COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/examples/sft_a/test.py WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples/sft_a )
