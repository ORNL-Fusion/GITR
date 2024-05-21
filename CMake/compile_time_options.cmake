# specify C++ standard
set( CMAKE_CXX_STANDARD 20 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

# options are "Debug" and "Release" and "RelWithDebInfo"
set( CMAKE_BUILD_TYPE "Release" )

# preprocessor definitions in source code are defined below:
set( description "(no description added yet)" )
set( GITR_USE_CUDA 1 CACHE STRING "${description}" FORCE )
set( GITR_USE_OPENMP 0 CACHE STRING "${description}" FORCE )
set( GITR_USE_MPI 0 CACHE STRING "${description}" FORCE )
set( GITR_USE_DOUBLE 1 CACHE STRING "${description}" FORCE )

# there are also deprecated options in this list - frozen at a particular value
add_compile_definitions( 
        USE_CUDA=${GITR_USE_CUDA}
        USE_MPI=${GITR_USE_MPI}
        USE_OPENMP=${GITR_USE_OPENMP}
        USE_DOUBLE=${GITR_USE_DOUBLE}
        FIELD_ALIGNED_VALUES=0
        LC_INTERP=0
        GENERATE_LC=0
        BIASED_SURFACE=0
        ODEINT=0 )
