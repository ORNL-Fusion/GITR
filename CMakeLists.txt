cmake_minimum_required(VERSION 3.1)
project(gitr)

set(CUDA 1)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMake/)
set(CMAKE_CXX_STANDARD 14) # Doesn't set NVCC C++ standard
set(CXX_STANDARD_REQUIRED TRUE)

# Set source files
file(GLOB SRC_FILES src/*.cpp)

find_package(Thrust REQUIRED)
find_package(NetCDF COMPONENTS CXX REQUIRED)
find_package(LibConfig)
find_package (Boost COMPONENTS REQUIRED filesystem timer)
find_package(Boost 1.61.0)
IF (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
    include_directories(${Boost_LIBRARY_DIR})
endif()

# Add includes
include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/include)
include_directories(${THRUST_INCLUDE_DIRS})
include_directories(${NETCDF_INCLUDE_DIRS})
include_directories(${LIBCONFIGPP_INCLUDE_DIR})
# Simulation defines
add_definitions(-DUSE_CUDA=1
        -DUSEMPI=0
        -DUSE_BOOST=1
        -DUSEIONIZATION=1
        -DUSERECOMBINATION=0
        -DUSEPERPDIFFUSION=0
        -DUSECOULOMBCOLLISIONS=0
        -DUSETHERMALFORCE=0
        -DUSESURFACEMODEL=0
        -DUSESHEATHEFIELD=1
        -DBIASED_SURFACE=1
        -DUSEPRESHEATHEFIELD=1
        -DBFIELD_INTERP=0
        -DLC_INTERP=3
        -DGENERATE_LC=1
        -DEFIELD_INTERP=0
        -DPRESHEATH_INTERP=1
        -DDENSITY_INTERP=2
        -DTEMP_INTERP=2
        -DFLOWV_INTERP=1
        -DGRADT_INTERP=2
        -DODEINT=0
        -DFIXEDSEEDS=1
        -DPARTICLESEEDS=1
        -DGEOM_TRACE=0
        -DGEOM_HASH=1
        -DGEOM_HASH_SHEATH=1
        -DPARTICLE_TRACKS=0
        -DPARTICLE_SOURCE=2
        -DSPECTROSCOPY=0
        -DUSE3DTETGEOM=1
        -DUSECYLSYMM=1
        -DCHECK_COMPATIBILITY=1)

# Setup for CUDA backend
if (${CUDA})
    Message("Enabling CUDA backend")
    find_package(CUDA REQUIRED)
    set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
    set(CUDA_SEPARABLE_COMPILATION TRUE)
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -DCUDA  --std=c++11 -g -G --cudart shared) #-O3 --expt-extended-lambda --expt-relaxed-constexpr -g -G --cudart shared
    set_source_files_properties(${SRC_FILES} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
    cuda_add_executable(GITR ${SRC_FILES})
endif ()

# setup for OpenMP backend
if (${OPENMP})
    Message("Enabling OpenMP backend")
    find_package(OpenMP REQUIRED)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set_source_files_properties(${SRC_FILES} PROPERTIES LANGUAGE CXX)
    add_executable(GITR ${SRC_FILES})
    add_definitions("-x c++ -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP")
    set_target_properties(GITR PROPERTIES LINKER_LANGUAGE CXX)
endif ()

# Setup for CPP backend
if (${CPP_PAR})
    Message("Enabling CPP backend")
    add_executable(GITR ${SRC_FILES})
    add_definitions("-x c++ -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP -DCPP_PAR -Wno-unused-local-typedef")
endif ()

# Add debug libraries
if (${CUDA_DEBUG})
    target_link_libraries(GITR ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libnvToolsExt.so)
endif ()

# Add libraries
target_link_libraries(GITR ${NETCDF_CXX_LIBRARIES})
target_link_libraries(GITR ${LIBCONFIGPP_LIBRARIES})
target_link_libraries(GITR Boost::system Boost::filesystem Boost::timer)
    #libboost_system};${lboost_timer};${lboost_filesystem})
# Install binaries
install(TARGETS GITR DESTINATION bin)
