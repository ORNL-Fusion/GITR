cmake_minimum_required(VERSION 3.1)
project(gitr)

set(OPENMP 1)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMake/)
set(CMAKE_CXX_STANDARD 14) # Doesn't set NVCC C++ standard
set(CXX_STANDARD_REQUIRED TRUE)

# Set source files
file(GLOB SRC_FILES src/*.cpp)

find_package(Thrust REQUIRED)
find_package(NetCDF COMPONENTS CXX REQUIRED)
find_package(LibConfig)

# Add includes
include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/include)
include_directories(${THRUST_INCLUDE_DIRS})
include_directories(${NETCDF_INCLUDE_DIRS})
include_directories(${LIBCONFIGPP_INCLUDE_DIR})

# Simulation defines
add_definitions(-DUSE_CUDA=0
        -DUSEMPI=0
        -DUSE_BOOST=0
        -DUSEIONIZATION=1
        -DUSERECOMBINATION=0
        -DUSEPERPDIFFUSION=0
        -DUSECOULOMBCOLLISIONS=0
        -DUSETHERMALFORCE=0
        -DUSESURFACEMODEL=0
        -DUSESHEATHEFIELD=1
        -DUSEPRESHEATHFIELD=0
        -DBFIELD_INTERP=2
        -DEFIELD_INTERP=1
        -DPRESHEATH_INTERP=2
        -DDENSITY_INTERP=0
        -DTEMP_INTERP=0
        -DFLOWV_INTERP=0
        -DGRADT_INTERP=2
        -DODEINT=0
        -DFIXEDSEEDS=1
        -DGEOM_TRACE=0
        -DPARTICLE_TRACKS=0
        -DPARTICLE_SOURCE=1
        -DUSECYLSYMM=1)

# Setup for CUDA backend
if (${CUDA})
    Message("Enabling CUDA backend")
    find_package(CUDA REQUIRED)
    set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
    set(CUDA_SEPARABLE_COMPILATION TRUE)
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -DCUDA --expt-extended-lambda --expt-relaxed-constexpr --std=c++11)
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

# Install binaries
install(TARGETS GITR DESTINATION bin)
