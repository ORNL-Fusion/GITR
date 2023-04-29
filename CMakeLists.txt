# CMakeLists to build GITR for CUDA, OPENMP, CPP

cmake_minimum_required(VERSION 3.1)
project(gitr)

set(CUDA 1) # Set initial argument based on backend desired
               # CUDA, OPENMP, or CPP_PAR (usually for debugging)
set(APPLE 0)   # If using OPENMP or CPP_PAR and running on a mac
set(_MPI 0)    # This is temporarily disabled, leave 0

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMake/)
set(CMAKE_CXX_STANDARD 17) # Doesn't set NVCC C++ standard
set(CXX_STANDARD_REQUIRED TRUE)

# Set source files
file(GLOB SRC_FILES src/*.cpp)
file(GLOB HEADER_FILES include/*.h include/*.hpp test/tests_general.hpp)
file(GLOB HEADER_DIRS include test)

include_directories( ${HEADER_DIRS})

find_package(Thrust REQUIRED)
find_package(NetCDF COMPONENTS CXX REQUIRED)
find_package(LibConfig)
IF (_MPI)
  find_package(MPI REQUIRED)
  IF (MPI_FOUND)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
  endif()
endif()
find_package(OpenMP REQUIRED)

# Add includes
include_directories(${NETCDF_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/include)
include_directories(${THRUST_INCLUDE_DIRS})
include_directories(${LIBCONFIGPP_INCLUDE_DIR})
add_definitions(-DUSE_PARALLEL -DHAS_NETCDF_PAR_H)

# Simulation defines
add_definitions(-DUSE_CUDA=${USE_CUDA}
        -DUSE_MPI=${USE_MPI}
        -DUSE_DOUBLE=${USE_DOUBLE}
        -DUSEFIELDALIGNEDVALUES=0
        -DLC_INTERP=0
        -DGENERATE_LC=0
        -DBIASED_SURFACE=0)

# Setup for CUDA backend
if (${CUDA})
    Message("Enabling CUDA backend")
    find_package(CUDA REQUIRED)
    set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
    set(CUDA_SEPARABLE_COMPILATION TRUE)
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -DCUDA -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA --std=c++17 -O3 --expt-relaxed-constexpr --expt-extended-lambda) #-O3 --expt-extended-lambda --expt-relaxed-constexpr -g -G --cudart shared

    set_source_files_properties(${SRC_FILES} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
    cuda_add_executable(GITR ${SRC_FILES})
endif ()

# setup for OpenMP backend
if(APPLE AND OPENMP)
    if(CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(OpenMP_C "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
      set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
      set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
      set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    endif()
endif()

if (${OPENMP})
    find_package(CUDA REQUIRED)
    Message("Enabling OpenMP backend")
    find_package(OpenMP REQUIRED)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++1z -O3") # -nostdinc++ -I/opt/local/libexec/llvm-devel/include/c++/v1 -L/opt/local/libexec/llvm-devel/lib/c++/v1 -Wl,-rpath,-I/opt/local/libexec/llvm-devel/include/c++/v1") # -Wno-c++11-narrowing -fstack-protector-all -fstack-check -W -Wall") #-O3 -g -G stdlib=libc++
    set_source_files_properties(${SRC_FILES} PROPERTIES LANGUAGE CXX)
    add_executable(GITR ${SRC_FILES})
    add_definitions("-x c++ -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP -DTHRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_OMP")
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
target_link_libraries(GITR ${NETCDF_C_LIBRARIES})
target_link_libraries(GITR ${LIBCONFIGPP_LIBRARIES})
IF (_MPI)
  target_link_libraries(GITR ${MPI_C_LIBRARIES})
  target_link_libraries(GITR ${MPI_CXX_LIBRARIES})
endif()

enable_testing()
#add_test(TestName sft_a)
##add_test(atomic_tests atomic_tests)
#
find_package (Python3 COMPONENTS Interpreter)
message( STATUS "python stuff " ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR} ${EXECUTABLE_OUTPUT_PATH})
add_test(NAME sample_test COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/examples/sft_a/test.py WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples/sft_a )
#add_test(TestName geom_test)
# end test comment

# Install binaries
install(TARGETS GITR DESTINATION bin)
get_directory_property( DirDefs DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS )
foreach( d ${DirDefs} )
    message( STATUS "Found Define: " ${d} )
endforeach()
message( STATUS "DirDefs: " ${DirDefs} )
