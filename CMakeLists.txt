cmake_minimum_required(VERSION 3.1)
project(gitr)

set(OPENMP 1)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMake/)
set(CMAKE_CXX_STANDARD 17) # Doesn't set NVCC C++ standard
set(CXX_STANDARD_REQUIRED TRUE)

# Set source files
file(GLOB SRC_FILES src/*.cpp)
file(GLOB TEST_FILES src/*.cpp test/tests_general.cpp)
file(GLOB TEST_FILE_IO test/file_io_tests.cpp)
file(GLOB TEST_COULOMB test/coulomb_tests.cpp)
file(GLOB TEST_FIELD test/field_tests.cpp)
#file(GLOB TEST_ATOMIC test/atomic_tests.cpp)
list (REMOVE_ITEM TEST_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/gitr.cpp)
message("${TEST_FILES}")
file(GLOB HEADER_FILES include/*.h include/*.hpp test/tests_general.hpp)
add_library(YourLib ${TEST_FILES})
include_directories( ${HEADER_FILES})
find_package(Thrust REQUIRED)
find_package(NetCDF COMPONENTS CXX REQUIRED)
find_package(LibConfig)
find_package(MPI)
IF (MPI_FOUND)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
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
        -DUSEIONIZATION=${USEIONIZATION}
        -DUSERECOMBINATION=${USERECOMBINATION}
        -DUSEPERPDIFFUSION=${USEPERPDIFFUSION}
        -DUSECOULOMBCOLLISIONS=${USECOULOMBCOLLISIONS}
        -DUSEFRICTION=${USEFRICTION}
        -DUSEANGLESCATTERING=${USEANGLESCATTERING}
        -DUSEHEATING=${USEHEATING}
        -DUSETHERMALFORCE=${USETHERMALFORCE}
        -DUSESURFACEMODEL=${USESURFACEMODEL}
        -DUSESHEATHEFIELD=${USESHEATHEFIELD}
        -DBIASED_SURFACE=${BIASED_SURFACE}
        -DUSEPRESHEATHEFIELD=${USEPRESHEATHEFIELD}
        -DBFIELD_INTERP=${BFIELD_INTERP}
        -DLC_INTERP=${LC_INTERP}
        -DGENERATE_LC=${GENERATE_LC}
        -DEFIELD_INTERP=${EFIELD_INTERP}
        -DPRESHEATH_INTERP=${PRESHEATH_INTERP}
        -DDENSITY_INTERP=${DENSITY_INTERP}
        -DTEMP_INTERP=${TEMP_INTERP}
        -DFLOWV_INTERP=${FLOWV_INTERP}
        -DGRADT_INTERP=${GRADT_INTERP}
        -DODEINT=${ODEINT}
        -DFIXEDSEEDS=${FIXEDSEEDS}
        -DPARTICLESEEDS=${PARTICLESEEDS}
        -DGEOM_TRACE=${GEOM_TRACE}
        -DGEOM_HASH=${GEOM_HASH}
        -DGEOM_HASH_SHEATH=${GEOM_HASH_SHEATH}
        -DPARTICLE_TRACKS=${PARTICLE_TRACKS}
        -DPARTICLE_SOURCE=${PARTICLE_SOURCE}
        -DPARTICLE_SOURCE_SPACE=${PARTICLE_SOURCE_SPACE}
        -DPARTICLE_SOURCE_ENERGY=${PARTICLE_SOURCE_ENERGY}
        -DPARTICLE_SOURCE_ANGLE=${PARTICLE_SOURCE_ANGLE}
        -DPARTICLE_SOURCE_FILE=${PARTICLE_SOURCE_FILE}
        -DSPECTROSCOPY=${SPECTROSCOPY}
        -DUSE3DTETGEOM=${USE3DTETGEOM}
        -DUSECYLSYMM=${USECYLSYMM}
        -DUSEFIELDALIGNEDVALUES=${USEFIELDALIGNEDVALUES}
        -DFLUX_EA=${FLUX_EA}
        -DFORCE_EVAL=${FORCE_EVAL}
        -DUSE_SORT=${USE_SORT}
        -DCHECK_COMPATIBILITY=${CHECK_COMPATIBILITY})

# Setup for CUDA backend
if (${CUDA})
    Message("Enabling CUDA backend")
    find_package(CUDA REQUIRED)
    set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
    set(CUDA_SEPARABLE_COMPILATION TRUE)
    #set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -DCUDA -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA --std=c++14 -O3 --expt-relaxed-constexpr) #-O3 --expt-extended-lambda --expt-relaxed-constexpr -g -G --cudart shared
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -DCUDA -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA --std=c++14 -O3 --expt-relaxed-constexpr --expt-extended-lambda) #-O3 --expt-extended-lambda --expt-relaxed-constexpr -g -G --cudart shared

    set_source_files_properties(${SRC_FILES} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
    #test comment
    set_source_files_properties(${TEST_FILES} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
    set_source_files_properties(${TEST_FILE_IO} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
    set_source_files_properties(${TEST_COULOMB} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
    set_source_files_properties(${TEST_ATOMIC} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
    set_source_files_properties(${TEST_FIELD} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
    CUDA_ADD_EXECUTABLE(tests ${TEST_FILE_IO} ${TEST_FILES})
    CUDA_ADD_EXECUTABLE(coulomb_tests ${TEST_COULOMB} ${TEST_FILES})
    CUDA_ADD_EXECUTABLE(field_tests ${TEST_FIELD} ${TEST_FILES})
    CUDA_ADD_EXECUTABLE(atomic_tests ${TEST_ATOMIC} ${TEST_FILES})
    cuda_add_executable(GITR ${SRC_FILES})
endif ()

# setup for OpenMP backend
if (${OPENMP})
    ADD_EXECUTABLE(tests test/file_io_tests.cpp)
    ADD_EXECUTABLE(coulomb_tests test/coulomb_tests.cpp)
    ADD_EXECUTABLE(field_tests test/field_tests.cpp)
    #ADD_EXECUTABLE(atomic_tests test/atomic_tests.cpp)
    #ADD_EXECUTABLE(geom_test test/geometry_tests.cpp)
    Message("Enabling OpenMP backend")
    find_package(OpenMP REQUIRED)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O3 -gline-tables-only -std=c++1z -std=c++1z -stdlib=libc++ -nostdinc++ -I/opt/local/libexec/llvm-devel/include/c++/v1 -L/opt/local/libexec/llvm-devel/lib/c++/v1 -Wl,-rpath,-I/opt/local/libexec/llvm-devel/include/c++/v1") # -Wno-c++11-narrowing -fstack-protector-all -fstack-check -W -Wall") #-O3 -g -G -fsanitize=address
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++17") # -nostdinc++ -I/opt/local/libexec/llvm-devel/include/c++/v1 -L/opt/local/libexec/llvm-devel/lib/c++/v1 -Wl,-rpath,-I/opt/local/libexec/llvm-devel/include/c++/v1") # -Wno-c++11-narrowing -fstack-protector-all -fstack-check -W -Wall") #-O3 -g -G stdlib=libc++
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
target_link_libraries(GITR ${MPI_C_LIBRARIES})
target_link_libraries(GITR ${MPI_CXX_LIBRARIES})
target_link_libraries(GITR YourLib)
#target_link_libraries(GITR "-lstdc++fs")
#test comment
target_link_libraries(tests YourLib)
target_link_libraries(tests ${NETCDF_CXX_LIBRARIES})
target_link_libraries(tests ${NETCDF_C_LIBRARIES})
target_link_libraries(tests ${LIBCONFIGPP_LIBRARIES})
target_link_libraries(tests ${MPI_C_LIBRARIES})
target_link_libraries(tests ${MPI_CXX_LIBRARIES})
#target_link_libraries(tests "-lstdc++fs")

#target_link_libraries(geom_test YourLib)
#target_link_libraries(geom_test ${NETCDF_CXX_LIBRARIES})
#target_link_libraries(geom_test ${NETCDF_C_LIBRARIES})
#target_link_libraries(geom_test ${LIBCONFIGPP_LIBRARIES})
#target_link_libraries(geom_test ${MPI_C_LIBRARIES})
#target_link_libraries(geom_test ${MPI_CXX_LIBRARIES})
##target_link_libraries(tests "-lstdc++fs")
target_link_libraries(coulomb_tests YourLib)
target_link_libraries(coulomb_tests ${NETCDF_CXX_LIBRARIES})
target_link_libraries(coulomb_tests ${NETCDF_C_LIBRARIES})
target_link_libraries(coulomb_tests ${LIBCONFIGPP_LIBRARIES})
target_link_libraries(coulomb_tests ${MPI_C_LIBRARIES})
target_link_libraries(coulomb_tests ${MPI_CXX_LIBRARIES})
#target_link_libraries(coulomb_tests "-lstdc++fs")
target_link_libraries(field_tests YourLib)
target_link_libraries(field_tests ${NETCDF_CXX_LIBRARIES})
target_link_libraries(field_tests ${NETCDF_C_LIBRARIES})
target_link_libraries(field_tests ${LIBCONFIGPP_LIBRARIES})
target_link_libraries(field_tests ${MPI_C_LIBRARIES})
target_link_libraries(field_tests ${MPI_CXX_LIBRARIES})
enable_testing()
add_test(TestName tests)
add_test(Coulomb_tests coulomb_tests)
add_test(field_tests field_tests)
#target_link_libraries(atomic_tests YourLib)
#target_link_libraries(atomic_tests ${NETCDF_CXX_LIBRARIES})
#target_link_libraries(atomic_tests ${NETCDF_C_LIBRARIES})
#target_link_libraries(atomic_tests ${LIBCONFIGPP_LIBRARIES})
#target_link_libraries(atomic_tests ${MPI_C_LIBRARIES})
#target_link_libraries(atomic_tests ${MPI_CXX_LIBRARIES})
###target_link_libraries(atomic_tests "-lstdc++fs")
enable_testing()
add_test(TestName tests)
#add_test(atomic_tests atomic_tests)

find_package (Python3 COMPONENTS Interpreter)
message( STATUS "python stuff " ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR} ${EXECUTABLE_OUTPUT_PATH})
add_test(NAME sample_test COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/sample_test.py WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test )
#add_test(TestName geom_test)
# end test comment
# Install binaries
install(TARGETS GITR DESTINATION bin)
get_directory_property( DirDefs DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS )
foreach( d ${DirDefs} )
    message( STATUS "Found Define: " ${d} )
endforeach()
message( STATUS "DirDefs: " ${DirDefs} )
