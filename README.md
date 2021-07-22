# GITR 
Global Impurity Transport Code

[![Build GITR](https://github.com/ORNL-Fusion/GITR/actions/workflows/cmake.yml/badge.svg)](https://github.com/ORNL-Fusion/GITR/actions/workflows/cmake.yml)

### Note for legacy GITR users:
For reference, please visit the archived copy of the GITR project - [GITR_legacy](https://github.com/ORNL-Fusion/GITR_legacy)

## Description
The GITR program takes background plasma profiles, equilibrium, geometry, and surface model 
and performs large scale simulation of plasma induced erosion, plasma transport of those 
impurities, self-sputtering, and redeposition.

The physics implemented in GITR is based on the trace-impurity assumption. i.e. 
The density of impurities created in the GITR code is negligible in its effect on the 
background plasma parameters.

Beginning from a set of initial conditions, the code steps each particle through a set of
operators until certain conditions on the particles are reached (crossing a 
boundary or timing out). The operators acting on the particles are dependent on the prescibed 
fields and profiles of the background.


![Trace Impurity Transport](images/TraceImp.png)

![Operator Loop and Equation of Motion](images/GITR_integration.png)


## Directory Layout
```
. Top Level: Top level build system file **CMakeLists.txt**, LICENSE file, README.md (this file)
├── CMake           ---> build system files
├── docs            ---> more detailed documentation 
├── dev_docs        ---> documentation for developers
├── images          ---> repo visuals for websites and presentations
├── include         ---> C++ header files
├── src             ---> C++ source files
├── test_include    ---> C++ unit test header files
└── test_src        ---> C++ unit test source files
```
## Dependencies 

- cmake version 3.13 or newer required
- **CUDA**
  - Enabled by default, disable with -DGITR_USE_CUDA=0
  - Requires existing installation. Set:
    - -DCMAKE_CUDA_COMPILER=/path/to/nvcc
- [libconfig](https://github.com/hyperrealm/libconfig)
  - consumes human-readable input config files
  - to use existing installation, set:
    - -DLIBCONFIG_INCLUDE_DIR=/path/to/libconfig/include (libconfig headers)
    - -DLIBCONFIG_LIBRARY=/path/to/libconfig.so
    - -DLIBCONFIGPP_INCLUDE_DIR=/path/to/libconfig++/include (libconfig headers)
    - -DLIBCONFIGPP_LIBRARY=/path/to/libconfig++.so
- [netcdf-c](https://github.com/Unidata/netcdf-c)
  - input/output raw data format
  - to use existing installation set:
    - -DNETCDF_LIBRARY=/path/to/libnetcdf.so
    - -DNETCDF_INCLUDE_DIR=/path/to/netcdf-c/include (netcdf-c headers)
- [netcdf-cxx4](https://github.com/Unidata/netcdf-cxx4)
  - C++ extensions to netcdf-c
  - to use an existing installation, set:
    - -DNETCDF_CXX_LIBRARY=/path/to/libnetcdf-cxx4.so
    - -DNETCDF_CXX_INCLUDE_DIR=/path/to/netcdf-cxx4/include (netcdf-c headers)
- [thrust](https://github.com/thrust/thrust)
  - header-only library
    - included in CUDA installation if gpu support enabled
  - to use existing installation, set:
    - -DTHRUST_INCLUDE_DIR=/path/to/thrust
    - this should only be necessary if CUDA is disabled

## Installation

### Configure
Configure build system with CMake. Physics operators can be activated via **-D**-style build-time
 options provided to CMake.

> cmake -S /path/to/GITR -B /path/to/build -D*option_name*

or

> cd GITR/build;

> cmake -D*option_name* ..

The list of options can be viewed in:

> CMake/user_options.cmake

### Build

Once the project is configured, compile:

> cd build

> make -j

### Run

GITR expects to be run in a directory containing directories **input** and **output**.
The **input** directory must contain a file called *gitrInput.cfg*. Reference 
docs/runtime_config.md for details about this file. These following options in the file must
be mirrored with their CMake **-D**-style counterpart build-time option.

Navigate to this directory and run:

> /path/to/build/GITR

## Canonical Example

The default configuration options in GITR are compatible with the input deck in:
[GITR_CPC_example](https://github.com/ORNL-Fusion/GITR-CPC-Example).

## Testing

Navigate to the user-created build directory and run:

> ctest

## Bugs/Issues

Create a new issue under GitHub's Issues tab.

## Forum

Search the GitHub discussions tab for existing threads or start a new one.

## Contribute

Fork this repository, branch off of *dev*, and open a merge request into *dev*.

## Release Notes

Navigate to

> GITR/docs/release_notes.md
