# GITR
Global Impurity Transport Code

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
├── images          ---> repo visuals for websites and presentations
├── include         ---> C++ header files
├── src             ---> C++ source files
├── test_include    ---> C++ unit test header files
└── test_src        ---> C++ unit test source files
```

## New build instructions:
This branch contains a new build system for GITR. It currently does not support GPUs.
To see the list of preprocessor macros configured and defined by CMake at build time,
examine GITR/CMake/define\_options.cmake. Default values can be overriden on the cmake
command line via -D style options. To build, create a "build\_directory", and run
cmake -S /path/to/GITR -B /path/to/build\_directory -D{whatever option you want to override}
The built files should be in that build\_directory

## Dependencies 

- cmake version 3.13 or newer required
- [libconfig](https://github.com/hyperrealm/libconfig)
  - consumes human-readable input config files
- [netcdf-c](https://github.com/Unidata/netcdf-c)
  - input/output raw data format
- [netcdf-xx4](https://github.com/Unidata/netcdf-cxx4)
  - C++ extensions to netcdf-c
- [thrust](https://github.com/thrust/thrust)
  - header-only library
    - included in CUDA installation if gpu support enabled
    - included explicitly if cpu-only build

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

### Manually-supplied dependencies
If the user has their own pre-installed set of dependencies, point GITR to them by setting the
following specific arguments on the command line:
*Not yet implemented - Put the same list of dependencies but include the variables needed to
point the find() command to it. Find minimum working set*

### Build

Once the project is configured, compile:

> cd build
> make -j

### Run

GITR expects to be run in a directory containing directories **input** and **output**.
The **input** directory must contain a file called *gitrInput.cfg*. Reference 
docs/runtime_config.md for details about this file. These following options in the file must
be mirrored with their CMake **-D**-style counterpart build-time option.

Captain! Here - add a list of the flags. Or just put it in the reference file?

Navigate to this directory and run:

> /path/to/build/GITR


## Testing

Unit test coverage not yet implemented

## Bugs/Issues

Create a new issue under GitHub's Issues tab.

## Forum

Search the GitHub discussions tab for existing threads or start a new one.

## Contribute

Fork this repository, branch off of *dev*, and open a merge request into *dev*.
