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
boundary or timing out). The operators acting on the particles are dependent on the prescribed
fields and profiles of the background.


![Trace Impurity Transport](images/TraceImp.png)

![Operator Loop and Equation of Motion](images/GITR_integration.png)


## Directory Layout
```
. Top Level: Top level build system file **CMakeLists.txt**, LICENSE file, README.md (this file)
├── CMake           ---> build system files
├── images          ---> repo visuals for websites and presentations
├── examples        ---> contains a post-processing example script
├── include         ---> C++ header files
├── src             ---> C++ source files
├── test_include    ---> C++ unit test header files
└── test_src        ---> C++ unit test source files
```
## Installation

### Ubuntu 20.04

1. Install g++ and associated utilities. This is best done by installing the "build-essential" package with the command:

> apt install build-essential

This also installs make and other useful utilities.

2. Install the newest version of CMake. This can be done with **apt** or by building from source.
   To build from source, visit the CMake website and download the newest stable release as a
   zipped tar archive. Unizip and extract with

> tar -xvf <cmake_file_name>.tar.gz

   Move the directory into your home directory. Create an out-of-source build directory.
   Navigate to the cmake source directory and run:

> ./bootstrap --prefix=/path/to/cmake_build_directory
> make -j
> make install

  Navigate to the build folder and verify you see bin, doc, and share - bin contains the cmake
  executable. Get the full filepath of the executable in the bin directory with:

> readlink -f <cmake-executable file>

  Running this file like

> ./<cmake-executable-file> --version

  from within the bin folder should print out the version. Now, to make this executable invokable
  from any directory, open the file .bashrc in your home directory or create it if it doesn't
  exist. Add this line to the end:

> alias cmake=<paste the filepath of the cmake executable immediately after the equals sign, no quotes or spaces>

To use this cmake binary, either open a new terminal window to re-parse the .bashrc file, or
re-parse it manually by running:

> source ~/.bashrc

3. To run large problems, you will need to leverage a GPU. This means that nvcc, the NVIDIA CUDA
   compiler, must be installed. Follow these instructions:
   [**CUDA install**](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html)
   You will likely need a restart after this to initiate CUDA drivers.

4. CUDA is now installed, but you must point your shell to it for it to find it. Add the PATH
   and LD\_LIBRARY\_PATH exports to your .bashrc file. Open a new terminal window or source it as
   explained in the CMake installation step for it to take effect.

5. For Netcdf to work properly, you must have the m4 package installed:

> sudo apt install m4

6. For Netcdf to work, you must also have HDF5 installed. Make sure it is installing
   a modern version from the package repo:

> sudo apt install libhdf5-dev

### Mac OSx
1. If you do not have the Homebrew package manager, get it here at: https://brew.sh/ 
2. For `brew` to work, you may need to run the following:
> source ~/.bashrc
3. You must install HDF5:
> brew install hdf5@1.12
4. You must install CMake if you do not already have it:
> brew install cmake
5. You may need to install m4 as well:
> brew install m4

### Configure
Configure build system with CMake. Physics operators can be activated via **-D**-style build-time
 options provided to CMake.

> cmake -S /path/to/GITR -B /path/to/build -D*option_name*

If you get an error that reads: `Compiler requires the CUDA toolkit. Please set the CUDAToolkit_ROOT variable.` use the option flag `-DGITR_USE_CUDA=0`
   
Alternatively:

> cd GITR/build
> cmake -D*option_name* ..

The list of options can be viewed in:

> CMake/user_options.cmake

### Build

Once the project is configured, compile:

> cd build
> make -j

### Run

GITR expects to be run in a directory containing subdirectories **input** and **output**.
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

### Adding a test:

1. Navigate to `CMake/define_test_components.cmake`.

2. Pick a name for the test and add it to
   the CMake variable `gpu_test_targets`
   if it can optionally use the GPU. Otherwise put it in
   `cpu_test_targets`.

3. Create the actual unit test file - it must be named **exactly** the name you picked with .cpp
   at the end, and in the directory `GITR/test_src`. You must include:

> \#include "test\_utils.hpp"

4. Link any libraries your test needs. Do this in `GITR/CMake/crosslink_components.cmake`.

### Adding a data file accessible to the unit tests:
1. Include this file in your test file:

> \#include "test\_data\_filepath.hpp"

   It contains preprocessor definitions for filepaths. This file is automatically generated by
   the build system. To use a data file in the tests, you will need to instruct the build system
   to create a new entry in that header file.

2. Copy your test file into `GITR/test_data/`.

3. Add the following lines anywhere after the macro definition in 
   `GITR/CMake/configure_test_data.cmake`:

> generate\_testing\_file( "test\_data/your\_test\_file.extension" )
> set( YOUR\_PREPROCESSOR\_DEFINE\_NAME ${destination\_path})

4. Navigate to `GITR/CMake/define_define_test_components.cmake`. Add a line:

> \#cmakedefine YOUR\_PREPROCESSOR\_DEFINE\_NAME "${YOUR\_PREPROCESSOR\_DEFINE\_NAME}"

## Bugs/Issues

Create a new issue under GitHub's Issues tab.

## Forum

Search the GitHub discussions tab for existing threads or start a new one.

## Contribute

Fork this repository, branch off of *dev*, and open a merge request into *dev*.

## Release Notes

Navigate to `GITR/docs/release_notes.md`

## Dependencies 

- cmake version 3.13 or newer required
- **CUDA**
  - Enabled by default, disable with -DGITR_USE_CUDA=0
  - Requires existing installation. Set:
    - -DCMAKE_CUDA_COMPILER=/path/to/nvcc
- [**libconfig**](https://github.com/hyperrealm/libconfig)
  - consumes human-readable input config files
  - to use an existing installation, set:
    - -DLIBCONFIG_INCLUDE_DIR=/path/to/libconfig/include (libconfig headers)
    - -DLIBCONFIG_LIBRARY=/path/to/libconfig.so
    - -DLIBCONFIGPP_INCLUDE_DIR=/path/to/libconfig++/include (libconfig headers)
    - -DLIBCONFIGPP_LIBRARY=/path/to/libconfig++.so
- [**netcdf-c**](https://github.com/Unidata/netcdf-c)
  - input/output raw data format
  - to use existing installation set:
    - -DNETCDF_LIBRARY=/path/to/libnetcdf.so
    - -DNETCDF_INCLUDE_DIR=/path/to/netcdf-c/include (netcdf-c headers)
- [**netcdf-cxx4**](https://github.com/Unidata/netcdf-cxx4)
  - C++ extensions to netcdf-c
  - to use an existing installation, set:
    - -DNETCDF_CXX_LIBRARY=/path/to/libnetcdf-cxx4.so
    - -DNETCDF_CXX_INCLUDE_DIR=/path/to/netcdf-cxx4/include (netcdf-c headers)
- [**thrust**](https://github.com/thrust/thrust)
  - header-only library
    - included in CUDA installation if gpu support enabled
  - to use existing installation, set:
    - -DTHRUST_INCLUDE_DIR=/path/to/thrust
    - this should only be necessary if CUDA is disabled

