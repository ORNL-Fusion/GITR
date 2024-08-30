# GITR 
Global Impurity Transport Code

### Note for legacy GITR users:
For reference, please visit the archived copy of the GITR project - [GITR_legacy](https://github.com/ORNL-Fusion/GITR_legacy)

## Description
The GITR program takes background plasma profiles, equilibrium, geometry, and surface model 
and performs large scale simulation of plasma induced erosion, plasma transport of those 
impurities, sputtering, and redeposition.

The physics implemented in GITR is based on the trace-impurity assumption. i.e. 
The density of impurities created in the GITR code is negligible in its effect on the 
background plasma parameters.

Beginning from a set of initial conditions, the code steps each particle through a set of
operators until certain conditions on the particles are reached (crossing a 
boundary or timing out). The operators acting on the particles are dependent on the prescribed
fields and profiles of the background.


![Trace Impurity Transport](images/TraceImp.png)

![Operator Loop and Equation of Motion](images/GITR_integration.png)

# Installing GITR

GITR is a C++ code that must be compiled on the target platform before it can be run. It can be configured to take advantage of CPU or GPU
parallelism if they are available.

There are 2 options for installation:

1. Container-based installation (recommended)
2. Bare-metal installation (discouraged)

**Note: Users are highly encouraged not to attempt a bare-metal installation unless a container installation is impossible. 
Steps for a bare-metal installation are unique to every platform, and vary depending on software versions throughout
a long stack of system dependencies**

## Container Installation Options

Container images can be built on any device using the Docker Engine and Dockerfile image specifications in the repo. 
Alternatively, pre-configured images can be pulled from the package registry in this repo and run immediately using Singularity, Podman, Docker,
or any other OCI-compilant container service.

Two types of containers exist:

1. Interactive development containers - GITR source code on host is volume mapped into a interactive container environment with all of GITR's dependencies.
   Inside this environment, the user can modify files, compile artifacts, configure hardware usage, and run problems. This is for developers who with to modify GITR.
   These containers can be built on the target systems or pulled pre-built from the package registry.
   
2. Noninteractive run-only containers - GITR binaries are pre-compiled within a non-interactive container image. It is not possible to modify or rebuild GITR in this image.
   This is for users who wish to run GITR on problems, but have no need to modify GITR.

## Building Containers

### container_0 instructions

### container_1 instructions

## Pulling Pre-built containers

### container_0 instructions
### container_1 instructions

## Running Containers:

### interactive container instructions
### noninteractive container instructions

## Running GITR in containers:

### interactive container instructions
### noninteractive container instructions

## Bare Metal Installation Instructions

*Note: content below being migrated into sections above*



## Docker Build (local container, no CUDA)

This is a straightforward option for compiling GITR. This build is recommended for using GITR on a local device. 

1. Install Docker Desktop

2. Open Docker Desktop and log in – this will start the Docker Daemon running

3. From your GITR source directory, run `bash containers/run_alpine.sh`

4. You should be inside a Docker container now, so `cd host` to navigate to the GITR source directory

5. Run the following commands to build and make GITR from inside the container. These commands must be repeated each time you initiate the container.
> rm -rf build
> 
> mkdir build
>
> cmake -S . -B build
>
> cd build
>
> make -j 4

6. Now GITR is built. You can run GITR by executing it from the directory above your input folder while still inside your container. We recommend copying your problem's input directory into `GITR/scratch/input` when using a Docker container. 

> path_to_GITR_build_directory/GITR


## Podman Build (remote container, with CUDA)

This is a straightforward option for compiling GITR. This build is recommended for using GITR on a remote device, such as a NERSC computer. 

1. Check that your remote device has podman by running `podman -v`. If your device does not have podman, we recommend using the Bare Metal Build below.

2. Run `podman pull stonecoldhughes/gitr:env_test` to download a binary container image to provide an environment containing all pre-configured GITR dependencies. GITR will later build and run inside this container. 

3. Open the GITR container from the GITR source directory with `bash containers/run_podman.sh`

4. You should be inside a podman container now, so `cd host` to navigate to the GITR source directory

5. Run the following commands to build and make GITR from inside the container. These commands must be repeated each time you initiate the container.
> rm -rf build
> 
> mkdir build
>
> cmake -S . -B build
>
> cd build
>
> make -j 4

6. Now GITR is built. You can run GITR by executing it from the directory above your input folder while still inside your container. We recommend copying your problem's input directory into `GITR/scratch/input` when using a Docker container. 


## Bare Metal Build (no containers)

This is a straightforward option for compiling GITR.

1. Install libconfig

2. Install netcdf

3. Edit `build.sh` to ensure that it points to the appropriate lib and include **directories** of libconfig and netcdf

4. If you are using a NERSC device or other shared high-performance computing resource:
> module load cmake
> 
> module unload cray-hdf5-parallel
>
> module swap PrgEnv-gnu PrgEnv-nvidia
>
> module load nvidia/23.9
> 
> module load cray-hdf5
> 
> module load cray-netcdf

5. `bash build.sh <path_to_GITR_directory> <path_to_GITR_build_directory>`

6. Check `build/cmake_output.txt` to ensure that "Build files have been written"

7. Check `build/build_output.txt` to ensure that you have "Built target GITR"

8. Now GITR is built. You can run GITR by executing it from the directory above your input folder.

> path_to_GITR_build_directory/GITR

## Environment

The GITR software relies on several other software installations to operate. These *dependencies* fall under *system dependencies* and *3rd-party dependencies*. GITR's advanced build system downloads and builds all 3rd-party dependencies automatically, but this is not the case with system dependencies. These must all be installed by the user prior to attempting to build GITR. Numerous approaches are available to the user for installing/updating/managing system dependencies. An approach using the **spack** utility is briefly described below, loosely following https://spack-tutorial.readthedocs.io/en/latest/.

### Ubuntu 20.04

0. Ensure that basic system dependencies like a working compiler are installed and discoverable on your system. If it is a blank system, you will need to install these with the native Ubuntu package manager *apt*:

> apt install build-essential

At this time, HDF5 must also be installed as a system dependency with the native system package manager. It cannot be installed with spack.

> apt install hdf5-hl

1. Download spack: 

> git clone https://github.com/spack/spack.git

2. Instantiate spack in your environment. This can optionally be placed in your .bashrc file if you want it done automatically upon every login:

> source ~/spack/share/spack/setup-env.sh

3. Direct spack to find a compiler to use:

> spack compilers

This command should produce non-empty output. The discovered compilers will be listed.

4. We may now begin using spack to install the rest of the system dependencies. Beginning with the newest version of gcc:

List packages matching *gcc*:
> spack list gcc

List versions of package:

> spack versions gcc

Install one (preferably the latest):

> spack install gcc @11.2.0

> spack load gcc@11.2.0

> spack compiler find

You are building a literal compiler. Expect this to take a

> Very.

> Long.

> Time...

5. Next, we will use the compiler we just built to build the rest of the dependencies:

> spack list cmake

> spack versions cmake

> spack install cmake @3.22.1 %gcc@11.2.0

6. **Optional**: for CUDA support, similarly install CUDA:

> spack list cuda

> spack versions cuda

> spack install cuda @11.5.1 %gcc@11.2.0

7. **Optional** for blazingly fast source compilation, similarly install Ninja build system:

> spack list ninja

> spack versions ninja

> spack install ninja @1.10.2 %gcc@11.2.0

8. Now that these softwares are made available to spack, they must be loaded into the current environment so that they are discoverable in the current environment. List packages and load them:

> spack find -x

> spack find -x --loaded

> spack load gcc

> spack load cuda

> spack load cmake

> spack load ninja

> spack find -x --loaded

This final command should print out all the loaded environments.

### Mac OSx

1. If you do not have the Homebrew package manager, get it here at: https://brew.sh/ 
2. For `brew` to work, you may need to run the following:
> source ~/.bashrc
3. You must install HDF5:
> brew install hdf5@1.12
4. You must install CMake if you do not already have it:
> brew install cmake
5. You must install Ninja:
> brew install ninja
6. You may need to install m4 as well:
> brew install m4

## Installation
 
### Hardware Configuration:
 
Configure build system with CMake. Physics operators can be activated via **-D**-style build-time
 options provided to CMake.

> cmake -S /path/to/GITR -B /path/to/build -D*option_name*

The list of options can be viewed in:

> CMake/user_options.cmake

### Build

Once the project is configured, compile:

> cd build

If using Unix Makefiles:

> make -j

If using Ninja:

> ninja -j 0

### Run

GITR expects to be run in a directory containing subdirectories **input** and **output**.
The **input** directory must contain a file called *gitrInput.cfg*.

Navigate to this directory and run:

> /path/to/build/GITR

### Configuration options

There are 32 options GITR expects to consume at runtime, in a block of the gitrInput.cfg
file called:

> flags

A list of the required options and a brief description of each is provided:

  - USESURFACEMODEL
    - Binary option - turn on or off for surface modeling
  - FLUX_EA
    - Binary option - turn on to collect energy/angle flux
  - SPECTROSCOPY
    - Ternary option for density histogram collection:
      - 0: do not collect density histograms
      - 2: create 2d histograms
      - 3: create 3d histograms
  - BIASED_SURFACE
    - Binary option, 0 or 1. Enabled by USESHEATHEFIELD
  - USE3DTETGEOM
    - binary variable
      - 0: 3d off, 2d
      - 1: 3d on, 3d
  - USECYLSYMM 
    - description
  - BFIELD_INTERP
    - description
  - GRADT_INTERP
    - description
  - FORCE_EVAL
    - description
  - SORT_PARTICLES
    - description
  - USE_ADAPTIVE_DT
    - description
  - GEOM_HASH
    - description
  - PARTICLE_SOURCE_FILE
    - description
  - PARTICLE_SOURCE_SPACE
    - description
  - PARTICLE_SOURCE_ENERGY
    - description
  - PARTICLE_SOURCE_ANGLE
    - description
  - PARTICLE_TRACKS
    - description
  - PRESHEATH_INTERP
    - description
  - EFIELD_INTERP
    - description
  - USE_SURFACE_POTENTIAL
    - description
  - FLOWV_INTERP
    - description
  - DENSITY_INTERP
    - description
  - TEMP_INTERP
    - description
  - GEOM_HASH_SHEATH
    - description
  - USETHERMALFORCE
    - description
  - USESHEATHEFIELD
    - description
  - USEPRESHEATHEFIELD
    - description
  - USE_IONIZATION
    - description
  - USECOULOMBCOLLISIONS
    - description
  - USEPERPDIFFUSION
    - description
  - USEFIELDALIGNEDVALUES
    - description

## Canonical Example

The default configuration options in GITR are compatible with the input deck in:
[GITR_CPC_example](https://github.com/ORNL-Fusion/GITR-CPC-Example).


## System Dependencies:

- cmake version 3.13 or newer required
 
- gcc
 
- Ninja

- HDF5
 
- **CUDA**
  - Enabled by default, disable with -DGITR_USE_CUDA=0
  - Requires existing installation. Set:
    - -DCMAKE_CUDA_COMPILER=/path/to/nvcc
