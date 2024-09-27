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

## Pulling Pre-built containers

If possible, avoid building the container images in favor of pulling pre-built images.

1. Navigate to https://hub.docker.com/r/stonecoldhughes/gitr
2. Navigate to the "tags" page.
3. Copy the "Docker Pull" command and run in the terminal: this will pull the image locally.

## Building Containers

If it is not possible to pull container images from Dockerhub, they can be built locally with Docker.

### Interactive GPU

1. bash containers/build_gpu_interactive.sh

### Interactive CPU

1. bash containers/build_alpine_interactive.sh

### Non-interactive GPU

bash containers/build_gpu_noninteractive.sh

## Running Containers:

### Interactive GPU

bash containers/run_gpu_interactive.sh

This will drop the user into the container environment. The directory /host in the container is the same directory
as the directory where bash/run_gpu_interactive.sh is invoked above. In this example, /host in the container is GITR on the host.


### Interactive CPU

bash containers/run_alpine_interactive.sh

This will drop the user into the container environment. The directory /host in the container is the same directory
as the directory where bash/run_alpine_interactive is invoked above. In this example, /host in the container is GITR on the host.

### Non-interactive GPU

bash containers/run_gpu_noninteractive.sh

note: this should be run in the same location that the "GITR" executable would be run.
The directory must include the input/output directories and the gitrInput.cfg file

## Bare Metal Installation Instructions

A convenience script "build.sh" is included to help build on bare metal. Each -D style CMake option in this file indicates
a filepath to either a library file or a header file. These options must be
modified to point to the correct files and directories on the host system. These file locations vary significantly from one system to the next.

bash build.sh /path/to/GITR/CMakeLists.txt /path/to/build_directory

The paths in this script point to third party libraries and header files GITR compiles with. All of them need to be manually installed:

libconfig
netcdf-c
netcdf-cxx
catch2
cuda
openmp
thrust
hdf5

Please reference any of the Dockerfiles in GITR/containers to see examples of these libraries installed on Linux. The Dockerfiles
in this directory have a *.df file extension.

GITR/containers/gpu_gitr_interactive.df for example.




*Note: content below being migrated into sections above*


## Build Examples:

### Docker Desktop Build (local container, no CUDA)

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

2. Run `podman pull stonecoldhughes/gitr:gpu_gitr_interactive` to download a binary container image to provide an environment containing all pre-configured GITR dependencies. GITR will later build and run inside this container. 

3. Open the GITR container from the GITR source directory with `bash containers/run_podman.sh`

4. You should be inside a podman container now, so `cd /host` to navigate to the GITR source directory

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


## Bare Metal Build on NERSC (no containers)

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


### Mac OSx, No CUDA

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

## Build Process

### Build step 1 (configure the build)
 
Configure build system with CMake. Physics operators can be activated via **-D**-style build-time
 options provided to CMake.

> cmake -S /path/to/GITR -B /path/to/build -D*option_name*

The list of options can be viewed in:

> CMake/user_options.cmake

### Build step 2 (run the build)

Once the project is configured, compile:

> cd build

If using Unix Makefiles:

> make -j

If using Ninja:

> ninja -j 0

### Running the GITR executable

GITR expects to be run in a directory containing subdirectories **input** and **output**.
The **input** directory must contain a file called *gitrInput.cfg*.

Navigate to this directory and run:

> /path/to/build/GITR

## Canonical Example

The default configuration options in GITR are compatible with the input deck in:
[GITR_CPC_example](https://github.com/ORNL-Fusion/GITR-CPC-Example).
