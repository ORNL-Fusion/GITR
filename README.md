<h1 align="center">
  <img src="gitrLogo.jpeg" alt="GITR Logo" width="200">
  <br>
  GITR - Global Impurity Transport Code
</h1>


---

## Overview

GITR is a state-of-the-art simulation program designed for modeling plasma behavior in fusion reactors. It provides comprehensive capabilities for simulating plasma-induced erosion, impurity transport, sputtering, and redeposition. GITR is optimized for performance and operates under the trace-impurity assumption, making it suitable for a wide range of applications in fusion research.

## Directory Structure
```
📦 gitr
┣ 📂 CMake -> Build system files 
┣ 📂 images -> Visuals for websites and presentations
┣ 📂 examples -> Post-processing example script
┣ 📂 include -> C++ header files
┣ 📂 src -> C++ source files
┣ 📂 test_include -> C++ unit test header files
┗ 📂 test_src -> C++ unit test source files
```
## Environment Setup

GITR relies on specific dependencies for operation, including system and 3rd-party dependencies. Here are instructions for setting up the environment on Ubuntu 20.04 and Mac OS X.

### Ubuntu 20.04

1. Install basic system dependencies using apt.
2. Download and configure Spack for managing dependencies.
3. Install required dependencies using Spack.
4. Optionally, install CUDA and Ninja for additional performance.

### Mac OS X

1. Install Homebrew package manager.
2. Install HDF5, CMake, Ninja, and m4 using Homebrew.

## Installation

1. Configure the build system with CMake, specifying build-time options if needed.
2. Compile the project using either Make or Ninja.

## Getting Started

1. Create a directory with subdirectories "input" and "output."
2. Place a "gitrInput.cfg" file in the "input" directory.
3. Run GITR from within this directory.

## Configuration Options

GITR provides various runtime configuration options listed in the "flags" section of "gitrInput.cfg."

## Testing

Use CTest to run tests from the user-created build directory. You can add tests by following the provided instructions.

## Reporting Issues

Did you find a bug or have an issue? Please report it in the [GitHub Issues section](https://github.com/ORNL-Fusion/GITR/issues).

## Contributing

We welcome contributions to GITR! Here's how you can contribute:

1. Fork the repository.
2. Create a branch of "dev."
3. Make your changes and improvements.
4. Open a merge request into "dev."

## Release Notes

For the latest release notes, please refer to [GITR/docs/release_notes.md](https://github.com/ORNL-Fusion/GITR/blob/multi_impurity_ions/release_notes.md).

## System Requirements

- CMake (version 3.13 or newer)
- GCC
- Ninja
- HDF5
- CUDA (optional, enabled by default)

---

<p align="center">
  &copy; 2023 GITR Development Team
</p>
