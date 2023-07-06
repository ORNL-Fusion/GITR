# Build instructions for conda build

## Creating the conda environment
The environment is created using the `git_env.yaml` file:
```
cd <root of repository>
mamba env create -f git_env.yaml
conda activate gitr_env
```

## Compiling
Make sure that you are on a computing node that has a gpu when compiling.

Building is as simple as:
```
rm -rf build; mkdir build; ./conda_build.sh $PWD/build
```