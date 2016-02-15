# GITR
Global Impurity Transport Code

GITR is currently in the prototype phase of development for large scale simulation of plasma induced erosion, plasma transport of those impurities, self-sputtering, and redeposition.

The physics implemented in GITR is based on the trace-impurity assumption. i.e. The density of impurities created in the GITR code is negligible in its effect on the background plasma parameters.

![Trace Impurity Transport](TraceImp.png)

The code is manifested in the time stepping of a particle with a set of initial conditions through a set of operators until certain conditions on the particles are reached (crossing a boundary or timing out). The operators acting on the particles are dependent on the prescibed fields and profiles of the background.

![Operator Loop and Equation of Motion](PlasmaTransp.png)

Inputs to the code are set up in the gitrInput.m file and include
1. Volume specification and gridding
2. Surface parameterization and gridding
3. Background fields and profiles
4. Number of particles and initial conditions
5. Time stepping information
6. Desired output

There also exists a folder for processing output data /matlab/PostProcessing and a postProcessing.m routine to be run in order to produce plots of desired output.
