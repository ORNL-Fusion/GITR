/*

necessary components:

1. particles need:

everything in impurityParticleSource.initialConditions gitrInput.cfg from protoMPex example.
but nothing other than that module.


2. boris

3. placeholder geometry


flags to set/unset:

one particle is all that should be needed for this test



boris captures drifts - which is caused by rotating particles in fields that make them rotate?
and
normal Lorentz motion - from force fields acting on it

but not

diffusion - particle movement from collisions without explicitly modeling the collision event

final answer of this test: will be drift velocity


this test essentially sets up a uniform E and B field, and tests that the boris operator
correctly simulates the drift motion of the particles. It can be a one particle test.
Simply measure the start and end position, divide by time, get drift velocity, plug into the
equations and see if it equals the analytical solution

geometry related, checking algorithms. Respecting geometry stuff.


Priorities:

coupled simulations with SOLPS-ITER,
modularize configuration and geometry components ---> how?
make it compile for mac and linux and pass tests
finish unit tests
come up with new ones
finish changing all the build time options into runtime options
What options can I take out? Which ones are defunct?

Alyssa's stuff:

1. get your version of the repo sorted out and incorporated into the new one
2. make sure GITR_legacy_python works without issue with GITR.

2015 had an IPS fusion milestone etc. Is IPS that hard to use? New solution needed.
Let's find some other alternatives for coupled simulation alternatives.

How to automate it? Talk to Wael Elwasif.

boris efield and geometry and cylindrical options.
efield test is done in matlab
















*/
