# GITR-CPC-Example
A published example case for the global impurity transport code.
The paper in which this example is documented can be found here: https://www.sciencedirect.com/science/article/abs/pii/S0010465521000369

<img src="images/DepTracks2.png " width="500">

## How to run
```
git clone https://github.com/ORNL-Fusion/GITR.git
cd GITR
mkdir build
cd build
cmake ..
make -j
cd ../../
git clone https://github.com/ORNL-Fusion/GITR-CPC-Example.git
cd GITR-CPC-Example
mkdir output
../GITR/build/GITR
```

## Example Details
An example case to demonstrate GITR and the code outputs is now presented. This example has some features that are relevant to fusion erosion scenarios but is also simplified for demonstration purposes and satisfies the description given previously in Section 5: Surface Model and Particle Balance. In this case a D plasma (mass=2 AMU charge=1) withTi =Te =20eV andni =ne =1e19m−3 fills the geometry made up of a plane angled at 30 degrees with respect to a y-z plane and passing through the origin. Connecting walls
at x=+-1.5cm, y=+-1.5cm and z=3cm complete the closed geometry. The mesh resolution is set to 3 mm, and the surface geometry is comprised of 1280 boundary and surface (triangle) elements in which all boundaries are W material surfaces. This 3D surface geometry is depicted in Figure 10. The magnetic field is fixed to 2 Tesla in the x-direction. Ionization by electron impact is switched on using coefficients from the ADAS database [9]. The sheath electric field is applied according to the model specified in reference [16]. A simple surface model for reflection (R=0.5) and sputtering yield (Y=0.1) are used where any sputtered/reflected particle leaves the surface with 8 eV of energy, an azimuthal angle of φ = 45 deg and θ = ξ∗360deg where ξ is a random number between 0 and 1. The initial source of eroded particles is a point source at the origin with energy in the z-direction of 10 eV with a magnitude set to 1 particle per second and with nP=1e5 particles. Therefore, each computational particle represents 1e-5 physical particles per second. A time-step of 1e-10 s is used for 1e7 time steps for a simulated time of 1 ms. The problem has particles that are launched vertically (in the z-direction), which are ionized and then follow gyromotion orbits before accelerating back to the surface where the particles sputter and reflect to be launched from the surface again. The particles are weighted according to each process, and the weights are accounted for in binning for densities of impurities as well as impacts on the material surface.
Typical results of such a simulation can
yield particle tracks, charge-resolved volumetric density histograms, surface erosion,
deposition, sputtering and reflection maps,
and energy-angle distributions at the surface. Returning to the analytical example of Surface Model and Particle balance
where grossDeposition = nP ∗ 1−R = 1−(Y +R)
1.25e5 the GITR result predicts a value of 1.24997e5. This value is in very close in agreement, in which the difference can be attributed to the accumulation of machine adding accumulated error.
 
 
<img src="images/dep.png " width="500"> 
Figure 10: GITR example case simulated geometry with side cut away. Gross deposited mass rate (colored surfaces) and example particle tracks (red lines) are shown.
 
When binning or summing to infinity as is done in Equation 9, if the difference in bin size (grossDeposition) and the current particle weight differs by more than the precision being used, the weight will not be effectively added to the bin, thus truncating the sum before infinity. In most cases, the accuracy of these operations can be alleviated by increasing the fidelity of the bin size (either surface or volume). This reduces the magnitude of the total size of the bin and reduces the error in adding small weights to the bins.
Additional sample output from the GITR simulation is shown in Figures 11 and 12.
Since all surface elements in the simulation are W material surfaces that are in contact with plasma, there exists a sheath electric field structure according to [16], which serves to accelerate simulated impurity ions toward the local surface. The total potential of the surface is set to be φ = 3kTe/e where k is the Boltzmann constant, Te is electron temperature and e is the electron charge. Therefore ions entering the sheath can experience accelerations adding kinetic energy to the particle of up to Qe(3kTe) where Q is the impurity ion charge number. Since some particles
are ionized within the sheath region, in some cases, only a fraction of this potential energy is given to the impurity ion traveling back to the surface. The summed ion energy-angle distribution (IEAD) for all of the surfaces in the example problem is shown in Figure 11. This result demonstrates that impurity ion charges of approximately 1 to 6 are impacting the material surfaces with peaks corresponding to Qe(3kTe) (60, 120, 180, 240, 300, 360 eV). The surface IEAD demonstrates a tail of lower energy impacts which represent particles that are ionized near the surface, only receiving a fraction of the sheath potential acceleration before striking the surface.

<img src="images/iead.png " width="500">
Figure 11: GITR example case ion-surface impact energy-angle distribution summed over all surfaces.

Another example of GITR simulation output is the W impurity ion density profile shown in Figure 12. This profile demonstrates the impurity density of all charge states averaged in the y-direction. The particle source location of z=0 and x=0 can be seen with respect to the surface plane plotted as a white line. It shows that the launch trajectory of the impurity atoms in the z-direction dominates the density profile with a peak several mm from the surface. There are tails of density in the +/- x-direction demonstrating the confinement of particles to parallel transport along the magnetic field.

<img src="images/dens.png " width="500">
  Figure 12: GITR example case in which the W density profile (summed over all charge states) is integrated over y-domain.
