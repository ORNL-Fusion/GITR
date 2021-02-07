# Straight Line Particle Track Tests/Examples
The straight line (neutral particle streaming) tests are formulated in order
to test geometry operations during development to make sure things don't break.
This makes up a handful of tests to see when/where particles intersec the surface
for different geometry types. These geometry types include:
1. 2Dgeom - 2D cube geometry bounded in the 3rd dimension
2. 2D cube geometry periodic in the 3rd dimension
3. 2DgeomCyl - square geometry cylindrical
4. 3Dgeom - 3D cube geometry
5. 3D geometry intersection with a vertex
6. 2D DIII-D geometry cylindrical
7. 2D DIII-D geometry cylindrical using hash
8. 3D DIII-D geometry using hash

## Input File and Particle Initial Conditions
with the flags for ionization off, most of the operators will not be applied, keep them off
the particle source energy angle file =0 so it is a point source
So our most relevant section of the input file is the initial conditions section:
```
#input/gitrInput.cfg
.
.
.
impurityParticleSource = 
{
        nP = 1;
        sourceStrength = 1E+19;
        Z = 74.0;
    source_material_Z = 74;
    source_material_SurfaceBindingEnergy = 11.75;

        initialConditions = 
        {       
                x_start = 1.5;
                y_start = 0.0;
                z_start = 0.0000000001;
                energy_eV_x_start = 6.0;
                energy_eV_y_start = 0.0;
                energy_eV_z_start = 6.0;
                impurity_amu = 12.0;
                impurity_Z = 6.0;
            charge = 0.0;
                energy_eV=10.0;
                phi = 0.0;
                theta = 0.0;
    }


.
.
.
timeStep = 
{
    dt = 1E-5;
        nPtsPerGyroOrbit = 10000.0;
        ionization_nDtPerApply = 1;
        collision_nDtPerApply = 5;
        nT = 1;
}
```

## 2D Geometry Setup and File Structure
The GITR 2D geometries are made up of a collection of lines in the x-z plane which are then 
applied in 3D space with one of the following:
1. Bounded in the y direction (2 planes in the x-z plane)
2. Periodic in the y direction
3. Periodic in the theta direction (where the lines are in the r-z plane)

An example geometry cfg file used for the 2D geometry is shown below:
```
#2Dgeom/input/gitrGeometry.cfg
geom = 
{
   x1 = [1.0,2.0,2.0,1.0] 
   z1 = [0.0,0.0,1.0,1.0] 
   x2 = [2.0,2.0,1.0,1.0] 
   z2 = [0.0,1.0,1.0,0.0]
   slope = [0.0,1000000000000.0,0.0,-1000000000000.0]
   intercept = [0.0,1000000000000.0,1.0,1000000000000.0]
   length = [1.0,1.0,1.0,1.0]   
   Z = [74.0,74.0,74.0,74.0,0.0]       
   y1 = -1.0;
   y2 = 1.0;
   periodic = 0;
}

```
The end points of the lines are shown to make up 4 lines:
(1,0) to (2,0) - slope and intercept of 0 length of 1
(2,0) to (2,1) - large slope and intercept (vertical line approximation)
(2,1) to (1,1) - zero slope and intercept of 1
(1,1) to (1,0) - large negative slope and intercept (vertical line approximation)

The material Z has one more element than lines to describe the material of bounding plane in the y direction.
y1 and y2 can describe the bounds in either y or theta depending on whether the geometry is being used in the cylindrical application or not.

periodic=1 is on and periodic=0 is off.

## 3D Geometry Setup and File Structure
The GITR 3D geometries are made up of a collection of triangles in 3D space. For the 3D geometries, several coefficients are pre-computed and stored in the geometry file.

An example geometry cfg file used for the 3D geometry is shown below:
```
#3Dgeom/input/gitrGeometry.cfg
geom = 
{
   x1 = [1.0,1.0,1.0,1.0,2.0,2.0,1.0,1.0,1.0,1.0,1.0,1.0 ] 
   y1 = [-1.0e-1,-1.0e-1,1.0e-1,1.0e-1,-1.0e-1,-1.0e-1,-1.0e-1,-1.0e-1,-1.0e-1,-1.0e-1,-1.0e-1,-1.0e-1 ] 
   z1 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0 ] 
   x2 = [2.0,1.0,2.0,1.0,2.0,2.0,1.0,1.0,2.0,1.0,2.0,1.0 ] 
   y2 = [-1.0e-1,-1.0e-1,1.0e-1,1.0e-1,1.0e-1,-1.0e-1,1.0e-1,-1.0e-1,-1.0e-1,1.0e-1,-1.0e-1,1.0e-1 ] 
   z2 = [0.0,1.0,0.0,1.0,0.0,1.0,0.0,0.0,1.0,1.0,0.0,0.0 ] 
   x3 = [2.0,2.0,2.0,2.0,2.0,2.0,1.0,1.0,2.0,2.0,2.0,2.0 ] 
   y3 = [-1.0e-1,-1.0e-1,1.0e-1,1.0e-1,1.0e-1,1.0e-1,1.0e-1,1.0e-1,1.0e-1,1.0e-1,1.0e-1,1.0e-1 ] 
   z3 = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0 ] 
   a = [0.0,0.0,0.0,0.0,2.0e-1,-2.0e-1,2.0e-1,0.0,0.0,0.0,0.0,0.0 ] 
   b = [-1.0,1.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 ] 
   c = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,2.0e-1,-2.0e-1,2.0e-1,-2.0e-1 ] 
   d = [-1.0e-1,1.0e-1,1.0e-1,-1.0e-1,-4.0e-1,4.0e-1,-2.0e-1,-0.0,-2.0e-1,2.0e-1,-0.0,-0.0 ] 
   plane_norm = [1.0,1.0,1.0,1.0,2.0e-1,2.0e-1,2.0e-1,0.0,2.0e-1,2.0e-1,2.0e-1,2.0e-1 ] 
   BCxBA = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,1.0 ] 
   CAxCB = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,1.0 ] 
   area = [5.0e-1,5.0e-1,5.0e-1,5.0e-1,1.0e-1,1.0e-1,1.0e-1,0.0,1.0e-1,1.0e-1,1.0e-1,1.0e-1 ] 
   Z = [23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0 ] 
   theta0=0.0;
   theta1=0.0;
   periodic=0;
}
}
```
The three vertices of each triangle are made up by (x1,y1,z1), (x2,y2,z2), (x3,y3,z3) corresponding to the vectors A,B,C. This makes up the geometry shown below.
3dgeom.png

The coefficients a,b,c, and d correspond to the plane equation for each triangle ax+by+cz+d=0 and plane_norm = sqrt(a^2 + b^2 + c^2)

BCxBA and CAxCB are the sign of the cross product of these vectors. This helps give the code a relative direction of the triangle normal and parallel directions.

periodic=1 is on and periodic=0 is off.
