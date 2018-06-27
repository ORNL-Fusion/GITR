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

## 2D Geometry Setup and File Structure
The GITR 2D geometries are made up of a collection of lines in the x-z plane which are then 
applied in 3D space with one of the following:
1. Bounded in the y direction (2 planes in the x-z plane)
2. Periodic in the y direction
3. Periodic in the theta direction (where the lines are in the r-z plane)

An example geometry cfg file used for the 2D geometry is shown below:
```
geom = 
{
   x1 = [1.0,2.0,2.0,1.0] 
   z1 = [0.0,0.0,1.0,1.0] 
   x2 = [2.0,2.0,1.0,1.0] 
   z2 = [0.0,1.0,1.0,0.0]
   slope = [0.000000,1000000000000.000000,0.000000,-1000000000000.000000]
   intercept = [0.00000,1000000000000.000000,1.0,1000000000000.000000]
   length = [1.0,1.0,1.0,1.0]   
   Z = [74.0,74.0,74.0,74.0,0.0]       
   y1 = -1.0;
   y2 = 1.0;
   periodic = 0;
}
```
