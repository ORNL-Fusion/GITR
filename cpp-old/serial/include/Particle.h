#include <random>

class Particle
{
   public:
      double x;   
      double y; 
      double z;
      double vx;   
      double vy; 
      double vz;
      double Z;
      double amu;
      double perpDistanceToSurface;
	std::mt19937 stream;	
      
      void BorisMove(double dt, double xMinV,double xMaxV,double yMin,double yMax,double zMin,double zMax);
      void Ionization(double dt);         
};
