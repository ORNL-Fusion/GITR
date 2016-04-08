#include "h1.h"
#include <mpi.h>
#include <math.h>
#include <iostream>
#include <random>
#include <iomanip>
#include <cstdlib>
#include <libconfig.h++>
#include <cmath>

//MESH

void MESH(int myID, int nWRs,double r[], double z[], double dr,double dz, int localMz, int Mrp1, double Rin, double Rout, double Ztop, double Ar[], double Az[], double Vij[])
{
double start = double(localMz)*double(myID-1);

 
  for(int i=1 ; i<(localMz + 3) ; i++)
	{
	z[i-1] = (2*(start+double(i))-1)*dz/2.0;
	}
	



  for(int i=1 ; i<(Mrp1) ; i++)
	{
	
	r[i-1] = (2*double(i)-1)*dr/2 +Rin;
		Ar[i-1] = 2.0*M_PI*(r[i-1]-dr/2.0)*dz;
			//std::cout << Ar[i] << std::endl;
			Az[i-1] = 2.0*M_PI*r[i-1]*dr;
			Vij[i-1] = 2.0*M_PI*r[i-1]*dr*dz;

	}
	

  Ar[Mrp1-1] = 2*M_PI*(r[Mrp1-2]+dr/2)*dz;
  

  
}




//INIT
void INIT(int myID, int nWRs, int nP, Particle p[])
{

	int ParticlesPerProcessor;
	
	ParticlesPerProcessor = nP/nWRs;
	//std::cout<< "np " << nP << "nWRs " << nWRs << "ppp " << ParticlesPerProcessor << std::endl;
	for(int i=0 ; i<ParticlesPerProcessor ; i++)
	{
	//std::cout<< "Processor: " << myID << " particle, " << i<< std::endl;
	p[i].x = 0.0;
	p[i].y = 0.0;
	p[i].z = 0.0;
	p[i].vx = -3225.7;
	p[i].vy = 0.0;
	p[i].vz = 0.0;
	p[i].Z = 0.0;
	p[i].amu = 184.0;
	p[i].perpDistanceToSurface = 0.0;
	//std::cout<< "amu: " << p[i].amu << std::endl;
		std::random_device rd;
	std::mt19937 mt(rd());
	p[i].stream = mt;
	}
	
}

