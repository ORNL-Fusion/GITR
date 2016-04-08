#define EInterpNumber 2

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <mpi.h>
#include "h1.h"
using namespace std;

void Efield(double E[], double perpDistanceToSurface)
{
#if EInterpNumber == 0

#elif EInterpNumber == 2
	double Emag;
	double Lc = 1.0e5;
	double s = Lc*0.5 - perpDistanceToSurface;//Distance from stagnation point in direction of surface along Bfield
	
	
	double surfaceDirection_unit[3] =  {0.866, 0, -0.50};//{0.8192, 0,  -0.5736};//;//{0.9063, 0, -0.4226};
	double Epara = 0.0;
	
	Epara = 20.0*(Lc/(2*s*sqrt(Lc*Lc*0.25 - s*s)) - 1.0/s);
	if(perpDistanceToSurface == 0.0){
	Epara = 0.0;
	}
	//std::cout << "s: " << s<< std::endl;
	//std::cout << "Epara: " << Epara << std::endl;
	Emag = 60*(0.8357/(2.1026e-05)*exp(-perpDistanceToSurface/2.1026e-05)+ (1.0 - 0.8357)/(3.64479e-4)*exp(-perpDistanceToSurface/3.64479e-4) );
	E[0] = Emag*surfaceDirection_unit[0];
	E[1] = Emag*surfaceDirection_unit[1];
	E[2] = Emag*surfaceDirection_unit[2] - Epara;
	//std::cout << "Ez: " << E[2] << std::endl;

#endif

}