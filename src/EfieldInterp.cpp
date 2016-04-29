#define EInterpNumber 2

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "h1.cuh"
using namespace std;

void Efield(double E[], double perpDistanceToSurface)
{
#if EInterpNumber == 0

#elif EInterpNumber == 2
	double Emag;
	double surfaceDirection_unit[3] =  {0.866, 0, -0.50};//{0.8192, 0,  -0.5736};//;//{0.9063, 0, -0.4226};
	
	Emag = 60*(0.8357/(2.1026e-05)*exp(-perpDistanceToSurface/2.1026e-05)+ (1.0 - 0.8357)/(3.64479-4)*exp(-perpDistanceToSurface/3.64479-4) );
	E[0] = Emag*surfaceDirection_unit[0];
	E[1] = Emag*surfaceDirection_unit[1];
	E[2] = Emag*surfaceDirection_unit[2];

#endif

}
