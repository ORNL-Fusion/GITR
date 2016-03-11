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
	double surfaceDirection[3] = {1.7321, 0, -1.00};
	double surfaceDirection_unit[3] = {0.866, 0, -0.50};
	
	Emag = 60*(0.85/(2.1e-5)*exp(-perpDistanceToSurface/2.1e-5)+ 0.15/(3.64479e-4)*exp(-perpDistanceToSurface/3.64479e-4) );
	E[0] = Emag*surfaceDirection_unit[0];
	E[1] = Emag*surfaceDirection_unit[1];
	E[2] = Emag*surfaceDirection_unit[2];

#endif

}