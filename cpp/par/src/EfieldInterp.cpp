#define EInterpNumber 2

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <mpi.h>
#include "h1.h"
using namespace std;

void Efield(double E[], double perpDistanceToSurface,double z, double x)
{
#if EInterpNumber == 0

#elif EInterpNumber == 2
	double Emag;
	double Lc = 1.0e5;
	double dl = perpDistanceToSurface;//fabs(fabs(z)-1.73205*fabs(x)); //distance from surface alon Bfield line
	double s = Lc*0.5 - dl;//Distance from stagnation point in direction of surface along Bfield
	//std::cout << "dl and s : " << dl << "  " << s << std::endl;
	
	double surfaceDirection_unit[3] =  {0.866, 0, -0.50};//{0.8192, 0,  -0.5736};//;//{0.9063, 0, -0.4226};
	double Epara = 0.0;
	double Erad = 0.0;
	double f_s;
	double Rfac = 1.0;
	
	//Epara = 20.0*(Lc/(2*s*sqrt(Lc*Lc*0.25 - s*s)) - 1.0/s);
	//Epara = 0.5*20.0*(Lc/((Lc - s)*sqrt(Lc*Lc - (Lc-s)*(Lc-s))) - 1.0/(Lc-s));
	       Epara   =   0.5 * 20.0 
                  * ( Lc/((Lc - dl)*sqrt(Lc*Lc 
                           - (Lc - dl)*(Lc - dl))) -1/(Lc - dl) );
	if(perpDistanceToSurface == 0.0){
	Epara = 0.0;
	}
	
// 	         if(4.0*Lc*Lc) -4.0*(Lc - fabs(s))*(Lc - fabs(s)) > 0.0)
// 	 {
//            f_s = log(2*Lc 
//                   + sqrt(pow(2*Lc,2.) -4.*pow((Lc - fabs(s)),2.)))
//                  - log(2*2*Lc );
//          }
//          else
//            f_s =log(0.5);
// 	Erad =  e_param->Rfac * f_s*(Te_local/nT_ptr->Te_abf);
	
	
	//std::cout << "s: " << s<< std::endl;
	//std::cout << "Epara: " << Epara << std::endl;
	Emag = 60*(0.8357/(2.0*1.05058e-05)*exp(-perpDistanceToSurface/(2.0*1.05058e-05))+ (1.0 - 0.8357)/(3.64479e-4)*exp(-perpDistanceToSurface/(3.64479e-4)) );
	E[0] = Emag*surfaceDirection_unit[0];
	E[1] = Emag*surfaceDirection_unit[1];
	E[2] = Emag*surfaceDirection_unit[2];
	//std::cout << "Ez: " << E[2] << std::endl;

#endif

}