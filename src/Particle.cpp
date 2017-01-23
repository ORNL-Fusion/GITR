#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include "h1.cuh"
using namespace std;


void Particle::BorisMove(double dt,double xMinV,double xMaxV,double yMin,double yMax,double zMin,double zMax)
{

	double v_minus[3]= {0, 0, 0};
	double v[3]= {0, 0, 0};
	double B[3] = {0, 0, -2};
	double E[3] = {0, 0, 0};
	double r[3] = {0, 0, 0};
	double surfaceDirection[3] = {1.7321, 0, -1.00};//{1.4281, 0, -1.00};//{2.1445, 0, -1.00};//{1.7321, 0, -1.00};
	double surfaceDirection_unit[3] = {0.866, 0, -0.50};//{0.9063, 0, -0.4226}; // {0.866, 0, -0.50};
	double t;
	double surface_dz_dx;
	
	surface_dz_dx = surfaceDirection[0];
	

	perpDistanceToSurface = ( -surfaceDirection[0]*x + z )/2.0;//1.7434;// 2.0 2.3662;
	//std::cout << "Perp dist" << perpDistanceToSurface << std::endl;
	//std::cout << " Particle position in xyz " << x << " " << y << " " << z << std::endl;
	//perpd = ( -0.5773*x + z )/1.1547;
Efield( E, perpDistanceToSurface);
//cout << Z <<  E[0] << E[1] << E[2] << endl;
	
	double Bmag = 2;
	
	           double     q_prime = Z*1.60217662e-19/(amu*1.6737236e-27)*dt*0.5;
              double  coeff = 2*q_prime/(1+(q_prime*Bmag)*(q_prime*Bmag));
              
              v[0] = vx;
              v[1] = vy;
              v[2] = vz;
              
              r[0] = x;
              r[1] = y;
              r[2] = z;	
              
              				  v_minus[0] = v[0] + q_prime*E[0];
                              v_minus[1] = v[1] + q_prime*E[1];
                              v_minus[2] = v[2] + q_prime*E[2];
                
                v[0] = v_minus[0] + q_prime*(v_minus[1]*B[2] - v_minus[2]*B[1]);
                v[1] = v_minus[1] + q_prime*(v_minus[2]*B[0] - v_minus[0]*B[2]);
                v[2] = v_minus[2] + q_prime*(v_minus[0]*B[1] - v_minus[1]*B[0]);
                
                v[0] = v_minus[0] + coeff*(v[1]*B[2] - v[2]*B[1]);
                v[1] = v_minus[1] + coeff*(v[2]*B[0] - v[0]*B[2]);
                v[2] = v_minus[2] + coeff*(v[0]*B[1] - v[1]*B[0]);
                
                v[0] = v[0] + q_prime*E[0];
                v[1] = v[1] + q_prime*E[1];
                v[2] = v[2] + q_prime*E[2];
                

                r[0] = x + v[0]*dt;
                r[1] = y + v[1]*dt;
                r[2] = z + v[2]*dt;
                
                perpDistanceToSurface = ( -surfaceDirection[0]*r[0] + r[2] )/2.0;
                
//                 			if (perpDistanceToSurface < 0)
// 			{
// 			            t = (surface_dz_dx*x - z)/(-surface_dz_dx*(x -r[0]) + (z -r[2]));
//                         x = x + (x - r[0])*t;
//                         y = y + (y - r[1])*t;
//                         z = z + (z - r[2])*t;
// 
// 			}
// 			else
// 			{
			x = r[0];
			y = r[1];
			z = r[2];
// 			}
			
			vx = v[0];
			vy = v[1];
			vz = v[2];
//    std::cout << "Pos: " << x << y << z << std::endl;
	
}

void Particle::Ionization(double dt)
{
	double tion;
	double P1;
	double Coeffs[10] = {3.0875e-13, 9.6970e-14,3.8631e-14,2.0649e-14,3.5021e-15,1.6037e-15,7.0230e-17,1.7442e-17,6.1966e-18,1.8790e-18};
	double density = 1e19;
	double Temp_eV = 20;
	tion = 1/(Coeffs[int(floor(Z+ 0.5f))]*density);
	
	P1 = 1-exp(-dt/tion);

	double r1=((double)rand()/(double)RAND_MAX);
	
		//cout << " Z " << Z << " P1 " << P1 << " r1 " << r1 << endl;
	
// 	if( Z == 0)
// 	{ P1 = 1.0;
// 	}
// 	else P1 = 0.0;
	
	if(r1 <= P1)
	{
		Z = Z+1;
		//std::cout << "ionization" << std::endl;
		//std::cout << " Particle position in xyz " << x << " " << y << " " << z << std::endl;
	} 
	
	//std::cout << "P1 and r1" << P1 << r1 << std::endl;
	}
