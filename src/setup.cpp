#include "h1.cuh"
#include <iostream>

#include <iomanip>
#include <cstdlib>
#include <libconfig.h++>
#include <cmath>


//INIT
void INIT(int nP, Particle p[], libconfig::Config &cfg)
{

double x_start = cfg.lookup("impurityParticleSource.initialConditions.x_start");
double y_start = cfg.lookup("impurityParticleSource.initialConditions.y_start");
double z_start = cfg.lookup("impurityParticleSource.initialConditions.z_start");

double energy_eV_x_start = cfg.lookup("impurityParticleSource.initialConditions.energy_eV_x_start");
double energy_eV_y_start = cfg.lookup("impurityParticleSource.initialConditions.energy_eV_y_start");
double energy_eV_z_start = cfg.lookup("impurityParticleSource.initialConditions.energy_eV_z_start");

double impurity_amu = cfg.lookup("impurityParticleSource.initialConditions.impurity_amu");
double impurity_Z = cfg.lookup("impurityParticleSource.initialConditions.impurity_Z");

double vx = energy_eV_x_start/std::abs(energy_eV_x_start)*sqrt(2.0*std::abs(energy_eV_x_start)*1.60217662e-19/(impurity_amu*1.6737236e-27));
double vy = energy_eV_y_start/std::abs(energy_eV_y_start)*sqrt(2.0*std::abs(energy_eV_y_start)*1.60217662e-19/(impurity_amu*1.6737236e-27));
double vz = energy_eV_z_start/std::abs(energy_eV_z_start)*sqrt(2.0*std::abs(energy_eV_z_start)*1.60217662e-19/(impurity_amu*1.6737236e-27));

if(energy_eV_x_start == 0.0) vx = 0.0;
if(energy_eV_y_start == 0.0) vy = 0.0;
if(energy_eV_z_start == 0.0) vz = 0.0;

	//std::cout<< "np " << nP << "nWRs " << nWRs << "ppp " << ParticlesPerProcessor << std::endl;
	for(int i=0 ; i<nP; i++)
	{
	//std::cout<< "Processor: " << myID << " particle, " << i<< std::endl;
	p[i].x = x_start;
	p[i].y = y_start;
	p[i].z = z_start;
	p[i].vx = vx;
	p[i].vy = vy;
	p[i].vz = vz;
	p[i].Z = impurity_Z;
	p[i].amu = impurity_amu;
	p[i].perpDistanceToSurface = 0.0;
	//std::cout<< "amu: " << p[i].amu << std::endl;
	}
	
}
