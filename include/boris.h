#ifndef _BORIS_
#define _BORIS_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#include "thrust/extrema.h"
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include <algorithm>
#include "Particle.h"
#include "Boundary.h"
#include "interp2d.hpp"
#include <boost/timer/timer.hpp>

using namespace boost::timer;

template <typename T>
CUDA_CALLABLE_MEMBER
int sgn(T val) {
            return (T(0) < val) - (val < T(0));
}

CUDA_CALLABLE_MEMBER

#ifdef __CUDACC__
double getE ( double x0, double y, double z, double E[], Boundary *boundaryVector, int nLines ) {
#else
double getE ( double x0, double y, double z, double E[], std::vector<Boundary> &boundaryVector, int nLines ) {
#endif
	double Emag;
	double fd = 0.0;
	double pot = 0.0;
    int minIndex = 0;
    double minDistance = 1e12;
    int direction_type;
    double tol = 1e12;
    double point1_dist;
    double point2_dist;
    double perp_dist;
    double directionUnitVector[3] = {0.0,0.0,0.0};
    double vectorMagnitude;
    double max = 0.0;
    double min = 0.0;
    double angle = 0.0;
    double Er = 0.0;
    double Et = 0.0;
#if USECYLSYMM > 0
    double x = sqrt(x0*x0 + y*y);
#else
    double x = x0;
#endif    
    for (int j=0; j< nLines; j++)
    {
        if (boundaryVector[j].Z != 0.0)
        {
            point1_dist = sqrt((x - boundaryVector[j].x1)*(x - boundaryVector[j].x1) + 
                    (z - boundaryVector[j].z1)*(z - boundaryVector[j].z1));
            point2_dist = sqrt((x - boundaryVector[j].x2)*(x - boundaryVector[j].x2) + 
                                        (z - boundaryVector[j].z2)*(z - boundaryVector[j].z2));
            perp_dist = (boundaryVector[j].slope_dzdx*x - z + boundaryVector[j].intercept_z)/
                sqrt(boundaryVector[j].slope_dzdx*boundaryVector[j].slope_dzdx + 1);   

            if (point1_dist > point2_dist)
            {
                max = point1_dist;
                min = point2_dist;
            }
            else
            {
                max = point2_dist;
                min = point1_dist;
            }
            if (boundaryVector[j].length*boundaryVector[j].length + perp_dist*perp_dist >=
                    max*max)
            {
                boundaryVector[j].distanceToParticle =fabs( perp_dist);
                boundaryVector[j].pointLine = 1;
            }
            else
            {
                boundaryVector[j].distanceToParticle = min;
                if (boundaryVector[j].distanceToParticle == point1_dist)
                {
                    boundaryVector[j].pointLine = 2;
                }
                else
                {
                    boundaryVector[j].pointLine = 3;
                }
            }

            if (boundaryVector[j].distanceToParticle < minDistance)
            {
                minDistance = boundaryVector[j].distanceToParticle;
                minIndex = j;
                direction_type = boundaryVector[j].pointLine;
            }
        }
        else
        {
            boundaryVector[j].distanceToParticle = tol;
        }
    if (direction_type == 1)
    {
        if (boundaryVector[minIndex].slope_dzdx == 0)
        {
            directionUnitVector[0] = 0.0;
            directionUnitVector[1] = 0.0;
            directionUnitVector[2] = 1.0*sgn(boundaryVector[minIndex].z1 - z);
        }
        else if (fabs(boundaryVector[minIndex].slope_dzdx)>= 0.75*tol)
        {
            
            directionUnitVector[0] = boundaryVector[minIndex].x1 - x;
            directionUnitVector[1] = 0.0;
            directionUnitVector[2] = 0.0;
        }
        else
        {
            directionUnitVector[0] = -1.0*sgn(boundaryVector[minIndex].slope_dzdx)*sgn(perp_dist);
            directionUnitVector[1] = 0.0;
            directionUnitVector[2] = -1.0*directionUnitVector[0]/(boundaryVector[minIndex].slope_dzdx);
        }
    }
    else if (direction_type == 2)
    {
        directionUnitVector[0] = (boundaryVector[minIndex].x1 - x);
        directionUnitVector[1] = 0.0;
        directionUnitVector[2] = (boundaryVector[minIndex].z1 - z);
    }
    else
    {
        directionUnitVector[0] = (boundaryVector[minIndex].x2 - x);
        directionUnitVector[1] = 0.0;
        directionUnitVector[2] = (boundaryVector[minIndex].z2 - z);
    }

    vectorMagnitude = sqrt(directionUnitVector[0]*directionUnitVector[0] + directionUnitVector[1]*directionUnitVector[1]
                                + directionUnitVector[2]*directionUnitVector[2]);
    directionUnitVector[0] = directionUnitVector[0]/vectorMagnitude;
    directionUnitVector[1] = directionUnitVector[1]/vectorMagnitude;
    directionUnitVector[2] = directionUnitVector[2]/vectorMagnitude;
    
    angle = boundaryVector[minIndex].angle;    
    fd  =  0.996480862464192 +8.78424468259523e-04  * angle     -
           4.674013060191755e-4  * pow(angle,2) +
           2.727826261148182e-5  * pow(angle,3) - 
           7.141202673279612e-7  * pow(angle,4) +
           8.56348440384227e-9   * pow(angle,5) -
           3.91580557074662e-11  * pow(angle,6);
    pot = 3.0*boundaryVector[minIndex].ti;
        Emag = pot*(fd/(2.0*boundaryVector[minIndex].debyeLength)*exp(-minDistance/(2.0*boundaryVector[minIndex].debyeLength))+ (1.0 - fd)/(boundaryVector[minIndex].larmorRadius)*exp(-minDistance/boundaryVector[minIndex].larmorRadius) );
        Er = Emag*directionUnitVector[0];
        Et = Emag*directionUnitVector[1];
        E[2] = Emag*directionUnitVector[2];

    }
    //std::cout << "pos " << x << " " << y << " "<< z << " min Dist" << minDistance << "Efield " << Emag << std::endl;
#if USECYLSYMM > 0
            //if cylindrical geometry
            double theta = atan(y/x0);
            if(x < 0.0)
            {
                if(y > 0.0)
                {
                    theta = theta + 3.14159265358979;
                }
                else
                {
                    theta = sqrt(theta*theta) + 3.14159265358979;
                }
            }
  
            E[0] = cos(theta)*Er - sin(theta)*Et;
            E[1] = sin(theta)*Er + cos(theta)*Et;
#else
            E[0] = Er;
            E[1] = Et;
#endif
    return minDistance;
}

struct move_boris { 
#ifdef __CUDACC__
    Boundary *boundaryVector;
#else
    std::vector<Boundary> &boundaryVector;
#endif
int nR_Bfield;
int nZ_Bfield;
double * BfieldGridRDevicePointer;
double * BfieldGridZDevicePointer;
double * BfieldRDevicePointer;
double * BfieldZDevicePointer;
double * BfieldTDevicePointer;

    const double span;
    const int nLines;
#ifdef __CUDACC__
    move_boris(double _span, Boundary *_boundaryVector,int _nLines,
            int _nR_Bfield, int _nZ_Bfield,
            double * _BfieldGridRDevicePointer,
            double * _BfieldGridZDevicePointer,
            double * _BfieldRDevicePointer,
            double * _BfieldZDevicePointer,
            double * _BfieldTDevicePointer
              ) : span(_span), boundaryVector(_boundaryVector), nLines(_nLines), nR_Bfield(_nR_Bfield), nZ_Bfield(_nZ_Bfield), BfieldGridRDevicePointer(_BfieldGridRDevicePointer), BfieldGridZDevicePointer(_BfieldGridZDevicePointer),
    BfieldRDevicePointer(_BfieldRDevicePointer), BfieldZDevicePointer(_BfieldZDevicePointer), BfieldTDevicePointer(_BfieldTDevicePointer) {}
#else
    move_boris(double _span, std::vector<Boundary> &_boundaryVector, int _nLines,
            int _nR_Bfield, int _nZ_Bfield,
            double * _BfieldGridRDevicePointer,
            double * _BfieldGridZDevicePointer,
            double * _BfieldRDevicePointer,
            double * _BfieldZDevicePointer,
            double * _BfieldTDevicePointer
            
            ) : span(_span), boundaryVector(_boundaryVector), nLines(_nLines), nR_Bfield(_nR_Bfield), nZ_Bfield(_nZ_Bfield), BfieldGridRDevicePointer(_BfieldGridRDevicePointer), BfieldGridZDevicePointer(_BfieldGridZDevicePointer),
    BfieldRDevicePointer(_BfieldRDevicePointer), BfieldZDevicePointer(_BfieldZDevicePointer), BfieldTDevicePointer(_BfieldTDevicePointer) {}
#endif    

CUDA_CALLABLE_MEMBER    
void operator()(Particle &p) const { 
double initTime = 0.0;
double interpETime = 0.0;
double interpBTime = 0.0;
double operationsTime = 0.0;
cpu_timer timer;
cpu_times initTime0 = timer.elapsed();
	    if(p.hitWall == 0.0)
        {
            double v_minus[3]= {0.0, 0.0, 0.0};
	        double v[3]= {0.0, 0.0, 0.0};
	        double E[3] = {0.0, 0.0, 0.0};
	        double B[3] = {0.0,0.0,0.0};
            double br;
            double bz;
            double bt;
	        double dt = span;
	        double Bmag = 0.0;
	        double q_prime = 0.0;
            double coeff = 0.0;
            int nSteps = floor( span / dt + 0.5);
            double minDist = 0.0;
#if ODEINT ==	0  
            for ( int s=0; s<nSteps; s++ ) 
            {
	//          minDist = getE(p.xprevious,p.yprevious,p.zprevious,E,boundaryVector,nLines);
interp2dVector(&B[0],p.xprevious,p.yprevious,p.zprevious,nR_Bfield,nZ_Bfield,
               BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer);        
B[0] = 0.0;
B[2] = 0.0;

        /*    br = interp2dCombined(p.xprevious,p.yprevious,p.zprevious,nR_Bfield,nZ_Bfield,
                BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer);
            bz = interp2dCombined(p.xprevious,p.yprevious,p.zprevious,nR_Bfield,nZ_Bfield,
                BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldZDevicePointer);
            bt = interp2dCombined(p.xprevious,p.yprevious,p.zprevious,nR_Bfield,nZ_Bfield,
                BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldTDevicePointer);
	       
#if USECYLSYMM > 0
            //if cylindrical geometry
            double theta = atan(p.yprevious/p.xprevious);
            if(p.xprevious < 0.0)
            {
                if(p.yprevious > 0.0)
                {
                    theta = theta + 3.1415;
                }
                else
                {
                    theta = sqrt(theta*theta) + 3.1415;
                }
            }
  
            B[0] = cos(theta)*br - sin(theta)*bt;
            B[2] = bz;
            B[1] = sin(theta)*br + cos(theta)*bt;
#else
            B[0] = br;
            B[2] = bz;
            B[1] = bt;
#endif*/
           // std::cout << "Particle Position " << p.xprevious << " " << p.yprevious << std::endl;
         //std::cout << "Bfield out " << B[0] << " " << B[1] <<" " << B[2] << std::endl;
            Bmag = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
	        q_prime = p.charge*1.60217662e-19/(p.amu*1.6737236e-27)*dt*0.5;
            coeff = 2*q_prime/(1+(q_prime*Bmag)*(q_prime*Bmag));
            
                v[0] = p.vx;
                v[1] = p.vy;
	            v[2] = p.vz;
                  
	            	  
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
	                                                                    
	            p.x = p.xprevious + v[0]*dt;
	            p.y = p.yprevious + v[1]*dt;
	            p.z = p.zprevious + v[2]*dt;
	                                                                    
	                             
	            						
	            p.vx = v[0];
	            p.vy = v[1];
	            p.vz = v[2];
    	    }
#endif

#if ODEINT == 1
        //RK4 integrator
        double m = p.amu*1.6737236e-27;
        double q_m = p.charge*1.60217662e-19/m;
        double r[3]= {0.0, 0.0, 0.0};
        double r2[3]= {0.0, 0.0, 0.0};
        double r3[3]= {0.0, 0.0, 0.0};
        double r4[3]= {0.0, 0.0, 0.0};
        double v2[3]= {0.0, 0.0, 0.0};
        double v3[3]= {0.0, 0.0, 0.0};
        double v4[3]= {0.0, 0.0, 0.0};
        double k1r[3]= {0.0, 0.0, 0.0};
        double k2r[3]= {0.0, 0.0, 0.0};
        double k3r[3]= {0.0, 0.0, 0.0};
        double k4r[3]= {0.0, 0.0, 0.0};
        double k1v[3]= {0.0, 0.0, 0.0};
        double k2v[3]= {0.0, 0.0, 0.0};
        double k3v[3]= {0.0, 0.0, 0.0};
        double k4v[3]= {0.0, 0.0, 0.0};
                v[0] = p.vx;
                v[1] = p.vy;
	            v[2] = p.vz;

                r[0] = p.xprevious;
                r[1] = p.yprevious;
	            r[2] = p.zprevious;
cpu_times initTime1 = timer.elapsed();
initTime = initTime + (initTime1.wall - initTime0.wall);
            for ( int s=0; s<nSteps; s++ ) 
            {
                //std::cout << "inside rk4 " << std::endl;
cpu_times operationsTime0 = timer.elapsed();
cpu_times interpETime0 = timer.elapsed();
	          minDist = getE(r[0],r[1],r[2],E,boundaryVector,nLines);
cpu_times interpBTime0 = timer.elapsed();
                interp2dVector(&B[0],r[0],r[1],r[2],nR_Bfield,nZ_Bfield,
               BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer);        
cpu_times interpBTime1 = timer.elapsed();
interpETime = interpETime + (interpBTime0.wall - interpETime0.wall);
interpBTime = interpBTime + (interpBTime1.wall - interpBTime0.wall);
                //B[0] = 0.0;
                //B[2] = 0.0;
                k1r[0] = v[0]*dt;
                k1r[1] = v[1]*dt;
                k1r[2] = v[2]*dt;

                k1v[0] = dt*q_m*(E[0] + (v[1]*B[2] - v[2]*B[1]));
                k1v[1] = dt*q_m*(E[1] + (v[2]*B[0] - v[0]*B[2]));
                k1v[2] = dt*q_m*(E[2] + (v[0]*B[1] - v[1]*B[0]));
            
                /*std::cout << " v0 " << v[0] << std::endl;
                std::cout << " dt " << dt << std::endl;
                std::cout << " q_m " << q_m << std::endl;
                std::cout << " B " << B[0] << " " << B[1] << " " << B[2] << std::endl;
                std::cout << " k1r " << k1r[0] << std::endl;
                std::cout << " k1v " << k1v[0] << std::endl;*/
                r2[0] = r[0] + k1r[0]*0.5;
                r2[1] = r[1] + k1r[0]*0.5;
                r2[2] = r[2] + k1r[0]*0.5;
               
                v2[0] = v[0] + k1v[0]*0.5;
                v2[1] = v[1] + k1v[0]*0.5;
                v2[2] = v[2] + k1v[0]*0.5;
interpETime0 = timer.elapsed();
	          minDist = getE(r2[0],r2[1],r2[2],E,boundaryVector,nLines);
interpBTime0 = timer.elapsed();
                interp2dVector(&B[0],r2[0],r2[1],r2[2],nR_Bfield,nZ_Bfield,
                        BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer);        
interpBTime1 = timer.elapsed();
interpETime = interpETime + (interpBTime0.wall - interpETime0.wall);
interpBTime = interpBTime + (interpBTime1.wall - interpBTime0.wall);
                k2r[0] = v2[0]*dt;
                k2r[1] = v2[1]*dt;
                k2r[2] = v2[2]*dt;
                k2v[0] = dt*q_m*(E[0] + (v2[1]*B[2] - v2[2]*B[1]));
                k2v[1] = dt*q_m*(E[1] + (v2[2]*B[0] - v2[0]*B[2]));
                k2v[2] = dt*q_m*(E[2] + (v2[0]*B[1] - v2[1]*B[0]));
                
                r3[0] = r[0] + k2r[0]*0.5;
                r3[1] = r[1] + k2r[0]*0.5;
                r3[2] = r[2] + k2r[0]*0.5;
               
                v3[0] = v[0] + k2v[0]*0.5;
                v3[1] = v[1] + k2v[0]*0.5;
                v3[2] = v[2] + k2v[0]*0.5;
interpETime0 = timer.elapsed();
	          minDist = getE(r3[0],r3[1],r3[2],E,boundaryVector,nLines);
interpBTime0 = timer.elapsed();
                interp2dVector(&B[0],r3[0],r3[1],r3[2],nR_Bfield,nZ_Bfield,
                        BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer);        
                
interpBTime1 = timer.elapsed();
interpETime = interpETime + (interpBTime0.wall - interpETime0.wall);
interpBTime = interpBTime + (interpBTime1.wall - interpBTime0.wall);
                k3r[0] = v3[0]*dt;
                k3r[1] = v3[1]*dt;
                k3r[2] = v3[2]*dt;
                k3v[0] = dt*q_m*(E[0] + (v3[1]*B[2] - v3[2]*B[1]));
                k3v[1] = dt*q_m*(E[1] + (v3[2]*B[0] - v3[0]*B[2]));
                k3v[2] = dt*q_m*(E[2] + (v3[0]*B[1] - v3[1]*B[0]));
                
                r4[0] = r[0] + k3r[0];
                r4[1] = r[1] + k3r[0];
                r4[2] = r[2] + k3r[0];
               
                v4[0] = v[0] + k3v[0];
                v4[1] = v[1] + k3v[0];
                v4[2] = v[2] + k3v[0];
interpETime0 = timer.elapsed();
            
	          minDist = getE(r4[0],r4[1],r4[2],E,boundaryVector,nLines);
interpBTime0 = timer.elapsed();
                interp2dVector(&B[0],r4[0],r4[1],r4[2],nR_Bfield,nZ_Bfield,
                        BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer);        
interpBTime1 = timer.elapsed();
interpETime = interpETime + (interpBTime0.wall - interpETime0.wall);
interpBTime = interpBTime + (interpBTime1.wall - interpBTime0.wall);
                k4r[0] = v4[0]*dt;
                k4r[1] = v4[1]*dt;
                k4r[2] = v4[2]*dt;
                k4v[0] = dt*q_m*(E[0] + (v4[1]*B[2] - v4[2]*B[1]));
                k4v[1] = dt*q_m*(E[1] + (v4[2]*B[0] - v4[0]*B[2]));
                k4v[2] = dt*q_m*(E[2] + (v4[0]*B[1] - v4[1]*B[0]));
        
                p.x = r[0] + (k1r[0] + 2*k2r[0] + 2*k3r[0] + k4r[0])/6;
                p.y = r[1] + (k1r[1] + 2*k2r[1] + 2*k3r[1] + k4r[1])/6;
                p.z = r[2] + (k1r[2] + 2*k2r[2] + 2*k3r[2] + k4r[2])/6;
                p.vx = v[0] + (k1v[0] + 2*k2v[0] + 2*k3v[0] + k4v[0])/6;
                p.vy = v[1] + (k1v[1] + 2*k2v[1] + 2*k3v[1] + k4v[1])/6;
                p.vz = v[2] + (k1v[2] + 2*k2v[2] + 2*k3v[2] + k4v[2])/6;
                //std::cout << " r0 " << r[0] << std::endl;
                //std::cout << " px " << p.x << std::endl;
cpu_times operationsTime1 = timer.elapsed();
operationsTime = operationsTime1.wall - operationsTime0.wall;

std::cout << "Operations Time: " << operationsTime <<std::endl;
std::cout << "Efield Interpolation Time: " << interpETime <<std::endl;
std::cout << "Bfield Interpolation Time: " << interpBTime <<std::endl;
std::cout << "Init Time: " << initTime <<std::endl;
            }
#endif
        }
    } 
};

#endif
