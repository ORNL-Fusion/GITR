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
#if USE_BOOST
#include <boost/timer/timer.hpp>
using namespace boost::timer;
#endif
template <typename T>
CUDA_CALLABLE_MEMBER
int sgn(T val) {
            return (T(0) < val) - (val < T(0));
}

CUDA_CALLABLE_MEMBER
void vectorAdd(float A[], float B[],float C[])
{
    C[0] = A[0] + B[0];
    C[1] = A[1] + B[1];
    C[2] = A[2] + B[2];
}

CUDA_CALLABLE_MEMBER
void vectorScalarMult(float a, float B[],float C[])
{
    C[0] = a*B[0];
    C[1] = a*B[1];
    C[2] = a*B[2];
}

CUDA_CALLABLE_MEMBER
void vectorAssign(float a, float b,float c, float D[])
{
    D[0] = a;
    D[1] = b;
    D[2] = c;
}

CUDA_CALLABLE_MEMBER
float vectorNorm(float A[])
{
    float norm = 0.0;
    norm = sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);

        return norm;
}

CUDA_CALLABLE_MEMBER
void vectorDotProduct(float A[], float B[], float C[])
{
    C[0] = A[0]*B[0];
    C[1] = A[1]*B[1];
    C[2] = A[2]*B[2];
}

CUDA_CALLABLE_MEMBER
void vectorCrossProduct(float A[], float B[], float C[])
{
    float tmp[3] = {0.0,0.0,0.0};
    tmp[0] = A[1]*B[2] - A[2]*B[1];
    tmp[1] = A[2]*B[0] - A[0]*B[2];
    tmp[2] = A[0]*B[1] - A[1]*B[0];

    C[0] = tmp[0];
    C[1] = tmp[1];
    C[2] = tmp[2];
}
CUDA_CALLABLE_MEMBER

float getE ( float x0, float y, float z, float E[], Boundary *boundaryVector, int nLines ) {
	float Emag = 0.0;
	float fd = 0.0;
	float pot = 0.0;
    int minIndex = 0;
    float minDistance = 1e12;
    int direction_type;
    float tol = 1e12;
    float point1_dist;
    float point2_dist;
    float perp_dist;
    float directionUnitVector[3] = {0.0,0.0,0.0};
    float vectorMagnitude;
    float max = 0.0;
    float min = 0.0;
    float angle = 0.0;
    float Er = 0.0;
    float Et = 0.0;
#if USECYLSYMM > 0
    float x = sqrt(x0*x0 + y*y);
#else
    float x = x0;
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
    //        std::cout << "p1dist p2dist perpDist " << point1_dist << " " << point2_dist << " " << perp_dist << std::endl;
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
            directionUnitVector[0] = 1.0*sgn((z - boundaryVector[minIndex].intercept_z)/(boundaryVector[minIndex].slope_dzdx) - x0);
            directionUnitVector[1] = 0.0;
            directionUnitVector[2] = 1.0*sgn(perp_dist)/(boundaryVector[minIndex].slope_dzdx);
        //std::cout << "sign boundarVec.slope  sign perp_dist " << sgn(boundaryVector[minIndex].slope_dzdx) << " " << sgn(perp_dist) << std::endl;
        }
        //std::cout << "direction_type 1 " << directionUnitVector[0] << " " << directionUnitVector[1] << " " << directionUnitVector[2] << std::endl;
    }
    else if (direction_type == 2)
    {
        directionUnitVector[0] = (boundaryVector[minIndex].x1 - x);
        directionUnitVector[1] = 0.0;
        directionUnitVector[2] = (boundaryVector[minIndex].z1 - z);
        //std::cout << "direction_type 2 " << directionUnitVector[0] << " " << directionUnitVector[1] << " " << directionUnitVector[2] << std::endl;
    }
    else
    {
        directionUnitVector[0] = (boundaryVector[minIndex].x2 - x);
        directionUnitVector[1] = 0.0;
        directionUnitVector[2] = (boundaryVector[minIndex].z2 - z);
        //std::cout << "direction_type 3 " << directionUnitVector[0] << " " << directionUnitVector[1] << " " << directionUnitVector[2] << std::endl;
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
    //std::cout << "potential and debye length " << pot << " " << boundaryVector[minIndex].debyeLength << " " << pot/boundaryVector[minIndex].debyeLength << std::endl;
        Emag = pot*(fd/(2.0*boundaryVector[minIndex].debyeLength)*exp(-minDistance/(2.0*boundaryVector[minIndex].debyeLength))+ (1.0 - fd)/(boundaryVector[minIndex].larmorRadius)*exp(-minDistance/boundaryVector[minIndex].larmorRadius) );
    if(minDistance == 0.0)
    {
        Emag = 0.0;
        directionUnitVector[0] = 0.0;
        directionUnitVector[1] = 0.0;
        directionUnitVector[2] = 0.0;

    }
        Er = Emag*directionUnitVector[0];
        Et = Emag*directionUnitVector[1];
        E[2] = Emag*directionUnitVector[2];
    //    std::cout << "Emag " << Emag << std::endl;
    //    std::cout << "Min dist " << minDistance << std::endl;
       // std::cout << "r " << x << "z " << z << std::endl;
       // std::cout << "E components " << Er << " " << Et << " " << E[2] << std::endl;
/*        if((Emag == 0) && (minDistance == 0.0)){
        std::cout << "direction Unit vec " << directionUnitVector[0] << " " << directionUnitVector[1] << " " << directionUnitVector[2] << std::endl;
        }*/
    
    //std::cout << "pos " << x << " " << y << " "<< z << " min Dist" << minDistance << "Efield " << Emag << std::endl;
#if USECYLSYMM > 0
            //if cylindrical geometry
            float theta = atan2(y,x0);
  
            E[0] = cos(theta)*Er - sin(theta)*Et;
            E[1] = sin(theta)*Er + cos(theta)*Et;
#else
            E[0] = Er;
            E[1] = Et;
#endif
      //      std::cout << "Ex and Ey " << E[0] << " " << E[1] << std::endl;
    return minDistance;
}

struct move_boris { 
    Boundary *boundaryVector;
    int nR_Bfield;
    int nZ_Bfield;
    float * BfieldGridRDevicePointer;
    float * BfieldGridZDevicePointer;
    float * BfieldRDevicePointer;
    float * BfieldZDevicePointer;
    float * BfieldTDevicePointer;
    int nR_Efield;
    int nZ_Efield;
    float * EfieldGridRDevicePointer;
    float * EfieldGridZDevicePointer;
    float * EfieldRDevicePointer;
    float * EfieldZDevicePointer;
    float * EfieldTDevicePointer;

    const float span;
    const int nLines;
   
    move_boris(float _span, Boundary *_boundaryVector,int _nLines,
            int _nR_Bfield, int _nZ_Bfield,
            float * _BfieldGridRDevicePointer,
            float * _BfieldGridZDevicePointer,
            float * _BfieldRDevicePointer,
            float * _BfieldZDevicePointer,
            float * _BfieldTDevicePointer,
            int _nR_Efield, int _nZ_Efield,
            float * _EfieldGridRDevicePointer,
            float * _EfieldGridZDevicePointer,
            float * _EfieldRDevicePointer,
            float * _EfieldZDevicePointer,
            float * _EfieldTDevicePointer) 
        
        : span(_span), boundaryVector(_boundaryVector), nLines(_nLines), nR_Bfield(_nR_Bfield), 
        nZ_Bfield(_nZ_Bfield), BfieldGridRDevicePointer(_BfieldGridRDevicePointer), 
        BfieldGridZDevicePointer(_BfieldGridZDevicePointer),
        BfieldRDevicePointer(_BfieldRDevicePointer), BfieldZDevicePointer(_BfieldZDevicePointer), 
        BfieldTDevicePointer(_BfieldTDevicePointer),
        nR_Efield(_nR_Efield), nZ_Efield(_nZ_Efield), 
        EfieldGridRDevicePointer(_EfieldGridRDevicePointer), 
        EfieldGridZDevicePointer(_EfieldGridZDevicePointer),
        EfieldRDevicePointer(_EfieldRDevicePointer), EfieldZDevicePointer(_EfieldZDevicePointer),
        EfieldTDevicePointer(_EfieldTDevicePointer) {}

CUDA_CALLABLE_MEMBER    
void operator()(Particle &p) const { 
#ifdef __CUDACC__
#else
float initTime = 0.0;
float interpETime = 0.0;
float interpBTime = 0.0;
float operationsTime = 0.0;
#if USE_BOOST
cpu_timer timer;
cpu_times initTime0 = timer.elapsed();
#endif
#endif
	    if(p.hitWall == 0.0)
        {
            float v_minus[3]= {0.0, 0.0, 0.0};
            float v_prime[3]= {0.0, 0.0, 0.0};
	        float v[3]= {0.0, 0.0, 0.0};
	        float E[3] = {0.0, 0.0, 0.0};
	        float PSE[3] = {0.0, 0.0, 0.0};
	        float B[3] = {0.0,0.0,0.0};
            float br;
            float bz;
            float bt;
	        float dt = span;
	        float Bmag = 0.0;
	        float q_prime = 0.0;
            float coeff = 0.0;
            int nSteps = floor( span / dt + 0.5);
            float minDist = 0.0;
#if ODEINT ==	0  
	        float qpE[3] = {0.0,0.0,0.0};
	        float vmxB[3] = {0.0,0.0,0.0};
	        float vpxB[3] = {0.0,0.0,0.0};
	        float qp_vmxB[3] = {0.0,0.0,0.0};
	        float c_vpxB[3] = {0.0,0.0,0.0};

            for ( int s=0; s<nSteps; s++ ) 
            {
#if USESHEATHEFIELD > 0
	          minDist = getE(p.xprevious,p.yprevious,p.zprevious,E,boundaryVector,nLines);
#endif

#if USEPRESHEATHEFIELD > 0
                 interp2dVector(&PSE[0],p.xprevious,p.yprevious,p.zprevious,nR_Efield,nZ_Efield,
                     EfieldGridRDevicePointer,EfieldGridZDevicePointer,EfieldRDevicePointer,
                     EfieldZDevicePointer,EfieldTDevicePointer);
                 
                 vectorAdd(E,PSE,E);
#endif              
                interp2dVector(&B[0],p.xprevious,p.yprevious,p.zprevious,nR_Bfield,nZ_Bfield,
                    BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,
                    BfieldZDevicePointer,BfieldTDevicePointer);        
            
                Bmag = vectorNorm(B);
	            q_prime = p.charge*1.60217662e-19/(p.amu*1.6737236e-27)*dt*0.5;
                coeff = 2*q_prime/(1+(q_prime*Bmag)*(q_prime*Bmag));
            
                vectorAssign(p.vx,p.vy,p.vz,v);
	            
                //v_minus = v + q_prime*E;
               vectorScalarMult(q_prime,E,qpE);
               vectorAdd(v,qpE,v_minus);

               //v_prime = v_minus + q_prime*(v_minus x B)
                vectorCrossProduct(v_minus,B,vmxB);
                vectorScalarMult(q_prime,vmxB,qp_vmxB);
                vectorAdd(v_minus,qp_vmxB,v_prime);       
                
                //v = v_minus + coeff*(v_prime x B)
                vectorCrossProduct(v_prime, B, vpxB);
                vectorScalarMult(coeff,vpxB,c_vpxB);
                vectorAdd(v_minus, c_vpxB, v);
                
                //v = v + q_prime*E
                vectorAdd(v,qpE,v);
	            
                p.x = p.xprevious + v[0]*dt;
	            p.y = p.yprevious + v[1]*dt;
	            p.z = p.zprevious + v[2]*dt;
	            						
	            p.vx = v[0];
	            p.vy = v[1];
	            p.vz = v[2];
    	    }
#endif

#if ODEINT == 1
        float m = p.amu*1.6737236e-27;
        float q_m = p.charge*1.60217662e-19/m;
        float r[3]= {0.0, 0.0, 0.0};
        float r2[3]= {0.0, 0.0, 0.0};
        float r3[3]= {0.0, 0.0, 0.0};
        float r4[3]= {0.0, 0.0, 0.0};
        float v2[3]= {0.0, 0.0, 0.0};
        float v3[3]= {0.0, 0.0, 0.0};
        float v4[3]= {0.0, 0.0, 0.0};
        float k1r[3]= {0.0, 0.0, 0.0};
        float k2r[3]= {0.0, 0.0, 0.0};
        float k3r[3]= {0.0, 0.0, 0.0};
        float k4r[3]= {0.0, 0.0, 0.0};
        float k1v[3]= {0.0, 0.0, 0.0};
        float k2v[3]= {0.0, 0.0, 0.0};
        float k3v[3]= {0.0, 0.0, 0.0};
        float k4v[3]= {0.0, 0.0, 0.0};
        float dtqm = dt*q_m;
        float vxB[3] = {0.0,0.0,0.0};
        float EplusvxB[3] = {0.0,0.0,0.0};
        float halfKr[3] = {0.0,0.0,0.0};
        float halfKv[3] = {0.0,0.0,0.0};
        float half = 0.5;
                v[0] = p.vx;
                v[1] = p.vy;
	            v[2] = p.vz;

                r[0] = p.xprevious;
                r[1] = p.yprevious;
	            r[2] = p.zprevious;
#ifdef __CUDACC__
#else
#if USE_BOOST
cpu_times initTime1 = timer.elapsed();
initTime = initTime + (initTime1.wall - initTime0.wall);
#endif
#endif
for ( int s=0; s<nSteps; s++ ) 
    {
#ifdef __CUDACC__
#else
#if USE_BOOST
    cpu_times operationsTime0 = timer.elapsed();
    cpu_times interpETime0 = timer.elapsed();
#endif
#endif
#if USESHEATHEFIELD > 0
    minDist = getE(r[0],r[1],r[2],E,boundaryVector,nLines);
#endif
#if USEPRESHEATHEFIELD > 0
    interp2dVector(&PSE[0],p.xprevious,p.yprevious,p.zprevious,nR_Efield,nZ_Efield,
          EfieldGridRDevicePointer,EfieldGridZDevicePointer,EfieldRDevicePointer,
          EfieldZDevicePointer,EfieldTDevicePointer);
                 
    vectorAdd(E,PSE,E);
#endif              
#ifdef __CUDACC__
#else
#if USE_BOOST
    cpu_times interpBTime0 = timer.elapsed();
#endif     
#endif
    interp2dVector(&B[0],r[0],r[1],r[2],nR_Bfield,nZ_Bfield,
               BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,
               BfieldZDevicePointer,BfieldTDevicePointer);        
#ifdef __CUDACC__
#else
#if USE_BOOST
    cpu_times interpBTime1 = timer.elapsed();
    interpETime = interpETime + (interpBTime0.wall - interpETime0.wall);
    interpBTime = interpBTime + (interpBTime1.wall - interpBTime0.wall);
#endif
#endif
    //k1r = dt*v
    vectorScalarMult(dt,v,k1r);
    /*
    k1r[0] = v[0]*dt;
    k1r[1] = v[1]*dt;
    k1r[2] = v[2]*dt;
    */
    //k1v = dt*q_m * (E + (v x B))
    vectorCrossProduct(v,B,vxB);
    vectorAdd(E,vxB,EplusvxB);
    vectorScalarMult(dtqm,EplusvxB,k1v);
    /*
    k1v[0] = dt*q_m*(E[0] + (v[1]*B[2] - v[2]*B[1]));
    k1v[1] = dt*q_m*(E[1] + (v[2]*B[0] - v[0]*B[2]));
    k1v[2] = dt*q_m*(E[2] + (v[0]*B[1] - v[1]*B[0]));
    */
    //r2 = r + 0.5*k1r
    vectorScalarMult(half,k1r,halfKr);
    vectorAdd(r,k1r,r2);
    /*
    r2[0] = r[0] + k1r[0]*0.5;
    r2[1] = r[1] + k1r[1]*0.5;
    r2[2] = r[2] + k1r[2]*0.5;
    */

    //v2 = v + 0.5*k1v
    vectorScalarMult(half,k1v,halfKv);
    vectorAdd(v, halfKv,v2);
        /*
    v2[0] = v[0] + k1v[0]*0.5;
    v2[1] = v[1] + k1v[1]*0.5;
    v2[2] = v[2] + k1v[2]*0.5;
    */
#ifdef __CUDACC__
#else
#if USE_BOOST
    interpETime0 = timer.elapsed();
#endif
#endif

#if USESHEATHEFIELD > 0	  
    minDist = getE(r2[0],r2[1],r2[2],E,boundaryVector,nLines);
#endif
#if USEPRESHEATHEFIELD > 0
    interp2dVector(&PSE[0],p.xprevious,p.yprevious,p.zprevious,nR_Efield,nZ_Efield,
               EfieldGridRDevicePointer,EfieldGridZDevicePointer,EfieldRDevicePointer,
               EfieldZDevicePointer,EfieldTDevicePointer);
    vectorAdd(E,PSE,E);
#endif              
#ifdef __CUDACC__
#else
#if USE_BOOST
    interpBTime0 = timer.elapsed();
#endif
#endif


    interp2dVector(&B[0],r2[0],r2[1],r2[2],nR_Bfield,nZ_Bfield,
             BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,
             BfieldZDevicePointer,BfieldTDevicePointer);        
#ifdef __CUDACC__
#else
#if USE_BOOST
    interpBTime1 = timer.elapsed();
    interpETime = interpETime + (interpBTime0.wall - interpETime0.wall);
    interpBTime = interpBTime + (interpBTime1.wall - interpBTime0.wall);
#endif
#endif
    //k2r = dt*v2
    vectorScalarMult(dt,v2,k2r);
    /*
    k2r[0] = v2[0]*dt;
    k2r[1] = v2[1]*dt;
    k2r[2] = v2[2]*dt;
    */
    //k2v = dt*q_m*(E + (v x B))
    vectorCrossProduct(v2,B,vxB);
    vectorAdd(E,vxB,EplusvxB);
    vectorScalarMult(dtqm,EplusvxB,k2v);
    /*
    k2v[0] = dt*q_m*(E[0] + (v2[1]*B[2] - v2[2]*B[1]));
    k2v[1] = dt*q_m*(E[1] + (v2[2]*B[0] - v2[0]*B[2]));
    k2v[2] = dt*q_m*(E[2] + (v2[0]*B[1] - v2[1]*B[0]));
    */
    //r3 = r + 0.5*k2r
    vectorScalarMult(half,k2r,halfKr);
    vectorAdd(r,k2r,r3);
    /*
    r3[0] = r[0] + k2r[0]*0.5;
    r3[1] = r[1] + k2r[1]*0.5;
    r3[2] = r[2] + k2r[2]*0.5;
    */
    //v3 = v + 0.5*k2v
    vectorScalarMult(half,k2v,halfKv);
    vectorAdd(v, halfKv,v3);
    /*
    v3[0] = v[0] + k2v[0]*0.5;
    v3[1] = v[1] + k2v[1]*0.5;
    v3[2] = v[2] + k2v[2]*0.5;
    */
#ifdef __CUDACC__
#else
#if USE_BOOST
    interpETime0 = timer.elapsed();
#endif
#endif

#if USESHEATHEFIELD > 0	  
    minDist = getE(r3[0],r3[1],r3[2],E,boundaryVector,nLines);
#endif
#if USEPRESHEATHEFIELD > 0
    interp2dVector(&PSE[0],p.xprevious,p.yprevious,p.zprevious,nR_Efield,nZ_Efield,
               EfieldGridRDevicePointer,EfieldGridZDevicePointer,EfieldRDevicePointer,
               EfieldZDevicePointer,EfieldTDevicePointer);
    vectorAdd(E,PSE,E);
#endif              

#ifdef __CUDACC__
#else
#if USE_BOOST
    interpBTime0 = timer.elapsed();
#endif
#endif
    interp2dVector(&B[0],r3[0],r3[1],r3[2],nR_Bfield,nZ_Bfield,
                 BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,
                 BfieldZDevicePointer,BfieldTDevicePointer);        
                
#ifdef __CUDACC__
#else
#if USE_BOOST
    interpBTime1 = timer.elapsed();
    interpETime = interpETime + (interpBTime0.wall - interpETime0.wall);
    interpBTime = interpBTime + (interpBTime1.wall - interpBTime0.wall);
#endif
#endif
    //k3r = dt*v3
    vectorScalarMult(dt,v3,k3r);
    /*
    k3r[0] = v3[0]*dt;
    k3r[1] = v3[1]*dt;
    k3r[2] = v3[2]*dt;
    */
    //k3v = dt*qm*(E + (v x B))
    vectorCrossProduct(v3,B,vxB);
    vectorAdd(E,vxB,EplusvxB);
    vectorScalarMult(dtqm,EplusvxB,k3v);
    /*
    k3v[0] = dt*q_m*(E[0] + (v3[1]*B[2] - v3[2]*B[1]));
    k3v[1] = dt*q_m*(E[1] + (v3[2]*B[0] - v3[0]*B[2]));
    k3v[2] = dt*q_m*(E[2] + (v3[0]*B[1] - v3[1]*B[0]));
    */
    //r4 = r + k3r
    vectorAdd(r, k3r,r4);
    /*
    r4[0] = r[0] + k3r[0];
    r4[1] = r[1] + k3r[1];
    r4[2] = r[2] + k3r[2];
    */
    //v4 = v + k3v
    vectorAdd(v, k3v, v4);
        /*
    v4[0] = v[0] + k3v[0];
    v4[1] = v[1] + k3v[1];
    v4[2] = v[2] + k3v[2];
    */
#ifdef __CUDACC__
#else
#if USE_BOOST
    interpETime0 = timer.elapsed();
#endif
#endif

#if USESHEATHEFIELD > 0            
	minDist = getE(r4[0],r4[1],r4[2],E,boundaryVector,nLines);
#endif
#if USEPRESHEATHEFIELD > 0
   interp2dVector(&PSE[0],p.xprevious,p.yprevious,p.zprevious,nR_Efield,nZ_Efield,
               EfieldGridRDevicePointer,EfieldGridZDevicePointer,EfieldRDevicePointer,
               EfieldZDevicePointer,EfieldTDevicePointer);
    vectorAdd(E,PSE,E);
#endif              
#ifdef __CUDACC__
#else
#if USE_BOOST
    interpBTime0 = timer.elapsed();
#endif
#endif

    interp2dVector(&B[0],r4[0],r4[1],r4[2],nR_Bfield,nZ_Bfield,
                        BfieldGridRDevicePointer,BfieldGridZDevicePointer,
                        BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer);        
#ifdef __CUDACC__
#else
#if USE_BOOST
    interpBTime1 = timer.elapsed();
    interpETime = interpETime + (interpBTime0.wall - interpETime0.wall);
    interpBTime = interpBTime + (interpBTime1.wall - interpBTime0.wall);
#endif
#endif

    //k4r = dt*v4
    vectorScalarMult(dt,v4,k4r);
    /*
   k4r[0] = v4[0]*dt;
   k4r[1] = v4[1]*dt;
   k4r[2] = v4[2]*dt;
   */
    //k4v = dt*q_m*(E + (v x B))
    vectorCrossProduct(v4,B,vxB);
    vectorAdd(E,vxB,EplusvxB);
    vectorScalarMult(dtqm,EplusvxB,k4v);
    /*
   k4v[0] = dt*q_m*(E[0] + (v4[1]*B[2] - v4[2]*B[1]));
   k4v[1] = dt*q_m*(E[1] + (v4[2]*B[0] - v4[0]*B[2]));
   k4v[2] = dt*q_m*(E[2] + (v4[0]*B[1] - v4[1]*B[0]));
   */
   p.x = r[0] + (k1r[0] + 2*k2r[0] + 2*k3r[0] + k4r[0])/6;
   p.y = r[1] + (k1r[1] + 2*k2r[1] + 2*k3r[1] + k4r[1])/6;
   p.z = r[2] + (k1r[2] + 2*k2r[2] + 2*k3r[2] + k4r[2])/6;
   p.vx = v[0] + (k1v[0] + 2*k2v[0] + 2*k3v[0] + k4v[0])/6;
   p.vy = v[1] + (k1v[1] + 2*k2v[1] + 2*k3v[1] + k4v[1])/6;
   p.vz = v[2] + (k1v[2] + 2*k2v[2] + 2*k3v[2] + k4v[2])/6;
#ifdef __CUDACC__
#else
#if USE_BOOST
    cpu_times operationsTime1 = timer.elapsed();
    operationsTime = operationsTime1.wall - operationsTime0.wall;
#endif
#endif
//std::cout << "Operations Time: " << operationsTime <<std::endl;
//std::cout << "Efield Interpolation Time: " << interpETime <<std::endl;
//std::cout << "Bfield Interpolation Time: " << interpBTime <<std::endl;
//std::cout << "Init Time: " << initTime <<std::endl;
            }
#endif
        }
    } 
};

#endif
