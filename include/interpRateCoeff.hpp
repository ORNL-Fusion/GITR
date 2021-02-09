#ifndef _INTERPRATECOEFF2D_
#define _INTERPRATECOEFF2D_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "interp2d.hpp"
#include <thrust/device_vector.h>
#include <vector>
#include <cmath>
using namespace std;

CUDA_CALLABLE_MEMBER


float rateCoeffInterp(int charge, float te, float ne,int nT, int nD, float* rateGrid_Tempp,float* rateGrid_Densp,float* Ratesp){

/*    std::vector<float>& rateGrid_Temp = *rateGrid_Tempp;
    std::vector<float>& rateGrid_Dens = *rateGrid_Densp;
    std::vector<float>& Rates = *Ratesp;
  */  
    int indT = 0;
    int indN = 0;
    float logT = std::log10(te);
    float logn = std::log10(ne);
    //std::cout << "Rategrid_temp in rateCoeffInterp " << rateGrid_Temp[1] << std::endl;
    float d_T = rateGrid_Tempp[1] - rateGrid_Tempp[0];
    float d_n = rateGrid_Densp[1] - rateGrid_Densp[0];
   // if (logT >= rateGrid_Tempp[0] && logT <= rateGrid_Tempp[nT-2])
   // {
        indT = std::floor((logT - rateGrid_Tempp[0])/d_T );//addition of 0.5 finds nearest gridpoint
    //}
    //if (logn >= rateGrid_Densp[0] && logn <= rateGrid_Densp[nD-2])
    //{
        indN = std::floor((logn - rateGrid_Densp[0])/d_n );
    //}
    //std::cout << "Indices density, temp " << indN << " " <<indT<<std::endl;
    //std::cout << "charge " << charge << std::endl;
    //std::cout << "Lower value " << Ratesp[charge*nT*nD + indT*nD + indN] << std::endl;
if(indT < 0 || indT > nT-2)
{indT = 0;}
if(indN < 0 || indN > nD-2)
{indN = 0;}
if(charge > 74-1)
{charge = 0;}
        float aT = std::pow(10.0f,rateGrid_Tempp[indT+1]) - te;
    float bT = te - std::pow(10.0f,rateGrid_Tempp[indT]);
    float abT = aT+bT;

    float aN = std::pow(10.0f,rateGrid_Densp[indN+1]) - ne;
    float bN = ne - std::pow(10.0f, rateGrid_Densp[indN]);
    float abN = aN + bN;

    //float interp_value = Rates[charge*rateGrid_Temp.size()*rateGrid_Dens.size()            + indT*rateGrid_Dens.size() + indN];

    float fx_z1 = (aN*std::pow(10.0f,Ratesp[charge*nT*nD + indT*nD + indN]) 
            + bN*std::pow(10.0f,Ratesp[charge*nT*nD            + indT*nD + indN + 1]))/abN;
    
    float fx_z2 = (aN*std::pow(10.0f,Ratesp[charge*nT*nD            + (indT+1)*nD + indN]) 
            + bN*std::pow(10.0f,Ratesp[charge*nT*nD            + (indT+1)*nD + indN+1]))/abN;
    float fxz = (aT*fx_z1+bT*fx_z2)/abT;
    //std::cout << "fxz1 and 2 " << fx_z1 << " " << fx_z2<< " "<< fxz << std::endl;
    //printf ("floats:%i %e  \n",charge,fxz);
    return fxz;    
}

CUDA_CALLABLE_MEMBER
float interpRateCoeff2d ( int charge, float x, float y, float z,int nx, int nz, float* tempGridxp,
       float* tempGridzp, float* Tempp,
       float* densGridxp,float* densGridzp,float* Densp,int nT_Rates, int nD_Rates,
       float* rateGrid_Temp,float* rateGrid_Dens,float* Rates ) {
//    std::cout << "rate test " << Tempp[0] << std::endl;
    /*std::vector<float>& Tdata = *Tempp;
    std::vector<float>& Tgridx = *tempGridxp;
    std::vector<float>& Tgridz = *tempGridzp;
    std::vector<float>& DensityData = *Densp;
    std::vector<float>& DensGridx = *densGridxp;
    std::vector<float>& DensGridz = *densGridzp;
*/
    //std::cout << "at tlocal interp routine " <<x << y << z<< " " << nx << nz<< std::endl;
    //std::cout << "Interpolating local temp at "<<x << " " << y << " " << z << std::endl;
    float tlocal = interp2dCombined(x,y,z,nx,nz,tempGridxp,tempGridzp,Tempp);
    //std::cout << "Interpolating local dens " << std::endl;
    float nlocal = interp2dCombined(x,y,z,nx,nz,densGridxp,densGridzp,Densp);
    //std::cout << "tlocal" << tlocal << std::endl;
    //std::cout << "nlocal" << nlocal << std::endl;
    //std::cout << "Interpolating RC " << std::endl;
    float RClocal = rateCoeffInterp(charge,tlocal,nlocal,nT_Rates,nD_Rates,rateGrid_Temp, rateGrid_Dens, Rates);
    float tion = 1/(RClocal*nlocal);
    if(tlocal == 0.0 || nlocal == 0.0) tion=1.0e12;
    //std::cout << "Returning " << std::endl;
    return tion;
}

#endif

