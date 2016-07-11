#ifndef _INTERPRATECOEFF2D_
#define _INTERPRATECOEFF2D_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include <thrust/device_vector.h>
#include <vector>
#include <math.h>
#include <cmath>
using namespace std;

//CUDA_CALLABLE_MEMBER


double rateCoeffInterp(int charge, double te, double ne,int nT, int nD, double* rateGrid_Tempp,double* rateGrid_Densp,double* Ratesp){

/*    std::vector<double>& rateGrid_Temp = *rateGrid_Tempp;
    std::vector<double>& rateGrid_Dens = *rateGrid_Densp;
    std::vector<double>& Rates = *Ratesp;
  */  
    int indT = 0;
    int indN = 0;
    double logT = log10(te);
    double logn = log10(ne);
    //std::cout << "Rategrid_temp in rateCoeffInterp " << rateGrid_Temp[1] << std::endl;
    double d_T = rateGrid_Tempp[1] - rateGrid_Tempp[0];
    double d_n = rateGrid_Densp[1] - rateGrid_Densp[0];
    if (logT >= rateGrid_Tempp[0])
    {
        indT = floor((logT - rateGrid_Tempp[0])/d_T );//addition of 0.5 finds nearest gridpoint
    }
    if (logn >= rateGrid_Densp[0])
    {
        indN = floor((logn - rateGrid_Densp[0])/d_n );
    }
    std::cout << "Indices density, temp " << indN << " " <<indT<<std::endl;
    std::cout << "charge " << charge << std::endl;
    std::cout << "Lower value " << Ratesp[charge*nT*nD + indT*nD + indN] << std::endl;
    double aT = pow(10.0,rateGrid_Tempp[indT+1]) - te;
    double bT = te - pow(10.0,rateGrid_Tempp[indT]);
    double abT = aT+bT;

    double aN = pow(10.0,rateGrid_Densp[indN+1]) - ne;
    double bN = ne - pow(10.0, rateGrid_Densp[indN]);
    double abN = aN + bN;

    //double interp_value = Rates[charge*rateGrid_Temp.size()*rateGrid_Dens.size()            + indT*rateGrid_Dens.size() + indN];

    double fx_z1 = (aN*pow(10.0,Ratesp[charge*nT*nD + indT*nD + indN]) 
            + bN*pow(10.0,Ratesp[charge*nT*nD            + indT*nD + indN + 1]))/abN;
    
    double fx_z2 = (aN*pow(10.0,Ratesp[charge*nT*nD            + (indT+1)*nD + indN]) 
            + bN*pow(10.0,Ratesp[charge*nT*nD            + (indT+1)*nD + indN+1]))/abN;
    double fxz = (aT*fx_z1+bT*fx_z2)/abT;
    std::cout << "fxz1 and 2 " << fx_z1 << " " << fx_z2<< " "<< fxz << std::endl;
    return fxz;    
}

double interpRateCoeff2d ( int charge, double x, double y, double z,int nx, int nz, double* tempGridxp,
       double* tempGridzp, double* Tempp,
       double* densGridxp,double* densGridzp,double* Densp,int nT_Rates, int nD_Rates,
       double* rateGrid_Temp,double* rateGrid_Dens,double* Rates ) {
//    std::cout << "rate test " << Tempp[0] << std::endl;
    /*std::vector<double>& Tdata = *Tempp;
    std::vector<double>& Tgridx = *tempGridxp;
    std::vector<double>& Tgridz = *tempGridzp;
    std::vector<double>& DensityData = *Densp;
    std::vector<double>& DensGridx = *densGridxp;
    std::vector<double>& DensGridz = *densGridzp;
*/
    double tlocal = interp2dCombined(x,y,z,nx,nz,tempGridxp,tempGridzp,Tempp);
    double nlocal = interp2dCombined(x,y,z,nx,nz,densGridxp,densGridzp,Densp);
    std::cout << "tlocal" << tlocal << std::endl;
    std::cout << "nlocal" << nlocal << std::endl;
    double RClocal = rateCoeffInterp(charge,tlocal,nlocal,nT_Rates,nD_Rates,rateGrid_Temp, rateGrid_Dens, Rates);
    return RClocal;
}

#endif

