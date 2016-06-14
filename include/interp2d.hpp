#ifndef _INTERP2D_
#define _INTERP2D_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include <vector>
#include <math.h>
using namespace std;

//CUDA_CALLABLE_MEMBER

double interp2d ( double x, double y, double z, std::vector<double> gridx,std::vector<double> gridz,
                    std::vector<double> data ) {
#if USECYLSYMM > 0
    double dim1 = sqrt(x*x + y*y);
#else
    double dim1 = x;
#endif    
    double d_dim1 = gridx[1] - gridx[0];
    double dz = gridz[1] - gridz[0];
    int i = floor((dim1 - gridx[0])/d_dim1 );//addition of 0.5 finds nearest gridpoint
    int j = floor((z - gridz[0])/dz);
    std::cout << "indices from interpolation " << i << " " << j << std::endl;
    double interp_value = data[i + j*gridz.size()];
    double fx_z1 = ((gridx[i+1]-dim1)*data[i+j*gridx.size()] + (dim1 - gridx[i])*data[i+1+j*gridx.size()])/d_dim1;
    double fx_z2 = ((gridx[i+1]-dim1)*data[i+(j+1)*gridx.size()] + (dim1 - gridx[i])*data[i+1+(j+1)*gridx.size()])/d_dim1; 
    std::cout << "four grid values" << data[ i+ j*gridx.size()] << " " << data[i+1+j*gridx.size()] << " " << data[i+(j+1)*gridx.size()] << "  " << data[i+1+(j+1)*gridx.size()] << std::endl;
    std::cout << "fx_zs " << fx_z1 << " " << fx_z2 << std::endl;
    double fxz = ((gridz[j+1]-z)*fx_z1+(z - gridz[j])*fx_z2)/dz;
    std::cout << "interp_value in function: " << interp_value << std::endl;
    return fxz;
}
#endif


