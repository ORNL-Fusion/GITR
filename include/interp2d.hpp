#ifndef _INTERP2D_
#define _INTERP2D_

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
using namespace std;

CUDA_CALLABLE_MEMBER

double interp2dCombined ( double x, double y, double z,int nx, int nz,
//#ifdef __CUDACC__
double* gridx,double* gridz,double* data ) {
/*#else        
std::vector<double>* gridxp,std::vector<double>* gridzp,std::vector<double>* datap ) {
    std::vector<double>& data = *datap;
    std::vector<double>& gridx = *gridxp;
    std::vector<double>& gridz = *gridzp;
#endif
*/
    double fxz = 0.0;
    if(nx*nz == 1)
    {
        fxz = data[0];
    }
    else{
#if USECYLSYMM > 0
    double dim1 = sqrt(x*x + y*y);
#else
    double dim1 = x;
#endif    
    double d_dim1 = gridx[1] - gridx[0];
    double dz = gridz[1] - gridz[0];
    int i = floor((dim1 - gridx[0])/d_dim1 );//addition of 0.5 finds nearest gridpoint
    int j = floor((z - gridz[0])/dz);
    double interp_value = data[i + j*nz];
    double fx_z1 = ((gridx[i+1]-dim1)*data[i+j*nx] + (dim1 - gridx[i])*data[i+1+j*nx])/d_dim1;
    double fx_z2 = ((gridx[i+1]-dim1)*data[i+(j+1)*nx] + (dim1 - gridx[i])*data[i+1+(j+1)*nx])/d_dim1; 
    fxz = ((gridz[j+1]-z)*fx_z1+(z - gridz[j])*fx_z2)/dz;
    }

    return fxz;
}
void interp2dVector (double* field, double x, double y, double z,int nx, int nz,
double* gridx,double* gridz,double* datar, double* dataz, double* datat ) {

   double Ar = interp2dCombined(x,y,z,nx,nz,gridx,gridz, datar);
   double At = interp2dCombined(x,y,z,nx,nz,gridx,gridz, datat);
   field[2] = interp2dCombined(x,y,z,nx,nz,gridx,gridz, dataz);

#if USECYLSYMM > 0
            //if cylindrical geometry
            double theta = atan(y/x);
            if(x < 0.0)
            {
                if(y > 0.0)
                {
                    theta = theta + 3.1415;
                }
                else
                {
                    theta = sqrt(theta*theta) + 3.1415;
                }
            }
  
            field[0] = cos(theta)*Ar - sin(theta)*At;
            field[1] = sin(theta)*Ar + cos(theta)*At;
#else
            field[0] = Ar;
            field[1] = At;
#endif

}
CUDA_CALLABLE_MEMBER_DEVICE

double interp2dcuda ( double x, double y, double z,int nx, int nz, double* gridx,
        double* gridz,double* data ) {

#if USECYLSYMM > 0
    double dim1 = sqrt(x*x + y*y);
#else
    double dim1 = x;
#endif    
    double d_dim1 = gridx[1] - gridx[0];
    double dz = gridz[1] - gridz[0];
    int i = floor((dim1 - gridx[0])/d_dim1 );//addition of 0.5 finds nearest gridpoint
    int j = floor((z - gridz[0])/dz);
    double interp_value = data[i + j*nz];
    double fx_z1 = ((gridx[i+1]-dim1)*data[i+j*nx] + (dim1 - gridx[i])*data[i+1+j*nx])/d_dim1;
    double fx_z2 = ((gridx[i+1]-dim1)*data[i+(j+1)*nx] + (dim1 - gridx[i])*data[i+1+(j+1)*nx])/d_dim1; 
    double fxz = ((gridz[j+1]-z)*fx_z1+(z - gridz[j])*fx_z2)/dz;
    return fxz;
}

double interp1dUnstructured(double samplePoint,int nx, double max_x, double* data)
{
    int done = 0;
    int low_index = 0;
    double interpolated_value = 0.0;

    for(int i=0;i<nx;i++)
    {
        if(done == 0)
        {
            if(samplePoint > data[i])
            {
                done = 1;
                low_index = i;
            }   
        }
    }

    interpolated_value = low_index*max_x/nx + max_x/nx*(samplePoint - data[low_index])/(data[low_index+1] - data[low_index]);
    return interpolated_value;
}
#endif

