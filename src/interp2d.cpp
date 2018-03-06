#if USE_BOOST
#endif
#include "interp2d.hpp"
#include <iostream>
#ifdef __CUDACC__
#endif

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

using namespace std;
CUDA_CALLABLE_MEMBER

float interp2d ( float x, float z,int nx, int nz,
    float* gridx,float* gridz,float* data ) {
    
    float fxz = 0.0;
    float fx_z1 = 0.0;
    float fx_z2 = 0.0; 
    if(nx*nz == 1)
    {
        fxz = data[0];
    }
    else{
    float dim1 = x;
    float d_dim1 = gridx[1] - gridx[0];
    float dz = gridz[1] - gridz[0];
    int i = floor((dim1 - gridx[0])/d_dim1);//addition of 0.5 finds nearest gridpoint
    int j = floor((z - gridz[0])/dz);
    
    //float interp_value = data[i + j*nx];
    if (i < 0) i =0;
    if (j< 0 ) j=0;
    if (i >=nx-1 && j>=nz-1)
    {
        fxz = data[nx-1+(nz-1)*nx];
    }
    else if (i >=nx-1)
    {
        fx_z1 = data[nx-1+j*nx];
        fx_z2 = data[nx-1+(j+1)*nx];
        fxz = ((gridz[j+1]-z)*fx_z1+(z - gridz[j])*fx_z2)/dz;
    }
    else if (j >=nz-1)
    {
        fx_z1 = data[i+(nz-1)*nx];
        fx_z2 = data[i+(nz-1)*nx];
        fxz = ((gridx[i+1]-dim1)*fx_z1+(dim1 - gridx[i])*fx_z2)/d_dim1;
        
    }
    else
    {
      fx_z1 = ((gridx[i+1]-dim1)*data[i+j*nx] + (dim1 - gridx[i])*data[i+1+j*nx])/d_dim1;
      fx_z2 = ((gridx[i+1]-dim1)*data[i+(j+1)*nx] + (dim1 - gridx[i])*data[i+1+(j+1)*nx])/d_dim1; 
      fxz = ((gridz[j+1]-z)*fx_z1+(z - gridz[j])*fx_z2)/dz;
    }
    }

    return fxz;
}
float interp2dCombined ( float x, float y, float z,int nx, int nz,
    float* gridx,float* gridz,float* data ) {
    
    float fxz = 0.0;
    float fx_z1 = 0.0;
    float fx_z2 = 0.0; 
    if(nx*nz == 1)
    {
        fxz = data[0];
    }
    else{
#if USECYLSYMM > 0
    float dim1 = sqrt(x*x + y*y);
#else
    float dim1 = x;
#endif    
    float d_dim1 = gridx[1] - gridx[0];
    float dz = gridz[1] - gridz[0];
    int i = floor((dim1 - gridx[0])/d_dim1);//addition of 0.5 finds nearest gridpoint
    int j = floor((z - gridz[0])/dz);
    
    //float interp_value = data[i + j*nx];
    if (i < 0) i=0;
    if (j < 0) j=0;
    if (i >=nx-1 && j>=nz-1)
    {
        fxz = data[nx-1+(nz-1)*nx];
    }
    else if (i >=nx-1)
    {
        fx_z1 = data[nx-1+j*nx];
        fx_z2 = data[nx-1+(j+1)*nx];
        fxz = ((gridz[j+1]-z)*fx_z1+(z - gridz[j])*fx_z2)/dz;
    }
    else if (j >=nz-1)
    {
        fx_z1 = data[i+(nz-1)*nx];
        fx_z2 = data[i+(nz-1)*nx];
        fxz = ((gridx[i+1]-dim1)*fx_z1+(dim1 - gridx[i])*fx_z2)/d_dim1;
        
    }
    else
    {
      fx_z1 = ((gridx[i+1]-dim1)*data[i+j*nx] + (dim1 - gridx[i])*data[i+1+j*nx])/d_dim1;
      fx_z2 = ((gridx[i+1]-dim1)*data[i+(j+1)*nx] + (dim1 - gridx[i])*data[i+1+(j+1)*nx])/d_dim1; 
      fxz = ((gridz[j+1]-z)*fx_z1+(z - gridz[j])*fx_z2)/dz;
    }
    }

    return fxz;
}

CUDA_CALLABLE_MEMBER

float interp3d ( float x, float y, float z,int nx,int ny, int nz,
    float* gridx,float* gridy, float* gridz,float* data ) {
    
    float fxyz = 0.0;

    float dx = gridx[1] - gridx[0];
    float dy = gridy[1] - gridy[0];
    float dz = gridz[1] - gridz[0];
    
    int i = floor((x - gridx[0])/dx);//addition of 0.5 finds nearest gridpoint
    int j = floor((y - gridy[0])/dy);
    int k = floor((z - gridz[0])/dz);
    //std::cout << "dxyz ijk " << dx << " "<<dy << " " << dz<< " " << i
    //    << " " << j << " " << k << std::endl;
    if(i <0 ) i=0;
    else if(i >=nx-1) i=nx-2;
    if(j <0 ) j=0;
    else if(j >=ny-1) j=ny-2;
    if(k <0 ) k=0;
    else if(k >=nz-1) k=nz-2;
    if(ny <=1) j=0;
    if(nz <=1) k=0;
    //if(j <0 || j>ny-1) j=0;
    //if(k <0 || k>nz-1) k=0;
    float fx_z0 = (data[i + j*nx + k*nx*ny]*(gridx[i+1]-x) + data[i +1 + j*nx + k*nx*ny]*(x-gridx[i]))/dx;
    float fx_z1 = (data[i + j*nx + (k+1)*nx*ny]*(gridx[i+1]-x) + data[i +1 + j*nx + (k+1)*nx*ny]*(x-gridx[i]))/dx;

    
    //std::cout << "fxz0 fxz1 " << fx_z0 << " "<<fx_z1 << std::endl;
    float fxy_z0 = (data[i + (j+1)*nx + k*nx*ny]*(gridx[i+1]-x) + data[i +1 + (j+1)*nx + k*nx*ny]*(x-gridx[i]))/dx;
    float fxy_z1 = (data[i + (j+1)*nx + (k+1)*nx*ny]*(gridx[i+1]-x) + data[i +1 + (j+1)*nx + (k+1)*nx*ny]*(x-gridx[i]))/dx;
    //std::cout << "fxyz0 fxyz1 " << fxy_z0 << " "<<fxy_z1 << std::endl;

    float fxz0 = (fx_z0*(gridz[k+1] - z) + fx_z1*(z-gridz[k]))/dz;
    float fxz1 = (fxy_z0*(gridz[k+1] - z) + fxy_z1*(z-gridz[k]))/dz;
    //std::cout << "fxz0 fxz1 " << fxz0 << " "<<fxz1 << std::endl;

    fxyz = (fxz0*(gridy[j+1] - y) + fxz1*(y-gridy[j]))/dy;
    if(ny <=1) fxyz=fxz0;
    if(nz <=1) fxyz=fx_z0;
    //std::cout <<"fxyz " << fxyz << std::endl;
    return fxyz;
}

CUDA_CALLABLE_MEMBER
void interp3dVector (float* field, float x, float y, float z,int nx,int ny, int nz,
        float* gridx,float* gridy,float* gridz,float* datar, float* dataz, float* datat ) {

    field[0] =  interp3d (x,y,z,nx,ny,nz,gridx, gridy,gridz,datar );
    field[1] =  interp3d (x,y,z,nx,ny,nz,gridx, gridy,gridz,datat );
    field[2] =  interp3d (x,y,z,nx,ny,nz,gridx, gridy,gridz,dataz );
}
CUDA_CALLABLE_MEMBER
void interp2dVector (float* field, float x, float y, float z,int nx, int nz,
float* gridx,float* gridz,float* datar, float* dataz, float* datat ) {

   float Ar = interp2dCombined(x,y,z,nx,nz,gridx,gridz, datar);
   float At = interp2dCombined(x,y,z,nx,nz,gridx,gridz, datat);
   field[2] = interp2dCombined(x,y,z,nx,nz,gridx,gridz, dataz);
#if USECYLSYMM > 0
            float theta = atan2f(y,x);   
            field[0] = cosf(theta)*Ar - sinf(theta)*At;
            field[1] = sinf(theta)*Ar + cosf(theta)*At;
#else
            field[0] = Ar;
            field[1] = At;
#endif

}
CUDA_CALLABLE_MEMBER
float interp1dUnstructured(float samplePoint,int nx, float max_x, float* data,int &lowInd)
{
    int done = 0;
    int low_index = 0;
    float interpolated_value = 0.0;

    for(int i=0;i<nx;i++)
    {
        if(done == 0)
        {
            if(samplePoint < data[i])
            {
                done = 1;
                low_index = i-1;
            }   
        }
    }
    //std::cout << " lowInd " << low_index << " " << data[low_index] << " " << data[low_index+1] <<
        //std::endl;
    interpolated_value =
        ((data[low_index+1] - samplePoint)*low_index*max_x/nx
        + (samplePoint - data[low_index])*(low_index+1)*max_x/nx)/(data[low_index+1]- data[low_index]);
      //(low_index+1)*max_x/nx
    lowInd = low_index;
    if(low_index < 0)
    {
        //std::cout << "WARNING: interpolated value is outside range of CDF " << std::endl;
        lowInd = 0;
        //interpolated_value = samplePoint*data[0];
        if(samplePoint > 0.0)
        {
          interpolated_value = samplePoint;
        }
        else{
          interpolated_value = 0.0;
        }
    }
    return interpolated_value;
}
CUDA_CALLABLE_MEMBER
float interp1dUnstructured2(float samplePoint,int nx, float *xdata, float* data)
{
    int done = 0;
    int low_index = 0;
    float interpolated_value = 0.0;

    for(int i=0;i<nx;i++)
    {
        if(done == 0)
        {
            if(samplePoint < data[i])
            {
                done = 1;
                low_index = i-1;
            }   
        }
    }
    interpolated_value =    xdata[low_index]  
        + (samplePoint - data[low_index])*(xdata[low_index + 1] - xdata[low_index])/(data[low_index+1]- data[low_index]);
    return interpolated_value;
}
CUDA_CALLABLE_MEMBER
float interp2dUnstructured(float x,float y,int nx,int ny, float *xgrid,float *ygrid, float* data)
{
    int doneX = 0;
    int doneY = 0;
    int xInd = 0;
    float xDiffUp = 0.0;
    float xDiffDown = 0.0;
    int yInd = 0;
    float dx;
    float yLowValue; 
    float yHighValue;
    float yDiffUp;
    float yDiffDown; 
    float dy;
    float fxy=0.0;
    float factor = 1.0;

    if(x >= xgrid[0] && x<= xgrid[nx-1])
    {
      for(int i=0;i<nx;i++)
      {
          if(!doneX)
          {
             if(x<xgrid[i])
               {
                  doneX = 1;
                  xInd = i-1;
               }
          }
      }
    }
    else
    {
        factor = 0.0;
    }
    if(y >= ygrid[0] && y<= ygrid[ny-1])
    {
      for(int i=0;i<ny;i++)
      {
          if(!doneY)
          {
             if(y<ygrid[i])
               {
                  doneY = 1;
                  yInd = i-1;
               }
          }
      }
    }
    else
    {
        factor = 0.0;
    }
   
    //std::cout << "x vals " << xgrid[xInd] << " " << xgrid[xInd+1];
    //std::cout << "y vals " << ygrid[yInd] << " " << ygrid[yInd+1];
    xDiffUp = xgrid[xInd+1] - x;
    xDiffDown = x-xgrid[xInd];
    dx = xgrid[xInd+1]-xgrid[xInd];
    //std::cout << "dx, data vals " << dx << " " << data[xInd + yInd*nx] << " " <<
                 //data[xInd+1 + yInd*nx] << " " << data[xInd + (yInd+1)*nx] << " " << 
                 //data[xInd+1 + (yInd+1)*nx] << std::endl;
    yLowValue = (xDiffUp*data[xInd + yInd*nx] + xDiffDown*data[xInd+1 + yInd*nx])/dx;
    yHighValue = (xDiffUp*data[xInd + (yInd+1)*nx] + xDiffDown*data[xInd+1 + (yInd+1)*nx])/dx;
    yDiffUp = ygrid[yInd+1]-y;
    yDiffDown = y - ygrid[yInd];
    dy = ygrid[yInd+1] - ygrid[yInd];
    fxy = factor*(yDiffUp*yLowValue + yDiffDown*yHighValue)/dy;

    return fxy;

}
