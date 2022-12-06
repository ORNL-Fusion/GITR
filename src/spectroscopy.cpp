#include "spectroscopy.h"

#if USE_CUDA >0
__device__ double atomicAdd1(gitr_precision* address, gitr_precision val)
{
    unsigned long long int* address_as_ull =
                        (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
      do {
             assumed = old;
             old = atomicCAS(address_as_ull, assumed,
                            __double_as_longlong(val + 
                                __longlong_as_double(assumed)));
                 // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
         } while (assumed != old);
                 
         return __longlong_as_double(old);
}
#endif

spec_bin::spec_bin(Flags* _flags, Particles *_particlesPointer, int _nBins,int _nX,int _nY, int _nZ, gitr_precision *_gridX,gitr_precision *_gridY,gitr_precision *_gridZ,
           gitr_precision* _bins, gitr_precision _dt, int cylsymm_, int spectroscopy_ ) : 
        flags(_flags), particlesPointer(_particlesPointer), nBins(_nBins),nX(_nX),nY(_nY), nZ(_nZ), gridX(_gridX),gridY(_gridY),gridZ(_gridZ), bins(_bins),
        dt(_dt), cylsymm( cylsymm_ ), spectroscopy( spectroscopy_ ) {}

CUDA_CALLABLE_MEMBER_DEVICE    
void spec_bin::operator()(std::size_t indx) const {
    gitr_precision dx = 0.0;
    gitr_precision dy = 0.0;
    gitr_precision dz = 0.0;
    gitr_precision x = particlesPointer->xprevious[indx];
    gitr_precision y = particlesPointer->yprevious[indx];
    gitr_precision z = particlesPointer->zprevious[indx];
    gitr_precision dt_particle = 0.0;

    gitr_precision dim1;
    if( spectroscopy > 2 )
    {
    dim1 = particlesPointer->xprevious[indx];
    }
    else
    {
     if( cylsymm > 0 )
     {
    dim1 = std::sqrt(x*x + y*y);
    }
    else
    {
    dim1 = x;
    }
    }

    if ((z > gridZ[0]) && (z < gridZ[nZ-1]))
        {
          if((dim1 > gridX[0]) && (dim1 < gridX[nX-1]))
          {
              dx = gridX[1] - gridX[0];
              dz = gridZ[1] - gridZ[0];

              int indx_X;
              int indx_Z;
              int indx_Y;
              int nnYY=1;

              if( spectroscopy < 3 )
              {
                indx_X = std::floor((dim1-gridX[0])/dx);
                indx_Z = std::floor((z-gridZ[0])/dz);
                indx_Y = 0;
                nnYY=1;
              if (indx_X < 0 || indx_X >= nX) indx_X = 0;

              if (indx_Z < 0 || indx_Z >= nZ) indx_Z = 0;

              int charge = std::floor(particlesPointer->charge[indx]);

              gitr_precision specWeight = 0.0;

              if(particlesPointer->hitWall[indx]== 0.0)
              {
                if (flags->USE_ADAPTIVE_DT) 
                {
	                if(particlesPointer->advance[indx])
                  {
                    dt_particle = particlesPointer->dt[indx];
                    specWeight = particlesPointer->weight[indx]*dt_particle/dt;
		              }
                }

                else
                {
                  specWeight = particlesPointer->weight[indx];
                }
#if USE_CUDA >0

	              int index = nBins*nX*nnYY*nZ + indx_Z*nX*nnYY +indx_Y*nX+ indx_X;

                atomicAdd1(&bins[index], specWeight);//0*nX*nZ + indx_Z*nZ + indx_X

              if(charge < nBins)
              {
                atomicAdd1( &bins[charge*nX*nnYY*nZ + indx_Z*nX*nnYY + indx_Y*nX+ indx_X],
                            1.0*specWeight);
              }

#else

              #pragma omp atomic
              bins[nBins*nX*nnYY*nZ + indx_Z*nX*nnYY  +indx_Y*nX +indx_X] += specWeight;

              if(charge < nBins)
              {
                #pragma omp atomic
                bins[charge*nX*nnYY*nZ + indx_Z*nX*nnYY +indx_Y*nX + indx_X] += specWeight;
              }
#endif
            }
              }

              else if((y > gridY[0]) && (y < gridY[nY-1]))
              { 
                indx_X = std::floor((dim1-gridX[0])/dx);
                indx_Z = std::floor((z-gridZ[0])/dz);

                dy = gridY[1] - gridY[0];

                indx_Y = std::floor((y-gridY[0])/dy);

                if (indx_Y < 0 || indx_Y >= nY) indx_Y = 0;

                nnYY = nY;

              if (indx_X < 0 || indx_X >= nX) indx_X = 0;

              if (indx_Z < 0 || indx_Z >= nZ) indx_Z = 0;

              int charge = std::floor(particlesPointer->charge[indx]);

              gitr_precision specWeight = 0.0;

              if(particlesPointer->hitWall[indx]== 0.0)
              {
                if (flags->USE_ADAPTIVE_DT) 
                {
	                if(particlesPointer->advance[indx])
                  {
                    dt_particle = particlesPointer->dt[indx];
                    specWeight = particlesPointer->weight[indx]*dt_particle/dt;
		              }
                }

                else
                {
                  specWeight = particlesPointer->weight[indx];
                }
#if USE_CUDA >0

	              int index = nBins*nX*nnYY*nZ + indx_Z*nX*nnYY +indx_Y*nX+ indx_X;

                atomicAdd1(&bins[index], specWeight);//0*nX*nZ + indx_Z*nZ + indx_X

              if(charge < nBins)
              {
                atomicAdd1( &bins[charge*nX*nnYY*nZ + indx_Z*nX*nnYY + indx_Y*nX+ indx_X],
                            1.0*specWeight);
              }

#else

              #pragma omp atomic
              bins[nBins*nX*nnYY*nZ + indx_Z*nX*nnYY  +indx_Y*nX +indx_X] += specWeight;

              if(charge < nBins)
              {
                #pragma omp atomic
                bins[charge*nX*nnYY*nZ + indx_Z*nX*nnYY +indx_Y*nX + indx_X] += specWeight;
              }
#endif
            }
          }
          }// end
       }
    }
