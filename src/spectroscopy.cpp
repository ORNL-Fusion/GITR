#include "spectroscopy.h"

spec_bin::spec_bin(class flags &f_init, Particles *_particlesPointer, int _nBins,int _nX,int _nY, int _nZ, gitr_precision *_gridX,gitr_precision *_gridY,gitr_precision *_gridZ,
           gitr_precision* _bins, gitr_precision _dt, int cylsymm_, int spectroscopy_,
           gitr_precision* _bins_vx, gitr_precision* _bins_vy, gitr_precision* _bins_vz,
           gitr_precision* _bins_E ) : 
        f(f_init), particlesPointer(_particlesPointer), nBins(_nBins),nX(_nX),nY(_nY), nZ(_nZ), gridX(_gridX),gridY(_gridY),gridZ(_gridZ), bins(_bins),
        dt(_dt), cylsymm( cylsymm_ ), spectroscopy( spectroscopy_ ),
        bins_vx(_bins_vx), bins_vy(_bins_vy), bins_vz(_bins_vz), bins_E(_bins_E) {}

CUDA_CALLABLE_MEMBER_DEVICE    
void spec_bin::operator()(std::size_t indx) const {

  if(particlesPointer->hitWall[indx]== 0.0)
  {
    gitr_precision dx = 0.0;
    gitr_precision dy = 0.0;
    gitr_precision dz = 0.0;

    gitr_precision x = particlesPointer->xprevious[indx];
    gitr_precision y = particlesPointer->yprevious[indx];
    gitr_precision z = particlesPointer->zprevious[indx];
    
    gitr_precision vx = particlesPointer->vx[indx];
    gitr_precision vy = particlesPointer->vy[indx];
    gitr_precision vz = particlesPointer->vz[indx];

    gitr_precision E0 = 0.5 * particlesPointer->amu[indx] * 1.6737236e-27 * (vx*vx + vy*vy + vz*vz) / 1.60217662e-19;

    int charge = std::floor(particlesPointer->charge[indx]);
    gitr_precision dt_particle = 0.0;
    gitr_precision specWeight = 0.0;

    gitr_precision dim1;
        
    int indx_X=0;
    int indx_Z=0;
    int indx_Y=0;
    int nnYY=1;

    bool add=false;

    // Determine particle dt and relative contribution
    //if (flags->USE_ADAPTIVE_DT)
    if (f.adaptive_dt)
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

    // Determine dimension 1 variable
    if( spectroscopy > 2 )
    {
      dim1 = x;
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
  
    // Determine indices
    if ((z > gridZ[0]) && (z < gridZ[nZ-1]))
    {
      if((dim1 > gridX[0]) && (dim1 < gridX[nX-1]))
      {
        dx = gridX[1] - gridX[0];
        dz = gridZ[1] - gridZ[0];

        indx_X = std::floor((dim1-gridX[0])/dx);
        indx_Z = std::floor((z-gridZ[0])/dz);
        indx_Y = 0;
        nnYY=1;

        if(spectroscopy < 3) add = true;
        
        if (indx_X < 0 || indx_X >= nX) indx_X = 0;

        if (indx_Z < 0 || indx_Z >= nZ) indx_Z = 0;

        if( spectroscopy > 2 )
        {
        
          if((y > gridY[0]) && (y < gridY[nY-1]))
          {
            dy = gridY[1] - gridY[0];
            
            add = true;

            indx_Y = std::floor((y-gridY[0])/dy);

            if (indx_Y < 0 || indx_Y >= nY) indx_Y = 0;

            nnYY = nY;

          }
        }
      }
    }

    gitr_precision vr = vx;
    gitr_precision vt = vy;

    if( spectroscopy == 2 && cylsymm > 0)
    {
          gitr_precision theta_position = std::atan2(y,x);
          gitr_precision theta_velocity = std::atan2(vy,vx);
          gitr_precision v0 = std::sqrt(vy*vy + vx*vx);

          vr = v0*std::cos(theta_position - theta_velocity);
          vt = -v0*std::sin(theta_position - theta_velocity);
    }

    if (add)
    {
    // Total density index
    int index = nBins*nX*nnYY*nZ + indx_Z*nX*nnYY +indx_Y*nX+ indx_X;
    int index_charge = charge*nX*nnYY*nZ + indx_Z*nX*nnYY + indx_Y*nX+ indx_X;
    
    int index_flux = indx_Z*nX*nnYY + indx_Y*nX+ indx_X;

#if USE_CUDA >0
          atomicAdd1(&bins[index], specWeight);


          if(charge < nBins)
          {
            atomicAdd1( &bins[index_charge],specWeight);
          }

          atomicAdd1(&bins_vx[index_flux], specWeight*vr);
          atomicAdd1(&bins_vy[index_flux], specWeight*vt);
          atomicAdd1(&bins_vz[index_flux], specWeight*vz);
          
          atomicAdd1(&bins_E[index_flux], specWeight*E0);
#else

          #pragma omp atomic
          bins[index] += specWeight;

          if(charge < nBins)
          {
            #pragma omp atomic
            bins[index_charge] += specWeight;
          }

          #pragma omp atomic
          bins_vx[index_flux] += specWeight*vr;
          #pragma omp atomic
          bins_vy[index_flux] += specWeight*vt;
          #pragma omp atomic
          bins_vz[index_flux] += specWeight*vz;
          
          #pragma omp atomic
          bins_E[index_flux] += specWeight*E0;
#endif
    }
  }
}
