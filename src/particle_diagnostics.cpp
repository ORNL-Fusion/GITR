#include "particle_diagnostics.h"
#include "spectroscopy.h"


particle_diagnostics::particle_diagnostics(Flags *_flags,
          Particles *_particlesPointer,
          Boundary *_boundaryVector,
          bool _times_logarithmic,
          gitr_precision _bin_edge_0_time,
          gitr_precision _bin_edge_1_time,
          gitr_precision _bin_edge_dt,
          int _n_bins_time,
          gitr_precision *_particle_time_histogram,
          bool _angle_logarithmic,
          gitr_precision _bin_edge_0_angle,
          gitr_precision _bin_edge_1_angle,
          gitr_precision _bin_edge_dtheta,
          int _n_bins_angle,
          gitr_precision *_particle_angle_histogram,
          int _nSurfaces)
      
        :

        flags(_flags), particlesPointer(_particlesPointer), boundaryVector(_boundaryVector),  times_logarithmic(_times_logarithmic),
                       bin_edge_0_time(_bin_edge_0_time), bin_edge_1_time(_bin_edge_1_time), bin_edge_dt(_bin_edge_dt),
                       n_bins_time(_n_bins_time), particle_time_histogram(_particle_time_histogram), angle_logarithmic(_angle_logarithmic), 
                       bin_edge_0_angle(_bin_edge_0_angle), bin_edge_1_angle(_bin_edge_1_angle), bin_edge_dtheta(_bin_edge_dtheta), 
                       n_bins_angle(_n_bins_angle), particle_angle_histogram(_particle_angle_histogram), nSurfaces(_nSurfaces)
{ }

CUDA_CALLABLE_MEMBER_DEVICE
void particle_diagnostics::operator()(std::size_t indx)
{
  if (particlesPointer->hitWall[indx] == 1.0) 
  {
    std::cout << "DEBUG TEST:\n" << std::endl;
    std::cout << particlesPointer->weight[indx] << std::endl;
            
    int wallHit = particlesPointer->surfaceHit[indx];
    int surfaceHit = boundaryVector[wallHit].surfaceNumber;

    gitr_precision p_time = particlesPointer->time[indx] - particlesPointer->transitTime[indx];
    particlesPointer->transitTime[indx] = particlesPointer->time[indx];

    int ind_time = std::floor((std::log10(p_time) - bin_edge_0_time)/bin_edge_dt);
    if (times_logarithmic != 1)  
    {   
      ind_time = std::floor((p_time - bin_edge_0_time)/bin_edge_dt);
    }

    if (ind_time < 0) ind_time = 0;
    if (ind_time >= n_bins_time) ind_time = n_bins_time - 1;

    int ind_2d_time = surfaceHit*n_bins_time + ind_time;
                        
    if (ind_2d_time >=0 && ind_2d_time < nSurfaces*n_bins_time)
    {
             #if USE_CUDA > 0
               atomicAdd1(&particle_time_histogram[ind_2d_time],particlesPointer->weight[indx]);
             #else      
               #pragma omp atomic
               particle_time_histogram[ind_2d_time] = particle_time_histogram[ind_2d_time] + particlesPointer->weight[indx];
             #endif
    }
  
    gitr_precision p_angle = particlesPointer->angle[indx] - particlesPointer->transitAngle[indx];
    particlesPointer->transitAngle[indx] = particlesPointer->angle[indx];

    int ind_angle = std::floor((p_angle - bin_edge_0_angle)/bin_edge_dtheta);
    if (angle_logarithmic == 1)
    {
      ind_angle = std::floor((std::log10(p_angle) - bin_edge_0_angle)/bin_edge_dtheta);
    }

    if (ind_angle < 0) ind_angle = 0;
    if (ind_angle >= n_bins_angle) ind_angle = n_bins_angle - 1;
            
    int ind_2d_angle = surfaceHit*n_bins_angle + ind_angle;

    if (ind_2d_angle >=0 && ind_2d_angle < nSurfaces*n_bins_angle)
    {
             #if USE_CUDA > 0
               atomicAdd1(&particle_angle_histogram[ind_2d_angle],particlesPointer->weight[indx]);
             #else      
               #pragma omp atomic
               particle_angle_histogram[ind_2d_angle] = particle_angle_histogram[ind_2d_angle] + particlesPointer->weight[indx];
             #endif
    }
  }
}

