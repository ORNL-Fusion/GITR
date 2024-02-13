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
          gitr_precision *_particle_time_histogram
          gitr_precision _bin_edge_0_angle,
          gitr_precision _bin_edge_1_angle,
          gitr_precision _bin_edge_dtheta,
          int _n_bins_angle,
          gitr_precision *_particle_angle_histogram)
      
        :

        flags(_flags), particlesPointer(_particlesPointer), boundaryVector(_boundaryVector),  times_logarithmic(_times_logarithmic),
                       bin_edge_0_time(_bin_edge_0_time), bin_edge_1_time(_bin_edge_1_time), bin_edge_dt(_bin_edge_dt),
                       n_bins_time(_n_bins_time), particle_time_histogram(_particle_time_histogram), bin_edge_0_angle(_bin_edge_0_angle), 
                       bin_edge_1_angle(_bin_edge_1_angle), bin_edge_dtheta(_bin_edge_dtheta), n_bins_angle(_n_bins_angle), 
                       particle_angle_histogram(_particle_angle_histogram)
{ }

CUDA_CALLABLE_MEMBER_DEVICE
void particle_diagnostics::operator()(std::size_t indx)
{
  if (particlesPointer->hitWall[indx] == 1.0) 
  {
    int wallHit = particlesPointer->surfaceHit[indx];
    int surfaceHit = boundaryVector[wallHit].surfaceNumber;

    gitr_precision p_time = particlesPointer->time[indx] - particlesPointer->transitTime[indx];
    particlesPointer->transitTime[indx] = particlesPointer->time[indx];

    int ind_time = std::floor((std::log10(p_time) - bin_edge_0_time)/bin_edge_dt);

    int ind_2d = surfaceHit*n_bins_time + ind_time;

    if (ind_time >=0 && ind_time < n_bins_time)
    {
             #if USE_CUDA > 0
               atomicAdd1(&particle_time_histogram[ind_2d],particlesPointer->weight[indx]);
             #else      
               #pragma omp atomic
               particle_time_histogram[ind_2d] = particle_time_histogram[ind_2d] + particlesPointer->weight[indx];
             #endif
    }
}

