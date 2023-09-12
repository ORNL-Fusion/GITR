#include "surfaceModel.h"
#include "processSurfacesModels.h"
#include "interp2d.hpp"
#include "materials.h"

    /* constructor */
    reflection::reflection(Particles* _particles, double _dt,
#if __CUDACC__
                            curandState *_state,
#else
                            std::mt19937 *_state,
#endif
    int _nLines,Boundary * _boundaryVector,
    Surfaces * _surfaces,
    int flux_ea_,
    int use_3d_geom_,
    int cylsymm_,
    int nspecies_) :
                             particles(_particles),
                             dt(_dt),
                             nLines(_nLines),
                             boundaryVector(_boundaryVector),
                             surfaces(_surfaces),
                             state(_state),
                             flux_ea( flux_ea_ ),
                             use_3d_geom( use_3d_geom_ ),
                             cylsymm( cylsymm_ ),
                              nspecies( nspecies_ )
                             { }

CUDA_CALLABLE_MEMBER_DEVICE
void reflection::operator()(std::size_t indx) const {
    
     if (particles->hitWall[indx] == 1.0) 
  {

  int wallHit = particles->surfaceHit[indx];
  int surfaceHit = boundaryVector[wallHit].surfaceNumber;
  int surface = boundaryVector[wallHit].surface;
  std::string targetMaterial = materialData[boundaryVector[wallHit].Z].name;  
  std::string incidentMaterial = materialData[particles->Z[indx]].name;

  printf("material: %s\n",targetMaterial.c_str());
  printf("incidentName: %s\n",incidentMaterial.c_str());
  
  auto surfaceDataResult = processSurfaceData(incidentMaterial, targetMaterial); 

  auto [nE_sputtRefCoeff, nA_sputtRefCoeff, E_sputtRefCoeff, A_sputtRefCoeff, Elog_sputtRefCoeff, energyDistGrid01, angleDistGrid01, spyl_surfaceModel, rfyl_surfaceModel,nE_sputtRefDistIn, nA_sputtRefDistIn,
        E_sputtRefDistIn, A_sputtRefDistIn, Elog_sputtRefDistIn, E_sputtRefDistOut, E_sputtRefDistOutRef, Aphi_sputtRefDistOut, Atheta_sputtRefDistOut, EDist_Y, AphiDist_Y, AthetaDist_Y,
      EDist_R, AphiDist_R, AthetaDist_R, nDistE_surfaceModel, nDistA_surfaceModel, nDistE_surfaceModelRef, nE_sputtRefDistOut, nA_sputtRefDistOut,
      energyDistGrid01Ref, EDist_CDF_Y, AphiDist_CDF_Y, AthetaDist_CDF_Y, EDist_CDF_R, AphiDist_CDF_R, AthetaDist_CDF_R, EDist_CDF_Y_regrid, AphiDist_CDF_Y_regrid,
        AthetaDist_CDF_Y_regrid, EDist_CDF_R_regrid, ADist_CDF_R_regrid, AthetaDist_CDF_R_regrid,
      spylAInterpVal, spylAthetaInterpVal, sputEInterpVal, rfylAInterpVal, rfylAthetaInterpVal, rflEInterpVal,
        nEdist, nAdist, E0dist, Edist, A0dist, Adist, nE_sputtRefDistOutRef
      ] = surfaceDataResult;

    gitr_precision E0_for_surface_model = 0.0;
    gitr_precision E0_for_flux_binning = 0.0;
    gitr_precision thetaImpact = 0.0;
    gitr_precision particleTrackVector[3] = {0.0};
    gitr_precision surfaceNormalVector[3] = {0.0};
    gitr_precision vSampled[3] = {0.0};
    gitr_precision norm_part = 0.0;
    int signPartDotNormal = 0;
    gitr_precision partDotNormal = 0.0;
    gitr_precision Enew = 0.0;
    gitr_precision angleSample = 0.0;
    int wallIndex = 0;
    gitr_precision tol = 1e12;
    gitr_precision Sr = 0.0;
    gitr_precision St = 0.0;
    gitr_precision Y0 = 0.0;
    gitr_precision R0 = 0.0;
    gitr_precision totalYR = 0.0;
    gitr_precision newWeight = 0.0;
    gitr_precision eInterpVal = 0.0;
    gitr_precision aInterpVal = 0.0;
    gitr_precision weight = particles->weight[indx];
    gitr_precision vx = particles->vx[indx];
    gitr_precision vy = particles->vy[indx];
    gitr_precision vz = particles->vz[indx];
    int AdistInd = 0;
    int EdistInd = 0;
    gitr_precision dEdist;
    gitr_precision dAdist;
    
    if( flux_ea > 0 )
    {
      dEdist = (Edist - E0dist) / static_cast<gitr_precision>(nEdist);
      dAdist = (Adist - A0dist) / static_cast<gitr_precision>(nAdist);
    }

    particles->firstCollision[indx] = 1;
    particleTrackVector[0] = vx;
    particleTrackVector[1] = vy;
    particleTrackVector[2] = vz;
    norm_part = std::sqrt(particleTrackVector[0] * particleTrackVector[0] + particleTrackVector[1] * particleTrackVector[1] + particleTrackVector[2] * particleTrackVector[2]);
    E0_for_surface_model = 0.5 * particles->amu[indx] * 1.6737236e-27 * (norm_part * norm_part) / 1.60217662e-19;
    E0_for_flux_binning = E0_for_surface_model;
    gitr_precision maxE_for_surface_model = std::pow(10.0,Elog_sputtRefCoeff[nE_sputtRefCoeff-1]);
    if (E0_for_surface_model > maxE_for_surface_model)
        E0_for_surface_model = maxE_for_surface_model;
     
    if (E0_for_flux_binning > Edist)
        E0_for_flux_binning = Edist;
      
    wallIndex = particles->wallIndex[indx];
    boundaryVector[wallHit].getSurfaceNormal(surfaceNormalVector, particles->y[indx], particles->x[indx], use_3d_geom, cylsymm );
    particleTrackVector[0] = particleTrackVector[0] / norm_part;
    particleTrackVector[1] = particleTrackVector[1] / norm_part;
    particleTrackVector[2] = particleTrackVector[2] / norm_part;

    int species_indx = particles->species[indx];
    partDotNormal = vectorDotProduct(particleTrackVector, surfaceNormalVector);
    thetaImpact = std::acos(partDotNormal);
    if (thetaImpact > 3.14159265359 * 0.5) {
      thetaImpact = std::abs(thetaImpact - (3.14159265359));
    }
    thetaImpact = thetaImpact * 180.0 / 3.14159265359;
    if (thetaImpact < 0.0)
      thetaImpact = 0.0;
    signPartDotNormal = std::copysign(1.0,partDotNormal);
    if (E0_for_surface_model == 0.0) {
      thetaImpact = 0.0;
    }
    if (boundaryVector[wallHit].Z > 0.0) 
    {
      Y0 = interp2d(thetaImpact, std::log10(E0_for_surface_model), nA_sputtRefCoeff, nE_sputtRefCoeff, A_sputtRefCoeff.data(), Elog_sputtRefCoeff.data(), spyl_surfaceModel.data()); 
      R0 = interp2d(thetaImpact, std::log10(E0_for_surface_model), nA_sputtRefCoeff, nE_sputtRefCoeff, A_sputtRefCoeff.data(), Elog_sputtRefCoeff.data(), rfyl_surfaceModel.data());
    } 
    else 
    {
      Y0 = 0.0;
      R0 = 0.0;
    }
    
    totalYR = Y0 + R0;

#ifdef __CUDACC__
    gitr_precision r7 = curand_uniform(&state[indx]);
    gitr_precision r8 = curand_uniform(&state[indx]);
    gitr_precision r9 = curand_uniform(&state[indx]);
    gitr_precision r10 = curand_uniform(&state[indx]);
#else
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    gitr_precision r7 = dist(state[indx]);
    gitr_precision r8 = dist(state[indx]);
    gitr_precision r9 = dist(state[indx]);
    gitr_precision r10 = dist(state[indx]);
#endif

    //particle either reflects or deposits
    gitr_precision sputtProb = Y0/totalYR;
    int didReflect = 0;
            
    if(totalYR > 0.0)
    {
      if(r7 > sputtProb) //reflects
      {
        didReflect = 1;
        aInterpVal = interp3d(r8, thetaImpact, std::log10(E0_for_surface_model),
                      nA_sputtRefDistOut, nA_sputtRefDistIn, nE_sputtRefDistIn,
                      angleDistGrid01.data(), A_sputtRefDistIn.data(),
                      E_sputtRefDistIn.data(), ADist_CDF_R_regrid.data());

        eInterpVal = interp3d(r9, thetaImpact, std::log10(E0_for_surface_model),
                      nE_sputtRefDistOutRef, nA_sputtRefDistIn, nE_sputtRefDistIn,
                      energyDistGrid01Ref.data(), A_sputtRefDistIn.data(),
                      E_sputtRefDistIn.data(), EDist_CDF_R_regrid.data());
      
      
  
         newWeight = weight*(totalYR);

        if( flux_ea > 0 )
        {
        EdistInd = std::floor((eInterpVal-E0dist)/dEdist);
        AdistInd = std::floor((aInterpVal-A0dist)/dAdist);
        if((EdistInd >= 0) && (EdistInd < nEdist) && 
           (AdistInd >= 0) && (AdistInd < nAdist))
        {
        #if USE_CUDA > 0
              atomicAdd1(&surfaces->reflDistribution[surfaceHit*nEdist*nAdist + EdistInd*nAdist + AdistInd],newWeight);
        #else      
              surfaces->reflDistribution[surfaceHit*nEdist*nAdist + EdistInd*nAdist + AdistInd] = 
                surfaces->reflDistribution[surfaceHit*nEdist*nAdist + EdistInd*nAdist + AdistInd] +  newWeight;
        #endif
        }
        }
        if(surface > 0)
        {

        #if USE_CUDA > 0
                atomicAdd1(&surfaces->grossDeposition[surfaceHit + species_indx],weight*(1.0-R0));
                atomicAdd1(&surfaces->grossErosion[surfaceHit + species_indx ],weight*Y0);
        #else
                #pragma omp atomic
                printf("Charge of the particle hitting surface %d is %g\n", surfaceHit, particles->Z[indx]);
                int species_indx = particles->species[indx];
                surfaces->grossDeposition[surfaceHit + species_indx ] += ( weight*(1.0-R0) );
                // surfaces->grossDeposition[surfaceHit] += ( weight*(1.0-R0) );
                #pragma omp atomic
                // surfaces->grossErosion[surfaceHit] += ( weight * Y0 );
                surfaces->grossErosion[surfaceHit + species_indx] += ( weight * Y0 );
        #endif
        }
      }

      else //sputters
      {
        aInterpVal = interp3d(r8,thetaImpact,std::log10(E0_for_surface_model),
                nA_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
                angleDistGrid01.data(),A_sputtRefDistIn.data(),
                E_sputtRefDistIn.data(),AphiDist_CDF_Y_regrid.data());
        eInterpVal = interp3d(r9,thetaImpact,std::log10(E0_for_surface_model),
                 nE_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
                 energyDistGrid01.data(),A_sputtRefDistIn.data(),E_sputtRefDistIn.data(),EDist_CDF_Y_regrid.data());
 
        newWeight=weight*totalYR;

         // Sputtering new particles 
            printf("Sputtering!");
            gitr_precision mass = materialData[boundaryVector[wallHit].Z].mass; 
            gitr_precision Eb = materialData[boundaryVector[wallHit].Z].surfaceBindingEnergy; 
            gitr_precision vTh = std::sqrt(2.0 * Eb * 11600. / (particles->amu[indx] * 1.66e-27));
            gitr_precision vSampled[3];
            vSampled[0] = vTh * std::sin(aInterpVal * 3.1415 / 180) * std::cos(2.0 * 3.1415 * r10);
            vSampled[1] = vTh * std::sin(aInterpVal * 3.1415 / 180) * std::sin(2.0 * 3.1415 * r10);
            vSampled[2] = vTh * std::cos(aInterpVal * 3.1415 / 180);
            // set particle properties
            particles->Z[indx] =  boundaryVector[wallHit].Z;
            particles->amu[indx] =  mass;
            particles->weight[indx] =  newWeight;
            // set species type: nspecies
            particles->species[indx] =  nspecies;

            // Transform velocity based on surface
            boundaryVector[wallHit].transformToSurface(vSampled, particles->y[indx], 
                                                        particles->x[indx], use_3d_geom, 
                                                        cylsymm );
            particles->vx[indx] = -static_cast<gitr_precision>(boundaryVector[wallHit].inDir) * vSampled[0];
            particles->vy[indx] = -static_cast<gitr_precision>(boundaryVector[wallHit].inDir) * vSampled[1];
            particles->vz[indx] = -static_cast<gitr_precision>(boundaryVector[wallHit].inDir) * vSampled[2];

            gitr_precision surface_buffer = 1.0e-4;
            particles->xprevious[indx] = particles->x[indx] - static_cast<gitr_precision>(boundaryVector[wallHit].inDir) * surfaceNormalVector[0] * surface_buffer;
            particles->yprevious[indx] = particles->y[indx] - static_cast<gitr_precision>(boundaryVector[wallHit].inDir) * surfaceNormalVector[1] * surface_buffer;
            particles->zprevious[indx] = particles->z[indx] - static_cast<gitr_precision>(boundaryVector[wallHit].inDir) * surfaceNormalVector[2] * surface_buffer;
            
             printf("particle being sputtered mass Z charge weight %g %g %g %g\n", particles->amu[indx],  particles->Z[indx],  particles->charge[indx],  particles->weight[indx]);

      if( flux_ea > 0 )
      {
        EdistInd = std::floor((eInterpVal-E0dist)/dEdist);
        AdistInd = std::floor((aInterpVal-A0dist)/dAdist);
        if((EdistInd >= 0) && (EdistInd < nEdist) && 
           (AdistInd >= 0) && (AdistInd < nAdist))
        {
        #if USE_CUDA > 0
              atomicAdd1(&surfaces->sputtDistribution[surfaceHit*nEdist*nAdist + EdistInd*nAdist + AdistInd],newWeight);
        #else      
              surfaces->sputtDistribution[surfaceHit*nEdist*nAdist + EdistInd*nAdist + AdistInd] = 
                surfaces->sputtDistribution[surfaceHit*nEdist*nAdist + EdistInd*nAdist + AdistInd] +  newWeight;
        #endif 
        }
      }
        if(sputtProb == 0.0) newWeight = 0.0;
        if(surface > 0)
        {
        #if USE_CUDA > 0
          atomicAdd1(&surfaces->grossDeposition[surfaceHit + species_indx ],weight*(1.0-R0));
          atomicAdd1(&surfaces->grossErosion[surfaceHit],weight*Y0);
          atomicAdd1(&surfaces->aveSputtYld[surfaceHit],Y0);
          if(weight > 0.0)
          {
              atomicAdd1(&surfaces->sputtYldCount[surfaceHit],1.0);
          }
        #else
          #pragma omp atomic
          surfaces->grossDeposition[surfaceHit + species_indx ] = surfaces->grossDeposition[surfaceHit + species_indx ]+weight*(1.0-R0);
          #pragma omp atomic
          surfaces->grossErosion[surfaceHit + species_indx ] = surfaces->grossErosion[surfaceHit + species_indx] + weight*Y0;
          surfaces->aveSputtYld[surfaceHit + species_indx ] = surfaces->aveSputtYld[surfaceHit + species_indx] + Y0;
          surfaces->sputtYldCount[surfaceHit + species_indx] = surfaces->sputtYldCount[surfaceHit + species_indx] + 1;
        #endif
        }
      }
    }
    else
    {       
      newWeight = 0.0;
      particles->hitWall[indx] = 2.0;
      if(surface > 0)
      {
        #if USE_CUDA > 0
                atomicAdd1(&surfaces->grossDeposition[surfaceHit],weight);
        #else
                #pragma omp atomic
                surfaces->grossDeposition[surfaceHit + species_indx] = surfaces->grossDeposition[surfaceHit + species_indx]+weight;
        #endif
	    }
    }
	    
    if(eInterpVal <= 0.0)
    {       
      newWeight = 0.0;
      particles->hitWall[indx] = 2.0;
      if(surface > 0)
      {
      #if USE_CUDA > 0
              atomicAdd1(&surfaces->grossDeposition[surfaceHit + species_indx],weight*R0);
              atomicAdd1(&surfaces->grossDeposition[surfaceHit + species_indx],-weight*Y0);
      #else
              surfaces->grossDeposition[surfaceHit + species_indx ] = surfaces->grossDeposition[surfaceHit + species_indx ]+weight;
      #endif
		  }
    }
    
    if(surface)
    {
    #if USE_CUDA > 0
        atomicAdd1(&surfaces->sumWeightStrike[surfaceHit + species_indx ],weight);
        atomicAdd1(&surfaces->sumParticlesStrike[surfaceHit + species_indx ],1.0);
    #else
        surfaces->sumWeightStrike[surfaceHit + species_indx ] =surfaces->sumWeightStrike[surfaceHit + species_indx ] +weight;
        surfaces->sumParticlesStrike[surfaceHit + species_indx ] = surfaces->sumParticlesStrike[surfaceHit + species_indx ]+1;
    #endif
    if( flux_ea > 0 )
    {
        EdistInd = std::floor((E0_for_flux_binning-E0dist)/dEdist);
        AdistInd = std::floor((thetaImpact-A0dist)/dAdist);
      
	      if((EdistInd >= 0) && (EdistInd < nEdist) && 
        (AdistInd >= 0) && (AdistInd < nAdist))
        {
#if USE_CUDA > 0
          atomicAdd1(&surfaces->energyDistribution[surfaceHit*nEdist*nAdist + 
                                               EdistInd*nAdist + AdistInd], weight);
#else

          surfaces->energyDistribution[surfaceHit*nEdist*nAdist + EdistInd*nAdist + AdistInd] = 
          surfaces->energyDistribution[surfaceHit*nEdist*nAdist + EdistInd*nAdist + AdistInd] +  weight;
#endif
        }
    }
      }
      //reflect with weight and new initial conditions
    if (boundaryVector[wallHit].Z > 0.0 && newWeight > 0.0)
    {
      printf("Reflecting!");
      particles->weight[indx] = newWeight;
      particles->hitWall[indx] = 0.0;
      particles->charge[indx] = 0.0;
      gitr_precision V0 = std::sqrt(2 * eInterpVal * 1.602e-19 / (particles->amu[indx] * 1.66e-27));
      particles->newVelocity[indx] = V0;
      vSampled[0] = V0 * std::sin(aInterpVal * 3.1415 / 180) * std::cos(2.0 * 3.1415 * r10);
      vSampled[1] = V0 * std::sin(aInterpVal * 3.1415 / 180) * std::sin(2.0 * 3.1415 * r10);
      vSampled[2] = V0 * std::cos(aInterpVal * 3.1415 / 180);
      boundaryVector[wallHit].transformToSurface(vSampled, particles->y[indx], particles->x[indx], use_3d_geom,cylsymm );
      particles->vx[indx] = -static_cast<gitr_precision>(boundaryVector[wallHit].inDir)  * vSampled[0];
      particles->vy[indx] = -static_cast<gitr_precision>(boundaryVector[wallHit].inDir)  * vSampled[1];
      particles->vz[indx] = -static_cast<gitr_precision>(boundaryVector[wallHit].inDir)  * vSampled[2];

      gitr_precision surface_buffer = 1.0e-4;
      particles->xprevious[indx] = particles->x[indx] - static_cast<gitr_precision>(boundaryVector[wallHit].inDir) * surfaceNormalVector[0] * surface_buffer;
      particles->yprevious[indx] = particles->y[indx] - static_cast<gitr_precision>(boundaryVector[wallHit].inDir) * surfaceNormalVector[1] * surface_buffer;
      particles->zprevious[indx] = particles->z[indx] - static_cast<gitr_precision>(boundaryVector[wallHit].inDir) * surfaceNormalVector[2] * surface_buffer;

      printf("particle getting reflected mass Z charge weight %g %g %g %g\n", particles->amu[indx],  particles->Z[indx],  particles->charge[indx],  particles->weight[indx]);
      
    } 
    else 
    {
      particles->hitWall[indx] = 2.0;
    } 
  }  
}