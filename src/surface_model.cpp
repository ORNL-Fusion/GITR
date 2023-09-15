#include "surfaceModel.h"
#include "processSurfacesModels.h"
#include "interp2d.hpp"
#include "materials.h"
#include "constants.h"


// initialize constructor
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
{}

CUDA_CALLABLE_MEMBER_DEVICE
void reflection::operator()(std::size_t indx) const {
    if (particles->hitWall[indx] == 1.0) {
        processParticleHit(indx);
    }
}

void reflection::processParticleHit(std::size_t indx) const {
    int wallHit = particles->surfaceHit[indx];
    int surfaceHit = boundaryVector[wallHit].surfaceNumber;
    int surface = boundaryVector[wallHit].surface;
    auto [incidentMaterial, targetMaterial] = getMaterials(wallHit, indx);
    printMaterials(incidentMaterial, targetMaterial);
    auto surfaceDataResult = processSurfaceData(incidentMaterial, targetMaterial); 
    // move this to constructor later
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
    gitr_precision surfaceNormalVector[3] = {0.0};
    gitr_precision Y0 = 0.0;
    gitr_precision R0 = 0.0;
    gitr_precision totalYR = 0.0;
    gitr_precision newWeight = 0.0;
    gitr_precision eInterpVal = 0.0;
    gitr_precision aInterpVal = 0.0;
    gitr_precision weight = particles->weight[indx];

    particles->firstCollision[indx] = 1;
    // get surface data 
    std::pair<gitr_precision, gitr_precision> incidentParticleEnergyAngle = computeIncidentParticleEnergyAngle(particles, indx, use_3d_geom, cylsymm, surfaceNormalVector);
    E0_for_surface_model = incidentParticleEnergyAngle.first;
    thetaImpact = incidentParticleEnergyAngle.second;
    E0_for_flux_binning = E0_for_surface_model;
    gitr_precision maxE_for_surface_model = std::pow(10.0,Elog_sputtRefCoeff[nE_sputtRefCoeff-1]);

    if (E0_for_surface_model > maxE_for_surface_model) E0_for_surface_model = maxE_for_surface_model; 
    if (E0_for_flux_binning > Edist) E0_for_flux_binning = Edist;

    const auto logE0Surface = std::log10(E0_for_surface_model);

    int species_indx = particles->species[indx];

    if (boundaryVector[wallHit].Z > 0.0) 
    {
      Y0 = interp2d(thetaImpact, logE0Surface, nA_sputtRefCoeff, nE_sputtRefCoeff, A_sputtRefCoeff.data(), Elog_sputtRefCoeff.data(), spyl_surfaceModel.data()); 
      R0 = interp2d(thetaImpact, logE0Surface, nA_sputtRefCoeff, nE_sputtRefCoeff, A_sputtRefCoeff.data(), Elog_sputtRefCoeff.data(), rfyl_surfaceModel.data());
    } 
    else 
    {
      Y0 = 0.0;
      R0 = 0.0;
    }
    totalYR = Y0 + R0;
    
    // get random values for device
    gitr_precision randomVals[4];
    for (int i = 0; i < 4; i++) {
        randomVals[i] = getRandomValueForDevice(indx);
    }
    //particle either reflects or sputters
    gitr_precision sputtProb = (totalYR > 0.0) ? Y0 / totalYR : 0.0;
    int didReflect = 0;
    
    if(totalYR > 0.0) {
       if (reflectionEvent( randomVals[0],  sputtProb,  totalYR)){
            didReflect = 1;
            aInterpVal = interp3d(randomVals[0], thetaImpact, logE0Surface,
                          nA_sputtRefDistOut, nA_sputtRefDistIn, nE_sputtRefDistIn,
                          angleDistGrid01.data(), A_sputtRefDistIn.data(),
                          E_sputtRefDistIn.data(), ADist_CDF_R_regrid.data());

            eInterpVal = interp3d(randomVals[2], thetaImpact, logE0Surface,
                          nE_sputtRefDistOutRef, nA_sputtRefDistIn, nE_sputtRefDistIn,
                          energyDistGrid01Ref.data(), A_sputtRefDistIn.data(),
                          E_sputtRefDistIn.data(), EDist_CDF_R_regrid.data());
            newWeight = weight*(totalYR);
            //reflect with weight and new initial conditions
            if (boundaryVector[wallHit].Z > 0.0 && newWeight > 0.0)
            {
              reflect(particles, indx, newWeight, eInterpVal, aInterpVal, boundaryVector, wallHit, use_3d_geom, cylsymm, randomVals[3], surfaceNormalVector);      
            } 
          } else if (sputteringEvent( randomVals[0],  sputtProb,  totalYR)){
            aInterpVal = interp3d(randomVals[0],thetaImpact,logE0Surface,
                    nA_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
                    angleDistGrid01.data(),A_sputtRefDistIn.data(),
                    E_sputtRefDistIn.data(),AphiDist_CDF_Y_regrid.data());
            eInterpVal = interp3d(randomVals[2],thetaImpact,logE0Surface,
                    nE_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
                    energyDistGrid01.data(),A_sputtRefDistIn.data(),E_sputtRefDistIn.data(),EDist_CDF_Y_regrid.data());
            newWeight=weight*totalYR;
            if (boundaryVector[wallHit].Z > 0.0 && newWeight > 0.0) {
                sputter( boundaryVector, wallHit, particles, indx, aInterpVal, randomVals[3], newWeight, nspecies, use_3d_geom, cylsymm, surfaceNormalVector);
                if(sputtProb == 0.0) newWeight = 0.0;
            }
         const bool validBoundary = boundaryVector[wallHit].Z > 0.0;
           if (!validBoundary || newWeight <= 0.0) {
             particles->hitWall[indx] = 2.0;  
            }
          if(eInterpVal <= 0.0)
            {       
              newWeight = 0.0;
              particles->hitWall[indx] = 2.0;
            } 
      }
  }
}

// utility functions
std::pair<std::string, std::string> reflection::getMaterials(int wallHit, std::size_t indx) const {
    std::string targetMaterial = materialData[boundaryVector[wallHit].Z].name;  
    std::string incidentMaterial = materialData[particles->Z[indx]].name;
    return {incidentMaterial, targetMaterial};
}

void reflection::printMaterials(const std::string& incidentMaterial, const std::string& targetMaterial) const {
    printf("material: %s\n", targetMaterial.c_str());
    printf("incidentName: %s\n", incidentMaterial.c_str());
}

void reflection::reflect(Particles* particles, int indx, gitr_precision newWeight, gitr_precision eInterpVal, gitr_precision aInterpVal, 
             Boundary* boundaryVector, int wallHit, bool use_3d_geom, bool cylsymm,
             gitr_precision randomSputterAngle, gitr_precision* surfaceNormalVector) const
{
    printf("Reflecting!");
    particles->weight[indx] = newWeight;
    particles->hitWall[indx] = 0.0;
    particles->charge[indx] = 0.0;

    gitr_precision V0 = std::sqrt(2 * eInterpVal * gitr_constants::e / (particles->amu[indx] *  gitr_constants::m_p));
    particles->newVelocity[indx] = V0;

    gitr_precision vSampled[3];
    vSampled[0] = V0 * std::sin(aInterpVal * M_PI / 180) * std::cos(2.0 * M_PI * randomSputterAngle);
    vSampled[1] = V0 * std::sin(aInterpVal * M_PI / 180) * std::sin(2.0 * M_PI * randomSputterAngle);
    vSampled[2] = V0 * std::cos(aInterpVal * M_PI / 180);

    boundaryVector[wallHit].transformToSurface(vSampled, particles->y[indx], particles->x[indx], use_3d_geom, cylsymm);

    particles->vx[indx] = -static_cast<gitr_precision>(boundaryVector[wallHit].inDir) * vSampled[0];
    particles->vy[indx] = -static_cast<gitr_precision>(boundaryVector[wallHit].inDir) * vSampled[1];
    particles->vz[indx] = -static_cast<gitr_precision>(boundaryVector[wallHit].inDir) * vSampled[2];

    gitr_precision surface_buffer = 1.0e-4;
    particles->xprevious[indx] = particles->x[indx] - static_cast<gitr_precision>(boundaryVector[wallHit].inDir) * surfaceNormalVector[0] * surface_buffer;
    particles->yprevious[indx] = particles->y[indx] - static_cast<gitr_precision>(boundaryVector[wallHit].inDir) * surfaceNormalVector[1] * surface_buffer;
    particles->zprevious[indx] = particles->z[indx] - static_cast<gitr_precision>(boundaryVector[wallHit].inDir) * surfaceNormalVector[2] * surface_buffer;

}


void reflection::sputter(
    Boundary* boundaryVector,
    int wallHit,
    Particles* particles,
    int indx,
    gitr_precision aInterpVal,
    gitr_precision randomSputterAngle,
    gitr_precision newWeight,
    int nspecies,
    bool use_3d_geom,
    bool cylsymm,
    const gitr_precision* surfaceNormalVector
) const
{
    // Sputtering new particles 
    printf("Sputtering!");
    gitr_precision mass = materialData[boundaryVector[wallHit].Z].mass; 
    gitr_precision Eb = materialData[boundaryVector[wallHit].Z].surfaceBindingEnergy; 
    gitr_precision vTh = std::sqrt(2 * aInterpVal * gitr_constants::e / (mass *  gitr_constants::m_p));
    gitr_precision vSampled[3];
    vSampled[0] = vTh * std::sin(aInterpVal * M_PI / 180) * std::cos(2.0 * M_PI * randomSputterAngle);
    vSampled[1] = vTh * std::sin(aInterpVal * M_PI/ 180) * std::sin(2.0 * M_PI * randomSputterAngle);
    vSampled[2] = vTh * std::cos(aInterpVal * M_PI / 180);
    
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
}


std::pair<gitr_precision, gitr_precision> reflection::computeIncidentParticleEnergyAngle(Particles* particles, int indx, int use_3d_geom, int cylsymm, gitr_precision* surfaceNormalVector) const
{
    gitr_precision particleTrackVector[3] = {0.0};
    int wallHit = particles->surfaceHit[indx];
    particleTrackVector[0] = particles->vx[indx];
    particleTrackVector[1] = particles->vy[indx];
    particleTrackVector[2] = particles->vz[indx];
    gitr_precision amu  = particles->amu[indx];

    gitr_precision norm_part = std::sqrt(particleTrackVector[0] * particleTrackVector[0] + 
                                         particleTrackVector[1] * particleTrackVector[1] + 
                                         particleTrackVector[2] * particleTrackVector[2]);
    
    
    gitr_precision E0_for_surface_model = 0.5 * amu * gitr_constants::m_p * (norm_part * norm_part) / gitr_constants::e;
    

    boundaryVector[wallHit].getSurfaceNormal(surfaceNormalVector, particles->y[indx], particles->x[indx], use_3d_geom, cylsymm);
    
    particleTrackVector[0] = particleTrackVector[0] / norm_part;
    particleTrackVector[1] = particleTrackVector[1] / norm_part;
    particleTrackVector[2] = particleTrackVector[2] / norm_part;

    int species_indx = particles->species[indx];
    gitr_precision partDotNormal = particleTrackVector[0] * surfaceNormalVector[0] + 
                               particleTrackVector[1] * surfaceNormalVector[1] + 
                               particleTrackVector[2] * surfaceNormalVector[2];

    gitr_precision thetaImpact = std::acos(partDotNormal);

    if (thetaImpact > M_PI * 0.5) {
      thetaImpact = std::abs(thetaImpact - M_PI);
    }
    thetaImpact = thetaImpact * 180.0 / M_PI;
    if (thetaImpact < 0.0)      thetaImpact = 0.0;
    if (E0_for_surface_model == 0.0) {thetaImpact = 0.0;}


    return std::make_pair(E0_for_surface_model, thetaImpact);
}


gitr_precision reflection::getRandomValueForDevice(int indx) const{
std::uniform_real_distribution<double> dist(0.0, 1.0);
#ifdef __CUDACC__
    return curand_uniform(&state[indx]);
#else
    return dist(state[indx]);
#endif
}

bool reflection::reflectionEvent(gitr_precision randomReflector, gitr_precision sputtProb, gitr_precision totalYR) const {
    return totalYR > 0.0 && randomReflector > sputtProb;
}

bool reflection::sputteringEvent(gitr_precision randomReflector, gitr_precision sputtProb, gitr_precision totalYR) const {
    return totalYR > 0.0 && randomReflector < sputtProb;
}
