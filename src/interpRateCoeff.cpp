#include "interpRateCoeff.hpp"

gitr_precision rateCoeffInterp(int charge, gitr_precision te, gitr_precision ne,int nT, int nD, gitr_precision* rateGrid_Tempp,gitr_precision* rateGrid_Densp,gitr_precision* Ratesp){
    
    // indices for temperature and density
    int indT = 0;
    int indN = 0;

    // Grids are log-linear, so we do interpolation in log scale
    gitr_precision logT = std::log10(te);
    gitr_precision logn = std::log10(ne);

    gitr_precision d_T = rateGrid_Tempp[1] - rateGrid_Tempp[0];
    gitr_precision d_n = rateGrid_Densp[1] - rateGrid_Densp[0];
    
    // Find lower index
    indT = std::floor((logT - rateGrid_Tempp[0])/d_T );
    indN = std::floor((logn - rateGrid_Densp[0])/d_n );

    // Safeguards against out of bounds
    if(indT < 0 || indT > nT-2)
    {indT = 0;}
    
    if(indN < 0 || indN > nD-2)
    {indN = 0;}
    
    if(charge > 74-1)
    {charge = 0;}
    
    // Linear interpolation
    gitr_precision aT = std::pow(10.0,rateGrid_Tempp[indT+1]) - te;
    gitr_precision bT = te - std::pow(10.0,rateGrid_Tempp[indT]);
    gitr_precision abT = aT+bT;

    gitr_precision aN = std::pow(10.0,rateGrid_Densp[indN+1]) - ne;
    gitr_precision bN = ne - std::pow(10.0, rateGrid_Densp[indN]);
    gitr_precision abN = aN + bN;

    gitr_precision fx_z1 = (aN*std::pow(10.0,Ratesp[charge*nT*nD + indT*nD + indN]) 
            + bN*std::pow(10.0,Ratesp[charge*nT*nD            + indT*nD + indN + 1]))/abN;
    
    gitr_precision fx_z2 = (aN*std::pow(10.0,Ratesp[charge*nT*nD            + (indT+1)*nD + indN]) 
            + bN*std::pow(10.0,Ratesp[charge*nT*nD            + (indT+1)*nD + indN+1]))/abN;
    gitr_precision fxz = (aT*fx_z1+bT*fx_z2)/abT;
    
    return fxz;
}

gitr_precision interpRateCoeff2d ( int charge, gitr_precision x, gitr_precision y, gitr_precision z,int nx, int nz, gitr_precision* tempGridxp,
       gitr_precision* tempGridzp, gitr_precision* Tempp,
       gitr_precision* densGridxp,gitr_precision* densGridzp,gitr_precision* Densp,int nT_Rates, int nD_Rates,
       gitr_precision* rateGrid_Temp,gitr_precision* rateGrid_Dens,gitr_precision* Rates,
       int cylsymm )
{
    
  gitr_precision tlocal = interp2dCombined( x,y,z,nx,nz,tempGridxp,tempGridzp,Tempp, cylsymm );
  gitr_precision nlocal = interp2dCombined(x,y,z,nx,nz,densGridxp,densGridzp,Densp, cylsymm );
  gitr_precision RClocal = rateCoeffInterp(charge,tlocal,nlocal,nT_Rates,nD_Rates,rateGrid_Temp, rateGrid_Dens, Rates);
  gitr_precision tion = 1.0/(RClocal*nlocal);
  if(tlocal == 0.0 || nlocal == 0.0) tion=1.0e12;
  return tion;
}
