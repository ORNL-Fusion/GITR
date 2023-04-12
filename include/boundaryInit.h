#ifndef _BOUNDARYINIT_
#define _BOUNDARYINIT_


#include "boris.h"
#include "interp2d.hpp"
#include "Particle.h"
#include "Boundary.h"
#ifdef __CUDACC__
#include <thrust/random.h>
#endif

#ifdef __GNUC__ 
#include <random>
#endif
#include <cmath>

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

struct boundary_init {
    gitr_precision background_Z;
    gitr_precision background_amu;
    int nR_Temp;
    int nZ_Temp;
    gitr_precision* TempGridr;
    gitr_precision* TempGridz;
    gitr_precision* ti;
    gitr_precision* te;
    int nx;
    int nz;
    gitr_precision* densityGridx;
    gitr_precision* densityGridz;
    gitr_precision* density;
    gitr_precision* ne;
    int nxB;
    int nzB;
    gitr_precision* bfieldGridr;
    gitr_precision* bfieldGridz;
    gitr_precision* bfieldR;
    gitr_precision* bfieldZ;
    gitr_precision* bfieldT;
    gitr_precision potential;
    int biased_surface;
    int surface_potential;
    int use_3d_geom;
    int cylsymm;
    
    boundary_init(gitr_precision _background_Z, gitr_precision _background_amu,int _nx, int _nz,
          gitr_precision* _densityGridx, gitr_precision* _densityGridz,gitr_precision* _density,gitr_precision* _ne,int _nxB,
          int _nzB, gitr_precision* _bfieldGridr, gitr_precision* _bfieldGridz,gitr_precision* _bfieldR,
          gitr_precision* _bfieldZ,gitr_precision* _bfieldT,int _nR_Temp, int _nZ_Temp,
          gitr_precision* _TempGridr, gitr_precision* _TempGridz, gitr_precision* _ti, gitr_precision* _te, gitr_precision _potential, int biased_surface_, int surface_potential_,
          int use_3d_geom_, int cylsymm_ )

     : background_Z(_background_Z),
        background_amu(_background_amu),
        nR_Temp(_nR_Temp),
        nZ_Temp(_nZ_Temp),
        TempGridr(_TempGridr),
        TempGridz(_TempGridz),
        ti(_ti),
        te(_te),
        nx(_nx),
        nz(_nz),
        densityGridx(_densityGridx),
        densityGridz(_densityGridz),
        density(_density),
        ne(_ne),
        nxB(_nxB),
        nzB(_nzB),
        bfieldGridr(_bfieldGridr),
        bfieldGridz(_bfieldGridz),
        bfieldR(_bfieldR),
        bfieldZ(_bfieldZ),
        bfieldT(_bfieldT),
        potential(_potential),
        biased_surface( biased_surface_ ),
        surface_potential( surface_potential_ ),
        use_3d_geom( use_3d_geom_ ),
        cylsymm( cylsymm_ )
        {}

    void operator()(Boundary &b) const {
        gitr_precision midpointx;
        gitr_precision midpointy;
        gitr_precision midpointz;
    if( use_3d_geom )
    {
        midpointx = b.x1 + 0.666666667*(b.x2 + 0.5*(b.x3-b.x2)-b.x1);
        midpointy = b.y1 + 0.666666667*(b.y2 + 0.5*(b.y3-b.y2)-b.y1);
        midpointz = b.z1 + 0.666666667*(b.z2 + 0.5*(b.z3-b.z2)-b.z1);
    }
    else
    {

        midpointx = 0.5*(b.x2 - b.x1)+ b.x1;
        midpointy = 0.0;
        midpointz = 0.5*(b.z2 - b.z1) + b.z1;
    }
        b.density = interp2dCombined(midpointx,midpointy,midpointz,nx,nz,densityGridx,densityGridz,density, cylsymm );
        b.ne = interp2dCombined(midpointx,midpointy,midpointz,nx,nz,densityGridx,densityGridz,ne, cylsymm );
        b.ti = interp2dCombined(midpointx,midpointy,midpointz,nR_Temp,nZ_Temp,TempGridr,TempGridz,ti, cylsymm );
        b.te = interp2dCombined(midpointx,midpointy,midpointz,nR_Temp,nZ_Temp,TempGridr,TempGridz,te, cylsymm );
        gitr_precision B[3] = {0.0,0.0,0.0};
interp2dVector(&B[0],midpointx,midpointy,midpointz,nxB,nzB,bfieldGridr,
                 bfieldGridz,bfieldR,bfieldZ,bfieldT, cylsymm );
        gitr_precision norm_B = vectorNorm(B);
        gitr_precision theta;
    if( use_3d_geom )
    {
        gitr_precision surfNorm[3] = {0.0,0.0,0.0};
        b.getSurfaceNormal(surfNorm,0.0,0.0, use_3d_geom, cylsymm );
        theta = std::acos(vectorDotProduct(B,surfNorm)/(vectorNorm(B)*vectorNorm(surfNorm)));
        if (theta > 3.14159265359*0.5)
        {
          theta = std::abs(theta - (3.14159265359));
        }
        b.unit_vec0 = b.inDir*b.a/b.plane_norm; //
        b.unit_vec1 = b.inDir*b.b/b.plane_norm; //
        b.unit_vec2 = b.inDir*b.c/b.plane_norm;
    }
    else
    {
        gitr_precision br = B[0];
        gitr_precision bt = B[1];
        gitr_precision bz = B[2];
        theta = std::acos((-br*b.slope_dzdx + bz)/(std::sqrt(br*br+bz*bz+bt*bt)*std::sqrt(b.slope_dzdx*b.slope_dzdx + 1.0)));
 
        if (theta > 3.14159265359*0.5)
        {
            theta = std::acos((br*b.slope_dzdx - bz)/(std::sqrt(br*br+bz*bz+bt*bt)*std::sqrt(b.slope_dzdx*b.slope_dzdx + 1.0)));
        }
    }
        b.angle = theta*180.0/3.14159265359;
        b.debyeLength = std::sqrt(8.854187e-12*b.te/(b.ne*std::pow(background_Z,2)*1.60217662e-19));
	//std::cout << "debyeLength " << b.debyeLength << std::endl;
	gitr_precision angle = b.angle;
  gitr_precision me = 1.0/2000.0;
  std::cout << "me " << me << " " << background_amu << " " << b.ti << " " << b.te << std::endl; 
	gitr_precision sheath_fac = std::abs(0.5*std::log((2*M_PI*me/background_amu)*(1+b.ti/b.te)));
	gitr_precision norm = std::acos(std::pow(std::exp(1),-sheath_fac));
	gitr_precision fd = 1.0+std::log(std::cos(angle/90.0*norm))/sheath_fac;
	if(fd < 0.0) fd = 0.0;
	if(b.te <= 0.0) fd = 0.0;
        if(b.te <= 0.0) sheath_fac = 0.0;
	b.fd = fd;
  std::cout << "fd " << fd << " " << sheath_fac << std::endl; 
	if(b.ne == 0.0) b.debyeLength = 1.0e12;
        b.larmorRadius = 1.44e-4*std::sqrt(background_amu*b.ti/2)/(background_Z*norm_B);
        b.flux = 0.25*b.density*std::sqrt(8.0*b.ti*1.602e-19/(3.1415*background_amu));
        b.impacts = 0.0;
        if( biased_surface )
        {
        b.potential = potential;
        //gitr_precision cs = std::sqrt(2*b.ti*1.602e-19/(1.66e-27*background_amu));
        //gitr_precision jsat_ion = 1.602e-19*b.density*cs;
        //b.ChildLangmuirDist = 2.0/3.0*std::pow(2*1.602e-19/(background_amu*1.66e-27),0.25)
        //*std::pow(potential,0.75)/(2.0*std::sqrt(3.1415*jsat_ion))*1.055e-5;
        if(b.te > 0.0)
        {
          b.ChildLangmuirDist = b.debyeLength*std::pow(std::abs(b.potential)/b.te,0.75);
        }
        else
        { b.ChildLangmuirDist = 1e12;
        }
        }
        else if( surface_potential <= 0 )
        {
        b.potential = sheath_fac*b.te;
        std::cout << "Surface number " << b.surfaceNumber << " has te and potential " 
                  << b.te << " " << b.potential << std::endl; 
        }
        //if(b.Z > 0.0)
        //{
        //std::cout << "Boundary ti density potensial and CLdist " <<b.ti << " " << 
        //    b.density << " " << b.potential << " " << b.ChildLangmuirDist << std::endl;   
        //}     
    }	
};

#endif
