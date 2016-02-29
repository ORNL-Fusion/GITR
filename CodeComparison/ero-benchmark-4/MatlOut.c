#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <stdarg.h>

#include "prapro.h"

#include "MatlOut.h"

#ifdef METAST
#include "Metastable.h"
#endif

#include "mathhlp.h"
#include "lim_go.h"
#include "distrib.h"
#include "nT.h"
#include "surface_simple.h"
#include "neutrals.h"
#include "emission.h"
#ifdef BE_CARBYDE
#include "SDTrimSP/trim_interface.h"
#endif
#ifdef JETLim
#include "sputtern.h"
#endif
#ifdef ITER_BM
#include "ITER_BM.h"
#endif
#ifndef UGEOMETRY_OLD
#include "geometry/snetCache.h"
#endif
#ifdef UGEOMETRY
#include "geometry/surface.h"
#endif
#include "local_cs.h"

/* some common definitions for SurfaceCells*/
#define SC_SPUTTERED 2
#define SC_TOTAL_RESOURCE 9
#define SC_DELETED_COMP 11
/* common definitions for CellProps*/
#ifdef TRIM
#define CP_DMDLX       1
#define CP_DMDLY       2
#define CP_ICLX        3
#define CP_ICLY        4
#define CP_DSURFT      5
#define CP_DAREA       6
#define CP_DCELLSIZE   7
#define CP_DSHPOT      8
#define CP_DHFLUC      9 /*WARNING: not written in SDTrimSP case - therefore same number as first SDTrimSP CP_...!*/
#endif
#ifdef BE_CARBYDE 
#define CP_BE2C 10
#endif
/*SDTrimSP definitions*/
#define CP_DDELTADEPTH 9
#define CP_DTOTDELTADEPTH 10
#define CP_QUX 5 /*-> output into three different files- therefore same number as before!*/
#define CP_XNO 6
#define CP_DNS 7
#define CP_PDEPTHBG 8
#define CP_PDEPTHFO 9
#define CP_BE2C_TRIM 10
#define CP_MAX 4 /* CP_MAX properties are written to all three files (redundant data)*/


#ifdef JETLim
#define JETSpotVAR\
  double dTanA = tan(58.4823/180*M_PI); \
  double dCosA = cos(58.4823/180*M_PI); \
  const double dX0 = -85.,  dZ0=51.5;   \
  const double dRSp2 = 3600.;           \
  double dX0_;                          \
  double dX, dY, dZ;                    \
  FILE *fJET;
#endif  /* JETLim */


/* -----------------------------------------------------
 Opens file as 'wt', if neccesary gives error message 
 Saves in it time, computer and task name.
--------------------------------------------------------*/
FILE *MFOpen (const char *FN,const char *sVerName, int iTStep)
{
	FILE *f;
	char sFN[PATH_STR_LEN];

#if defined(JUMP)
    char sVer[]="Version for JUMP";
	char sComp[]="JUMP";
#elif defined(NAGC)
	char sVer[]="Version using NagC library";
	char sComp[]="AIX NagC";
#elif defined(VCpp)
    char sVer[]="Version for MS Visual C++";
	char sComp[]="VC++";
#else
	char sVer[]="Standard Version";
	char sComp[]="STD";
#endif

	time_t aclock;
	struct tm *timeptr; 
    char *sTime;    


	/* char *asctime( const struct tm *timeptr ); */
	time(&aclock);
	timeptr = localtime( &aclock );
	sTime =asctime(timeptr);

	if (iTStep<0)
	  sprintf(sFN, "%s%s.m", sVerName, FN );
	else 
      sprintf(sFN, "%s%s_St%i.m", sVerName, FN, iTStep);

	printf ("Creating ASCII file '%s' . . . ",sFN);
	if ((f=fopen(sFN, "wt"))!=NULL){
		printf("Complete!\n");
		fprintf (f,"%%============================================================\n");
        fprintf (f,"%% Created by ERO (%s)\n%% for the task '%s'.\n%% %s",
			sVer, sVerName, sTime);
		fprintf (f,"%%============================================================\n");
		fprintf (f,"%%\n");
		  fprintf (f, "TskName= '%s';\n",sVerName);
          fprintf (f, "CompName= '%s';\n%%\n",sComp);

#ifdef PISCES 
		fprintf (f, "SpVer= 'PISCES';\n");
#ifdef PISCES1
		fprintf (f, "SpVerSub = 'PSC1.0';\n%%\n");
#else
        fprintf (f, "SpVerSub = 'PSC3.0';\n%%\n");
#endif
#endif

#ifdef UGEOMETRY
	   fprintf (f, "iUGEOMETRY= 1;\n");
#endif

#ifdef ITER_BM 
		fprintf (f, "SpVer= 'ITER_BM';\n");
#ifdef JETLim
		fprintf (f, "SpVerSub= 'JETLim';\n");
#endif 
#endif

#ifdef DIVERTOR 
		fprintf (f, "SpVer= 'DIVERTOR';\n%%\n");
#endif

#ifdef EMC3_BKG /* Markus Airila 310810 */
		fprintf (f, "SpVer= 'EMC3_BKG';\n%%\n");
#endif
		return f;
    }
	else 
	  ERROR_FATAL("Unable to open output file.",sFN); 
	free(timeptr);
	return NULL;
}

#ifdef SURF_RAD_STAT
void RadStatInit(struct StrSURF_RS *StrRS, struct SurfaceCell  **s_net)
{
/*#define SRD0_R       0
#define SRD1_PhysEr  1
#define SRD_MAX      2*/
/*
struct StrSURF_RS{
  int iNRp;       Number of radial points 
  int iNCell;    /* Number of surface Cells 
  int *iCI;      /* Cells radial indexes . . . 
  double *RP[SRD_MAX];
};*/

  /* int iCI[10]; */
}
#endif

/* -------------------------------------------------------------
   Plasma parameter 2D output . . .
   PP saved in the form of "spectroscopy" m-file.
   density - is Te, ne maximum values,
   spec    - is Te, ne averaged values
   ------------------------------------------------------------- */
void MatlPlPar2D(struct VolumeAll  *v_all,			
#ifdef HGAMMA
        struct nTdata               *nT,     
		struct Numeric             *num,        
	    void *pHGDt,
#endif
		int iFMC,
		struct Path	  *pfade,
		struct TimeCell   *time_c,
		int                step,
		struct strGenOptions *GenOpt)
{
  double dx          = v_all->diff[0];
  double dy          = v_all->diff[1];
  double dz          = v_all->diff[2];
  double dMinX       = v_all->adRange[0];
  double dMinY       = v_all->adRange[1];    

  int iCPt, iNPt       = v_all->anz[0];    /* points in toroidal dir. */
  int iCPp, iNPp       = v_all->anz[1];    /* points in poloidal dir. */
  int iCPz, iNPz       = v_all->anz[2];    /* points in vertical dir. */
  int iMVt, iMVp, iMVz;                    /* BCmt: number of cell for cross-section */
#ifdef HGAMMA
  int iCPar, iNPar     = 3;                /* Number of pars: n_e, T_e, H-Gamma, ... */
#else
  int iCPar, iNPar     = 2;                /* Number of pars: n_e, T_e, ... */
#endif

  char sFN_NTrv[DEF_STR_LEN], sBrkts[DEF_STR_LEN];
  char sElN[] ="Iz___Nm(Ch)";
  char sElN_[10];
  char *merged;
  FILE *fMOut;

  double *Buf; 

#ifdef HGAMMA
  double dRedF = 1.0;
#ifdef PISCES
  double dRedF_ = num->dRedF_;
  double dNDens = nT->dDensNtrl / 20;
  /* BCmt: according to 
           D.Whyte et.al., J. Vac. Sci. Technol. A, Vol. 17, No. 5, Sep/Oct 1999
		   N(D) =  N(D2)/20 */
#else
  double dNDens = nT->dDensNtrl;
  double    dRedF_ = 0.;
#endif
  double *BufDg;
  double dTeCur, dNeCur; /* BCmt: current pl.par. */
#endif
  double *adBxy, *adBxz, *adByz; 
  double *adBMxy, *adBMxz, *adBMyz;
  int iZCor, iCZ;  /* Corrected Z coordinates */
  double dCVal;
  char sCP[10];  /* Current prameter */
  /* ------------ vars . . . -------------*/

  if(dMinX<0)
    iMVt = (int)floor(-dMinX/dx);
  else iMVt = 0;
  if(dMinY<0)
    iMVp = (int)floor(-dMinY/dy);
  else iMVp = 0;
#ifdef PISCES
  iMVz = (int)floor(150/dz);
#else
  iMVz = (int)floor(iNPz/2.0);
#endif

  Buf=(double*)calloc(iNPt*iMax2(iNPp,iNPz),sizeof(double));
  if(Buf == NULL) 
	ERROR_FATAL_I ("Problems during memory allocation", iNPt*iMax2(iNPp,iNPz)*sizeof(double) );

/* ---------------  BCmt: H-Gamma preparation of spec. buffer ------------- */
#ifdef HGAMMA
  if (pHGDt == NULL)
	  ERROR_FATAL ("H-Gamma data is missing . . . ", "MatlPlPar2D");

  BufDg = (double*)calloc(iNPt*iNPp*iNPz,sizeof(double));
  if(BufDg == NULL) 
	ERROR_FATAL_I ("Problems during memory allocation", iNPt*iNPp*iNPz*sizeof(double));

   for (iCPt=0; iCPt<iNPt; iCPt++)
	for (iCPp=0; iCPp<iNPp; iCPp++)
      for (iCPz=0; iCPz<iNPz; iCPz++){
		dNeCur = v_all->b[iCPt][iCPp].vcell[iCPz].dNe; 
		dTeCur = v_all->b[iCPt][iCPp].vcell[iCPz].dTe; 

#ifdef PISCES
  if (dRedF_ > 0){
     dRedF = 1.0 - dNeCur/nT->ne * dRedF_;
	 if (dRedF < 0) dRedF = 0.0;
  }else{
     dRedF = 1;
  }
#else
     dRedF = 1;
#endif
		BufDg[iCPt + iCPp*iNPt + iCPz*iNPt*iNPp]=
		   dRedF * dNDens * dNeCur * dNewEmission(pHGDt, dTeCur, dNeCur);
	  }
#endif

  /* ---- Opening file, writing trivial things . . . ---- */
  if (iFMC==0)
    sprintf (sFN_NTrv,"PlPar2D_step%i",step);
  else
    sprintf (sFN_NTrv,"PlPar2D_fm%i_step%i",iFMC,step);
  fMOut= MFOpen (sFN_NTrv, pfade->RESULT, -1);
   
  fprintf (fMOut, "iFMC = %i;\n",iFMC);
  if (iFMC==0)
    fprintf (fMOut,"iflFM=0;  %% Normal mesh\n%%\n");
  else
    fprintf (fMOut,"iflFM=1;  %% Fine mesh\n%%\n");

  fprintf (fMOut, "FileType = 'PlPar2D';\n%%\n");

#ifdef HGAMMA
  fprintf (fMOut, "NtrlDens = %.3e;\n", dNDens);
  fprintf (fMOut, "dRedF    = %.3f;\n%%\n", dRedF_);
#endif

  fprintf (fMOut, "%% %i parameters as elements\n",iNPar);
  fprintf (fMOut, "iNEl = %i;\n%%\n",iNPar);
	
  fprintf (fMOut, "dx= %.2e;\n",dx);
  fprintf (fMOut, "dy= %.2e;\n",dy);
  fprintf (fMOut, "dz= %.2e;\n",dz);
  fprintf (fMOut, "%%\n");

  fprintf (fMOut, "iNPx= %i;\n",iNPt);
  fprintf (fMOut, "iNPy= %i;\n",iNPp);
  fprintf (fMOut, "iNPz= %i;\n",iNPz);
  fprintf (fMOut, "%%\n");

  fprintf (fMOut, "%% Cross-section at:\n");
  fprintf (fMOut, "iMVx= %i;  %%  %.2e mm\n", iMVt, dMinX+iMVt*dx);
  fprintf (fMOut, "iMVy= %i;  %%  %.2e mm\n", iMVp, dMinY+iMVp*dy);
  fprintf (fMOut, "iMVz= %i;  %%  %.2e mm\n", iMVz, iMVz*dz);
  fprintf (fMOut, "%%\n");

  fprintf(fMOut, "iTimeStepNo=%d;\n",step);
  fprintf(fMOut, "dStepTime=%f;\n",time_c->step); 
  fprintf(fMOut, "dTotalTime=%f;\n",time_c->total);
  fprintf(fMOut, "dRealTime=%f;\n",time_c->dRealTm);
  fprintf (fMOut, "%%\n");

  /* ---- Writing arrays of X,Y,Z . . . ---- */
  for (iCPt=0; iCPt<iNPt; iCPt++) Buf[iCPt]=dMinX+iCPt*dx;
  MlMatrOut (fMOut, "dPx", Buf, iNPt, 1, 'd');
	
  for (iCPp=0; iCPp<iNPp; iCPp++) Buf[iCPp]=dMinY+iCPp*dy;
  MlMatrOut (fMOut, "dPy", Buf, iNPp, 1, 'd');

  for (iCPz=0; iCPz<iNPz; iCPz++) Buf[iCPz]=iCPz*dz;
  MlMatrOut (fMOut, "dPz", Buf, iNPz, 1, 'd');

  for (iCPp=0; iCPp<iNPp; iCPp++)
  for (iCPt=0; iCPt<iNPt; iCPt++)
	Buf[iCPt+iCPp*iNPt] = v_all->b[iCPt][iCPp].z_min;
  MlMatrOut (fMOut, "dLimZxy", Buf, iNPt, iNPp, 'd');

  /* ---------- Maximum values -----------*/
  for(iCPar=1; iCPar <= iNPar; iCPar++){
    switch (iCPar){
		case 1: strcpy(sCP,"n_e"); break;
		case 2: strcpy(sCP,"T_e"); break;
#ifdef HGAMMA
		case 3: strcpy(sCP,"HGam"); break;
#endif
		default: ERROR_FATAL_I("Illegal plasma prameter index", iCPar);
	}
    strcpy (sElN,"___________");
    sprintf(sElN_, "PP_%s", sCP);
	strncpy (sElN, sElN_, strlen (sElN_));
	fprintf (fMOut, "%%\n%% -------------------- Parameter #%i '%s' -------------------\n", iCPar, sCP);
    fprintf (fMOut, "ElName(%i,:) = '%s';\n%%\n",iCPar,sElN);

	sprintf (sBrkts,"(%i,:,:)",iCPar);

    /* --------- 2D averaging ------------*/
	adBxz=(double*)calloc(iNPt*iNPz,sizeof(double));
    adByz=(double*)calloc(iNPp*iNPz,sizeof(double));
	adBxy=(double*)calloc(iNPt*iNPp,sizeof(double));
	adBMxz=(double*)calloc(iNPt*iNPz,sizeof(double));
    adBMyz=(double*)calloc(iNPp*iNPz,sizeof(double));
	adBMxy=(double*)calloc(iNPt*iNPp,sizeof(double));
    for (iCPt=0; iCPt<iNPt; iCPt++)
	for (iCPp=0; iCPp<iNPp; iCPp++) { 
        iZCor = (int) ceil (v_all->b[iCPt][iCPp].z_min/dz);
        for (iCPz=0; iCPz<iNPz; iCPz++){

		  switch(iCPar){
		  case 1: dCVal = v_all->b[iCPt][iCPp].vcell[iCPz].dNe; break;
		  case 2: dCVal = v_all->b[iCPt][iCPp].vcell[iCPz].dTe; break;
#ifdef HGAMMA
		  case 3: dCVal = BufDg[iCPt + iCPp*iNPt + iCPz*iNPt*iNPp]; break;
#endif
		  }

          iCZ = iCPz + iZCor;
		  if (iCPar < 3) /* BCmt: averaging for Ne, Te . . . */
	        adBxy[iCPt+iNPt*iCPp]+= dCVal/iNPz;
		  else
            adBxy[iCPt+iNPt*iCPp]+= dCVal;
          /* BCmt: maximum (absolete . . .)
		  if (dCVal>adBMxy[iCPt+iNPt*iCPp])  adBMxy[iCPt+iNPt*iCPp] = dCVal; */
		  if (iMVz == iCPz)  adBMxy[iCPt+iNPt*iCPp] = dCVal;

		    if (iCZ<iNPz){
		  if (iCPar < 3){ /* BCmt: averaging for Ne, Te . . . */
		    adBxz[iCPt+iNPt*iCZ] += dCVal/iNPp;
            adByz[iCPp+iNPp*iCZ] += dCVal/iNPt;
		  }else{
		    adBxz[iCPt+iNPt*iCZ] += dCVal;
            adByz[iCPp+iNPp*iCZ] += dCVal;
		  }
          /* BCmt: maximum (absolete . . .)
		  if (dCVal>adBMxz[iCPt+iNPt*iCZ])  adBMxz[iCPt+iNPt*iCZ] = dCVal;
		  if (dCVal>adBMyz[iCPp+iNPp*iCZ])  adBMyz[iCPp+iNPp*iCZ] = dCVal; */
		  if (iMVp == iCPp) adBMxz[iCPt+iNPt*iCZ] = dCVal;
		  if (iMVt == iCPt) adBMyz[iCPp+iNPp*iCZ] = dCVal; 
			}
        }
	}

	/* --------- 2D averaging ------------*/
    merged=sMerge("Dens2Dxy",sBrkts);
	MlMatrOut (fMOut, merged, adBxy, iNPt, iNPp, 'd');
	free(merged);
  merged=sMerge("Dens2Dxz",sBrkts);
	MlMatrOut (fMOut, merged, adBxz, iNPt, iNPz, 'd');
	free(merged);
  merged=sMerge("Dens2Dyz",sBrkts);
	MlMatrOut (fMOut, merged, adByz, iNPp, iNPz, 'd');
	free(merged);
    free(adBxz);
	free(adBxy);
    free(adByz);
    merged=sMerge("Spec2Dxy",sBrkts);
	MlMatrOut (fMOut, merged, adBMxy, iNPt, iNPp, 'd');
	free(merged);
    merged=sMerge("Spec2Dxz",sBrkts);
	MlMatrOut (fMOut, sMerge("Spec2Dxz",sBrkts), adBMxz, iNPt, iNPz, 'd');
	free(merged);
    merged=sMerge("Spec2Dyz",sBrkts);
	MlMatrOut (fMOut, merged, adBMyz, iNPp, iNPz, 'd');
	free(merged);
    free(adBMxz);
	free(adBMxy);
    free(adBMyz);
	/* --------- 2D averaging (end) --------*/

  }/* iCPar for */
  
#ifdef HGAMMA
  free(BufDg);
#endif
  free(Buf);
  fclose(fMOut);
}

#ifdef TR_ION_ACT
/* -------------------------------------------------------------
    aAdditional file for volume inforamtion output
	intergrated 2D planes for ionization
    averaged 2D planes for plasma parameters . . .
   ------------------------------------------------------------- */
void VolAddOut(struct VolumeAll  *v_all,
		int iFMC,
		struct Path	  *pfade,
		struct TimeCell   *time_c,
		int                step,
		struct strGenOptions *GenOpt)
{
  double dx          = v_all->diff[0];
  double dy          = v_all->diff[1];
  double dz          = v_all->diff[2];
  double dMinX       = v_all->adRange[0];
  double dMinY       = v_all->adRange[1];    

  int iCPt, iNPt       = v_all->anz[0];    /* points in toroidal dir. */
  int iCPp, iNPp       = v_all->anz[1];    /* points in poloidal dir. */
  int iCPz, iNPz       = v_all->anz[2];    /* points in vertical dir. */
  int iCEl, iNEl       = iNElHES();        /* number of elements (H.E.S.)*/
  int iCQ,   iNQ       = TR_ION_ACT;        /* number of elements (H.E.S.)*/

#ifdef PilotPSI
  int iCPr, iNPr; 						   /* for printing output in r/z coordinates */
  int iCTheta, iNTheta;                    /* for printing output in r/z coordinates */
#endif

  char sFN_NTrv[DEF_STR_LEN], sBrkts[DEF_STR_LEN];
  char sElN[] ="Iz___Nm(Ch)";
  char sElN_[10];

  FILE *fMOut;

  double *Buf;
  double *adBxy, *adBxz, *adByz; 
#ifdef PilotPSI
  double *adBrz;
#endif
  int *aiBxy, *aiBxz, *aiByz;
  int iZCor, iCZ;  /* Corrected Z coordinates */
  int iCVal;
  double dCVal;
  int iCElQ;
  char *sCEl;  /* Current element */
  char *merged;
#ifdef TR_DR_PHOT
	int iCHx;
#endif
  /* ------------ vars . . . -------------*/


  Buf=(double*)calloc(iNPt*iMax2(iNPp,iNPz),sizeof(double));
  
  /* ---- Opening file, writing trivial things . . . ---- */
  if (iFMC==0)
    sprintf (sFN_NTrv,"VolumeAdd_step%i",step);
  else
    sprintf (sFN_NTrv,"VolumeAdd_fm%i_step%i",iFMC,step);
  fMOut= MFOpen (sFN_NTrv, pfade->RESULT, -1);
   
  fprintf (fMOut, "iFMC = %i;\n",iFMC);
  if (iFMC==0)
    fprintf (fMOut,"iflFM=0;  %% Normal mesh\n%%\n");
  else
    fprintf (fMOut,"iflFM=1;  %% Fine mesh\n%%\n");

  fprintf (fMOut, "FileType = 'VolAdd';\n%%\n");

  fprintf (fMOut, "%% %i lements, %i charges\n",iNEl, iNQ);
  fprintf (fMOut, "iNEl = %i;\n%%\n",iNEl*iNQ);
	
  fprintf (fMOut, "dx= %.2e;\n",dx);
  fprintf (fMOut, "dy= %.2e;\n",dy);
  fprintf (fMOut, "dz= %.2e;\n",dz);
  fprintf (fMOut, "%%\n");

  fprintf (fMOut, "iNPx= %i;\n",iNPt);
  fprintf (fMOut, "iNPy= %i;\n",iNPp);
  fprintf (fMOut, "iNPz= %i;\n",iNPz);
  fprintf (fMOut, "%%\n");

  fprintf(fMOut, "iTimeStepNo=%d;\n",step);
  fprintf(fMOut, "dStepTime=%f;\n",time_c->step); 
  fprintf(fMOut, "dTotalTime=%f;\n",time_c->total);
  fprintf(fMOut, "dRealTime=%f;\n",time_c->dRealTm);
  fprintf (fMOut, "%%\n");

  /* ---- Writing arrays of X,Y,Z . . . ---- */
  for (iCPt=0; iCPt<iNPt; iCPt++) Buf[iCPt]=dMinX+iCPt*dx;
  MlMatrOut (fMOut, "dPx", Buf, iNPt, 1, 'd');
	
  for (iCPp=0; iCPp<iNPp; iCPp++) Buf[iCPp]=dMinY+iCPp*dy;
  MlMatrOut (fMOut, "dPy", Buf, iNPp, 1, 'd');

  for (iCPz=0; iCPz<iNPz; iCPz++) Buf[iCPz]=iCPz*dz;
  MlMatrOut (fMOut, "dPz", Buf, iNPz, 1, 'd');

  for (iCPp=0; iCPp<iNPp; iCPp++)
  for (iCPt=0; iCPt<iNPt; iCPt++)
	Buf[iCPt+iCPp*iNPt] = v_all->b[iCPt][iCPp].z_min;
  MlMatrOut (fMOut, "dLimZxy", Buf, iNPt, iNPp, 'd');

  /* ---------- Ionization acts . . . -----------*/
  for(iCEl=0; iCEl < iNEl; iCEl++)
  for(iCQ =0; iCQ < iNQ; iCQ++){ 
	sCEl = sHSeqNm(iCEl);
	iCElQ = iCEl*iNQ + iCQ +1; /* currect "Elem + Charge" specie  */
    strcpy (sElN,"___________");
    sprintf(sElN_, "Iz_%s(+%i)", sCEl, iCQ);
	strncpy (sElN, sElN_, strlen (sElN_));
	fprintf (fMOut, "%%\n%% -------------------- Element #%i '%s' -------------------\n",iCElQ,sElN);
    fprintf (fMOut, "ElName(%i,:) = '%s';\n%%\n",iCElQ,sElN);

	sprintf (sBrkts,"(%i,:,:)",iCElQ);

    /* --------- 2D averaging ------------*/
	aiBxz=(int*)calloc(iNPt*iNPz,sizeof(int));
    aiByz=(int*)calloc(iNPp*iNPz,sizeof(int));
	aiBxy=(int*)calloc(iNPt*iNPp,sizeof(int));
    for (iCPt=0; iCPt<iNPt; iCPt++)
	for (iCPp=0; iCPp<iNPp; iCPp++) { 
        iZCor = (int) ceil (v_all->b[iCPt][iCPp].z_min/dz);
        for (iCPz=0; iCPz<iNPz; iCPz++){

		  iCVal = v_all->b[iCPt][iCPp].vcell[iCPz]
		          .iNActIon[iCQ][iCEl];

          iCZ = iCPz + iZCor;
	      aiBxy[iCPt+iNPt*iCPp]+= iCVal;
		    if (iCZ<iNPz){
		  aiBxz[iCPt+iNPt*iCZ] += iCVal;
          aiByz[iCPp+iNPp*iCZ] += iCVal;
			}
        }
	}

	/* --------- 2D averaging ------------*/
    merged = sMerge("Dens2Dxy",sBrkts);
	MlMatrOut (fMOut, merged, aiBxy, iNPt, iNPp, 'i');
	free(merged);
    merged = sMerge("Dens2Dxz",sBrkts);
	MlMatrOut (fMOut, merged, aiBxz, iNPt, iNPz, 'i');
	free(merged);
    merged = sMerge("Dens2Dyz",sBrkts);
	MlMatrOut (fMOut, merged, aiByz, iNPp, iNPz, 'i');
	free(merged);
    free(aiBxz);
	free(aiBxy);
    free(aiByz);
    /* --------- 2D averaging (end) --------*/

	adBxz=(double*)calloc(iNPt*iNPz,sizeof(double));
    adByz=(double*)calloc(iNPp*iNPz,sizeof(double));
	adBxy=(double*)calloc(iNPt*iNPp,sizeof(double));
    for (iCPt=0; iCPt<iNPt; iCPt++)
	for (iCPp=0; iCPp<iNPp; iCPp++) { 
        iZCor = (int) ceil (v_all->b[iCPt][iCPp].z_min/dz);
        for (iCPz=0; iCPz<iNPz; iCPz++){

	  dCVal = v_all->b[iCPt][iCPp].vcell[iCPz]
	  		          .dIzAtoms[iCQ][iCEl];

          iCZ = iCPz + iZCor;
	      adBxy[iCPt+iNPt*iCPp]+= dCVal;
		    if (iCZ<iNPz){
		  adBxz[iCPt+iNPt*iCZ] += dCVal;
          adByz[iCPp+iNPp*iCZ] += dCVal;
			}
        }
	}
    merged=sMerge("Spec2Dxy",sBrkts);
	MlMatrOut (fMOut, merged, adBxy, iNPt, iNPp, 'd');
	free(merged);
	merged=sMerge("Spec2Dxz",sBrkts);
	MlMatrOut (fMOut, merged, adBxz, iNPt, iNPz, 'd');
	free(merged);
	merged=sMerge("Spec2Dyz",sBrkts);
	MlMatrOut (fMOut, merged, adByz, iNPp, iNPz, 'd');
	free(merged);
    free(adBxz);
	free(adBxy);
    free(adByz);
	/* --------- 2D averaging (end) --------*/

  }/* iCQ for */

  fclose(fMOut);

#ifdef TR_DR_PHOT
	if (iFMC==0)
    sprintf (sFN_NTrv,"DREvents_step%i",step);
  else
    sprintf (sFN_NTrv,"DREvents_fm%i_step%i",iFMC,step);
  fMOut= MFOpen (sFN_NTrv, pfade->RESULT, -1);

   fprintf (fMOut, "iFMC = %i;\n",iFMC);
  if (iFMC==0)
    fprintf (fMOut,"iflFM=0;  %% Normal mesh\n%%\n");
  else
    fprintf (fMOut,"iflFM=1;  %% Fine mesh\n%%\n");

  fprintf (fMOut, "FileType = 'DissRec';\n%%\n");

  fprintf (fMOut, "iNEl = 4;\n%%\n");
	
  fprintf (fMOut, "dx= %.2e;\n",dx);
  fprintf (fMOut, "dy= %.2e;\n",dy);
  fprintf (fMOut, "dz= %.2e;\n",dz);
  fprintf (fMOut, "%%\n");

  fprintf (fMOut, "iNPx= %i;\n",iNPt);
  fprintf (fMOut, "iNPy= %i;\n",iNPp);
  fprintf (fMOut, "iNPz= %i;\n",iNPz);
  fprintf (fMOut, "%%\n");

  fprintf(fMOut, "iTimeStepNo=%d;\n",step);
  fprintf(fMOut, "dStepTime=%f;\n",time_c->step); 
  fprintf(fMOut, "dTotalTime=%f;\n",time_c->total);
  fprintf(fMOut, "dRealTime=%f;\n",time_c->dRealTm);
  fprintf (fMOut, "%%\n");

  /* ---- Writing arrays of X,Y,Z . . . ---- */
  for (iCPt=0; iCPt<iNPt; iCPt++) Buf[iCPt]=dMinX+iCPt*dx;
  MlMatrOut (fMOut, "dPx", Buf, iNPt, 1, 'd');
	
  for (iCPp=0; iCPp<iNPp; iCPp++) Buf[iCPp]=dMinY+iCPp*dy;
  MlMatrOut (fMOut, "dPy", Buf, iNPp, 1, 'd');

  for (iCPz=0; iCPz<iNPz; iCPz++) Buf[iCPz]=iCPz*dz;
  MlMatrOut (fMOut, "dPz", Buf, iNPz, 1, 'd');

  for (iCPp=0; iCPp<iNPp; iCPp++)
  for (iCPt=0; iCPt<iNPt; iCPt++)
	Buf[iCPt+iCPp*iNPt] = v_all->b[iCPt][iCPp].z_min;
  MlMatrOut (fMOut, "dLimZxy", Buf, iNPt, iNPp, 'd');

  
  // Write "actual" output:
  for(iCHx=1; iCHx <= 4; iCHx++){  /* 4 things to trace: CH2+, CH3+ and CH4+ as well as total DR photons */
	switch(iCHx) {
		case 1:
		sCEl = "_ch2";
		break;
		case 2:
		sCEl = "_ch3";
		break;
		case 3:
		sCEl = "_ch4";
		break;
		case 4:
		sCEl = "_SUM";
		break;
		default:
		printf("Fatal error (MatlOut.c)! Cannot print more than 4 DR types!\n");
		exit(-1);
		}
    strcpy (sElN,"___________");
    sprintf(sElN_, "DR%s(+1)", sCEl);
	strncpy (sElN, sElN_, strlen (sElN_));
	fprintf (fMOut, "%%\n%% -------------------- Element #%i '%s' -------------------\n",iCHx,sElN);
    fprintf (fMOut, "ElName(%i,:) = '%s';\n%%\n",iCHx,sElN);

	sprintf (sBrkts,"(%i,:,:)",iCHx);

    /* --------- 2D averaging ------------*/
	aiBxz=(int*)calloc(iNPt*iNPz,sizeof(int));
    aiByz=(int*)calloc(iNPp*iNPz,sizeof(int));
	aiBxy=(int*)calloc(iNPt*iNPp,sizeof(int));
    for (iCPt=0; iCPt<iNPt; iCPt++)
	for (iCPp=0; iCPp<iNPp; iCPp++) { 
        iZCor = (int) ceil (v_all->b[iCPt][iCPp].z_min/dz);
        for (iCPz=0; iCPz<iNPz; iCPz++){

		  if(iCHx==4) //print sum of DR rates
			  iCVal = v_all->b[iCPt][iCPp].vcell[iCPz]
					  .iNActDR[0];
		  else
			  iCVal = v_all->b[iCPt][iCPp].vcell[iCPz]
					  .iNActDR[iCHx];

          iCZ = iCPz + iZCor;
	      aiBxy[iCPt+iNPt*iCPp]+= iCVal;
		    if (iCZ<iNPz){
		  aiBxz[iCPt+iNPt*iCZ] += iCVal;
          aiByz[iCPp+iNPp*iCZ] += iCVal;
			}
        }
	}

	/* --------- 2D averaging ------------*/
    merged = sMerge("Dens2Dxy",sBrkts);
	MlMatrOut (fMOut, merged, aiBxy, iNPt, iNPp, 'i');
	free(merged);
    merged = sMerge("Dens2Dxz",sBrkts);
	MlMatrOut (fMOut, merged, aiBxz, iNPt, iNPz, 'i');
	free(merged);
    merged = sMerge("Dens2Dyz",sBrkts);
	MlMatrOut (fMOut, merged, aiByz, iNPp, iNPz, 'i');
	free(merged);
    free(aiBxz);
	free(aiBxy);
    free(aiByz);
    /* --------- 2D averaging (end) --------*/

	adBxz=(double*)calloc(iNPt*iNPz,sizeof(double));
    adByz=(double*)calloc(iNPp*iNPz,sizeof(double));
	adBxy=(double*)calloc(iNPt*iNPp,sizeof(double));
    for (iCPt=0; iCPt<iNPt; iCPt++)
	for (iCPp=0; iCPp<iNPp; iCPp++) { 
        iZCor = (int) ceil (v_all->b[iCPt][iCPp].z_min/dz);
        for (iCPz=0; iCPz<iNPz; iCPz++){

	  if (iCHx==4) //print sum of DR rates
		  dCVal = v_all->b[iCPt][iCPp].vcell[iCPz]
						  .dDRAtoms[0] / time_c->step;
	  else
		  dCVal = v_all->b[iCPt][iCPp].vcell[iCPz]
						  .dDRAtoms[iCHx] / time_c->step;

          iCZ = iCPz + iZCor;
	      adBxy[iCPt+iNPt*iCPp]+= dCVal;
		    if (iCZ<iNPz){
		  adBxz[iCPt+iNPt*iCZ] += dCVal;
          adByz[iCPp+iNPp*iCZ] += dCVal;
			}
        }
	}
    merged=sMerge("Spec2Dxy",sBrkts);
	MlMatrOut (fMOut, merged, adBxy, iNPt, iNPp, 'd');
	free(merged);
	merged=sMerge("Spec2Dxz",sBrkts);
	MlMatrOut (fMOut, merged, adBxz, iNPt, iNPz, 'd');
	free(merged);
	merged=sMerge("Spec2Dyz",sBrkts);
	MlMatrOut (fMOut, merged, adByz, iNPp, iNPz, 'd');
	free(merged);
    free(adBxz);
	free(adBxy);
    free(adByz);
	/* --------- 2D averaging (end) --------*/


#ifdef PilotPSI
	// Output also in r/z coordinates, so that it can be compared with Abel-inverted values.
	if(GenOpt->cfl_print_rz_emission==1){
		iNTheta = 300;  /* number of integration points along theta direction */
		printf("Writing emission in rz coordinates... This only works for linear machines!!\n");
		iNPr = (iNPp<iNPt ? iNPp/2 : iNPt/2) - 1;
		adBrz=(double*)calloc( iNPr * iNPz,sizeof(double));
		for (iCPz=0; iCPz<iNPz; iCPz++) {
			for (iCPr=0; iCPr<iNPr; iCPr++) { 
			   for (iCTheta=0; iCTheta<iNTheta; iCTheta++){
					iCPt = floor( 0.5*iNPt + iCPr * sin(iCTheta * 2 * M_PI / iNTheta ) );
					iCPp = floor( 0.5*iNPp + iCPr * cos(iCTheta * 2 * M_PI / iNTheta ) );
					if (iCHx==4) //print sum of DR rates
						dCVal = v_all->b[iCPt][iCPp].vcell[iCPz]
									  .dDRAtoms[0] / (time_c->step * iNTheta);
					else
						dCVal = v_all->b[iCPt][iCPp].vcell[iCPz]
									  .dDRAtoms[iCHx] / (time_c->step * iNTheta);
           		
					adBrz[iCPr+iNPr*iCPz] += dCVal;
				}
			}
		}
		merged=sMerge("Spec2Drz",sBrkts);
		MlMatrOut (fMOut, merged, adBrz, iNPr, iNPz, 'd');
		free(adBrz);
		free(merged);
	}
#endif

  }
#endif /* TR_DR_PHOT */

  free(Buf);
}
#endif /* TR_ION_ACT */

/* -------------------------------------------------------------
    Main file for   density / emission output . . .
	intergrated 2D planes.
   ------------------------------------------------------------- */
void MatlEmissOut (struct VolumeAll  *v_all,
		int iFMC,
		struct Path	  *pfade,
		struct TimeCell   *time_c,
		int                step,
		struct LostParticle    *looser,
		struct NTRL_COUNT      *strNtrlC,
		struct strGenOptions *GenOpt)
{
  double dx          = v_all->diff[0];
  double dy          = v_all->diff[1];
  double dz          = v_all->diff[2];
  double dMinX       = v_all->adRange[0];
  double dMinY       = v_all->adRange[1];    

  /*  1/cm^3 (Warning: diff[i] should be in mm! */
  double dNDens3      = 1000. / (dx * dy * dz * time_c->step);

  int iCPt, iNPt       = v_all->anz[0];    /* points in toroidal dir. */
  int iCPp, iNPp       = v_all->anz[1];    /* points in poloidal dir. */
  int iCPz, iNPz       = v_all->anz[2];    /* points in vertical dir. */
#ifdef PilotPSI
  int iCPr, iNPr; 						   /* for printing output in r/z coordinates */
  int iCTheta, iNTheta;                    /* for printing output in r/z coordinates */
#endif
  int iCEl, iNEl       = v_all->anz[3];    /* number of elements */


  char sFN_NTrv[DEF_STR_LEN], sBrkts[DEF_STR_LEN];
  char sElN[] ="____";

  FILE *fMOut;

  double *Buf;
  double *adBxy, *adBxz, *adByz;
#ifdef PilotPSI
  double *adBrz;
#endif
  int iZCor, iCZ;  /* Corrected Z coordinates */
  double dCVal;
  char *merged;

#ifdef JETLim
  double dJI_Tot, dJI_Loc, dJD_Loc, dJI_LocOS /* opposite side . . . */;
  double dJI_Ne, dJI_Te;
  int    iJI_NC;
JETSpotVAR
#endif  /* JETLim */

  Buf=(double*)calloc(iNPt*iMax2(iNPp,iNPz),sizeof(double));
  
  /* ---- Opening file, writing trivial things . . . ---- */
  if (iFMC==0)
    sprintf (sFN_NTrv,"Emission_step%i",step);
  else
    sprintf (sFN_NTrv,"Emission_fm%i_step%i",iFMC,step);
  fMOut= MFOpen (sFN_NTrv, pfade->RESULT, -1);

#ifdef JETLim
  /* ---- Opening file, writing trivial things . . . ---- */
  if (iFMC==0)
    sprintf (sFN_NTrv,"EmissJET_step%i",step);
  else
    sprintf (sFN_NTrv,"EmissJET_fm%i_step%i",iFMC,step);
  fJET = MFOpen (sFN_NTrv, pfade->RESULT, -1);
#endif
   
  fprintf (fMOut, "iFMC = %i;\n",iFMC);
  if (iFMC==0)
    fprintf (fMOut,"iflFM=0;  %% Normal mesh\n%%\n");
  else
    fprintf (fMOut,"iflFM=1;  %% Fine mesh\n%%\n");

  fprintf (fMOut, "iNEl = %i;\n%%\n",iNEl);

  fprintf (fMOut, "%% Statistics for collisions with neutrals . . .\n");
  fprintf (fMOut, "%% ---------------------------------------------\n");
  fprintf (fMOut, "Ntrl.NoCol=   %6i;\n",strNtrlC->NoCol);
  fprintf (fMOut, "Ntrl.OneCol=  %6i;\n",strNtrlC->OneCol);
  fprintf (fMOut, "Ntrl.L10Col=  %6i;\n",strNtrlC->L10Col);
  fprintf (fMOut, "Ntrl.L100Col= %6i;\n",strNtrlC->L100Col);
  fprintf (fMOut, "Ntrl.ManyCol= %6i;\n",strNtrlC->ManyCol);
  fprintf (fMOut, "%%\n");
	
  fprintf (fMOut, "dx= %.2e;\n",dx);
  fprintf (fMOut, "dy= %.2e;\n",dy);
  fprintf (fMOut, "dz= %.2e;\n",dz);
  fprintf (fMOut, "%%\n");

  fprintf (fMOut, "iNPx= %i;\n",iNPt);
  fprintf (fMOut, "iNPy= %i;\n",iNPp);
  fprintf (fMOut, "iNPz= %i;\n",iNPz);
  fprintf (fMOut, "%%\n");

  fprintf(fMOut, "iTimeStepNo=%d;\n",step);
  fprintf(fMOut, "dStepTime=%f;\n",time_c->step); 
  fprintf(fMOut, "dTotalTime=%f;\n",time_c->total);
  fprintf(fMOut, "dRealTime=%f;\n",time_c->dRealTm);
  fprintf (fMOut, "%%\n");


  /* ---- Writing arrays of X,Y,Z . . . ---- */
  for (iCPt=0; iCPt<iNPt; iCPt++) Buf[iCPt]=dMinX+iCPt*dx;
  MlMatrOut (fMOut, "dPx", Buf, iNPt, 1, 'd');
	
  for (iCPp=0; iCPp<iNPp; iCPp++) Buf[iCPp]=dMinY+iCPp*dy;
  MlMatrOut (fMOut, "dPy", Buf, iNPp, 1, 'd');

  for (iCPz=0; iCPz<iNPz; iCPz++) Buf[iCPz]=iCPz*dz;
  MlMatrOut (fMOut, "dPz", Buf, iNPz, 1, 'd');

  for (iCPp=0; iCPp<iNPp; iCPp++)
  for (iCPt=0; iCPt<iNPt; iCPt++)
	Buf[iCPt+iCPp*iNPt] = v_all->b[iCPt][iCPp].z_min;
  MlMatrOut (fMOut, "dLimZxy", Buf, iNPt, iNPp, 'd');
  /* ----  X,Y,Z ---- */

  /* for by Element */
  for (iCEl=0; iCEl<iNEl; iCEl++){
    strcpy (sElN,"____");
	strncpy (sElN, v_all->name[iCEl], strlen (v_all->name[iCEl]));
	printf ("Writing emission of '%s' . . . ", sElN);
	fprintf (fMOut, "%%\n%% -------------------- Element '%s' -------------------\n"
		                                             ,sElN);
    fprintf (fMOut, "ElName(%i,:) = '%s (+%i)';\n",iCEl+1,sElN, v_all->charge[iCEl]);
	fprintf (fMOut, "El_NmPure{%i} = '%s';\n",iCEl+1,      v_all->name[iCEl]);
    fprintf (fMOut, "El_Charge(%i) = %i; %% [A]\n",iCEl+1, v_all->charge[iCEl]);
	fprintf (fMOut, "El_WL(%i) = %.2f; %% [A]\n",iCEl+1,     v_all->dWL[iCEl]);
#ifdef METAST
	fprintf (fMOut, "El_MS(%i) = %i;\n%%\n",iCEl+1, v_all->iMS[iCEl]);
#endif


#ifdef JETLim
   dJI_Tot=0.; dJI_Loc=0.; dJI_LocOS=0.; dJD_Loc=0.;
   dJI_Te =0;  dJI_Ne =0; 
   iJI_NC=0;
	/* --------- Special integration in the line of sight (JET Be Limiter) ------------*/
  if (GenOpt->cflMatlSpecOut){
	dJI_Tot =0;  dJI_Loc =0;
    for (iCPt=0; iCPt<iNPt; iCPt++)
	for (iCPp=0; iCPp<iNPp; iCPp++) { 
          iZCor = (int) ceil (v_all->b[iCPt][iCPp].z_min/dz);
          for (iCPz=0; iCPz<iNPz; iCPz++){
	    dCVal = v_all->b[iCPt][iCPp].vcell[iCPz].dAtomEmiss[iCEl];
            iCZ = iCPz + iZCor;
	    dJI_Tot += dCVal;
	  
        dX = dMinX+iCPt*dx;
	    dY = dMinY+iCPp*dy;
            dZ=iCPz*dz;
            /* Making a projection on perpendicular cylinder crossscection . . . */
        dX0_ = dX0 - (dZ-dZ0)*dTanA;
	    dX = (dX - dX0_ )*dCosA;
	    if ((dX*dX + dY*dY)<dRSp2){                    /* Cylindric spot with R=60mm */
               dJI_Loc += dCVal;
			   dJD_Loc += v_all->b[iCPt][iCPp].vcell[iCPz].atomps[iCEl];
			if (dZ<=77.){   /* Bcmt: most BeI/BeII does not penetrate deeper . . . abitrary a bit . . . */
               iJI_NC++;
			   dJI_Ne  += v_all->b[iCPt][iCPp].vcell[iCPz].dNe;
		       dJI_Te  += v_all->b[iCPt][iCPp].vcell[iCPz].dTe;
			}
		}
	    /* BCmt: only for testing . . .  
	    else 
	      v_all->b[iCPt][iCPp].vcell[iCPz].dAtomEmiss[iCEl] =0.; 
            */

		dX = -(dMinX+iCPt*dx);  /* just symmetric case . . . */
        dX0_ = dX0 - (dZ-dZ0)*dTanA;
	    dX = (dX - dX0_ )*dCosA;
	    if ((dX*dX + dY*dY)<3600.)
               dJI_LocOS += dCVal;
          }
        }
	fprintf(fMOut, "%% -------- JET -------\n");
	fprintf(fMOut, "JETIntI_Tot(%i) = %.3e;\n", iCEl+1, dJI_Tot);
	fprintf(fMOut, "JETIntI_Loc(%i) = %.3e;\n", iCEl+1, dJI_Loc);
	fprintf(fMOut, "JETIntD_Loc(%i) = %.3e;\n", iCEl+1, dJD_Loc);
	fprintf(fMOut, "JETIntI_LocOS(%i) = %.3e;\n", iCEl+1, dJI_LocOS);

	if(iCEl==0){
        fprintf (fJET, "iVNCell   = %i;\n", iJI_NC);
		fprintf (fJET, "dTe_VAver = %.3f;  %% [eV]\n", dJI_Te/iJI_NC);
		fprintf (fJET, "dNe_VAver = %.2e;  %% [cm-3]\n",  dJI_Ne/iJI_NC);
	}

	fprintf (fJET, "%%\n%% -------------------- Element '%s' -------------------\n",sElN);
    fprintf (fJET, "ElName(%i,:) = '%s (+%i)';\n",iCEl+1,sElN, v_all->charge[iCEl]);
	fprintf (fJET, "El_NmPure{%i} = '%s';\n",iCEl+1,      v_all->name[iCEl]);
    fprintf (fJET, "El_Charge(%i) = %i; %% [A]\n",iCEl+1, v_all->charge[iCEl]);
	fprintf (fJET, "El_WL(%i) = %.2f; %% [A]\n",iCEl+1,     v_all->dWL[iCEl]);
    fprintf (fJET, "JETIntI_Tot(%i) = %.3e;\n", iCEl+1, dJI_Tot);
	fprintf (fJET, "JETIntI_Loc(%i) = %.3e;\n", iCEl+1, dJI_Loc);
	fprintf (fJET, "JETIntD_Loc(%i) = %.3e;\n", iCEl+1, dJD_Loc);
	fprintf (fJET, "JETIntI_LocOS(%i) = %.3e;\n", iCEl+1, dJI_LocOS);
	fprintf (fJET, "%% -----------------------------------------------------------\n");
  }
#endif  /* JETLim */


    /* --------- 2D averaging ------------*/
	adBxz=(double*)calloc(iNPt*iNPz,sizeof(double));
    adByz=(double*)calloc(iNPp*iNPz,sizeof(double));
	adBxy=(double*)calloc(iNPt*iNPp,sizeof(double));
    for (iCPt=0; iCPt<iNPt; iCPt++)
	for (iCPp=0; iCPp<iNPp; iCPp++) { 
        iZCor = (int) ceil (v_all->b[iCPt][iCPp].z_min/dz);
        for (iCPz=0; iCPz<iNPz; iCPz++){
		  dCVal = v_all->b[iCPt][iCPp].vcell[iCPz].atomps[iCEl]
			  * dNDens3; /* [1/cm^3] */
          iCZ = iCPz + iZCor;
	      adBxy[iCPt+iNPt*iCPp]+= dCVal;
		    if (iCZ<iNPz){
		  adBxz[iCPt+iNPt*iCZ] += dCVal;
          adByz[iCPp+iNPp*iCZ] += dCVal;
			}
        }
	}
	sprintf (sBrkts,"(%i,:,:)",iCEl+1); 
	merged=sMerge("Dens2Dxy",sBrkts);
	MlMatrOut (fMOut, merged, adBxy, iNPt, iNPp, 'd');
	free(merged);
	merged=sMerge("Dens2Dxz",sBrkts);
	MlMatrOut (fMOut, merged, adBxz, iNPt, iNPz, 'd');
	free(merged);
	merged=sMerge("Dens2Dyz",sBrkts);
	MlMatrOut (fMOut, merged, adByz, iNPp, iNPz, 'd');
	free(merged);
    free(adBxz);
	free(adBxy);
    free(adByz);

#ifdef PilotPSI
	// Output particle density in r/z coordinates, so that it can be compared with Abel-inverted values.
	if(GenOpt->cfl_print_rz_emission==1){
		iNTheta = 300;  /* number of integration points along theta direction */
		printf("Writing emission in rz coordinates... This only works for linear machines!!\n");
		iNPr = (iNPp<iNPt ? iNPp/2 : iNPt/2) - 1;
		adBrz=(double*)calloc( iNPr * iNPz,sizeof(double));
		for (iCPz=0; iCPz<iNPz; iCPz++) {
			for (iCPr=0; iCPr<iNPr; iCPr++) { 
			   for (iCTheta=0; iCTheta<iNTheta; iCTheta++){
					iCPt = floor( 0.5*iNPt + iCPr * sin(iCTheta * 2 * M_PI / iNTheta ) );
					iCPp = floor( 0.5*iNPp + iCPr * cos(iCTheta * 2 * M_PI / iNTheta ) );
					dCVal = v_all->b[iCPt][iCPp].vcell[iCPz].atomps[iCEl] * dNDens3 / iNTheta;
					adBrz[iCPr+iNPr*iCPz] += dCVal;
				}
			}
		}
		merged=sMerge("Dens2Drz",sBrkts);
		MlMatrOut (fMOut, merged, adBrz, iNPr, iNPz, 'd');
		free(adBrz);
		free(merged);
	}
#endif
    /* --------- 2D averaging (end) --------*/

	/* --------- 2D averaging (Spectroscopy)------------*/
  if (GenOpt->cflMatlSpecOut){
	adBxz=(double*)calloc(iNPt*iNPz,sizeof(double));
    adByz=(double*)calloc(iNPp*iNPz,sizeof(double));
	adBxy=(double*)calloc(iNPt*iNPp,sizeof(double));
    for (iCPt=0; iCPt<iNPt; iCPt++)
	for (iCPp=0; iCPp<iNPp; iCPp++) { 
        iZCor = (int) ceil (v_all->b[iCPt][iCPp].z_min/dz);
        for (iCPz=0; iCPz<iNPz; iCPz++){
	  dCVal = v_all->b[iCPt][iCPp].vcell[iCPz].dAtomEmiss[iCEl];
          iCZ = iCPz + iZCor;
	      adBxy[iCPt+iNPt*iCPp]+= dCVal;
		    if (iCZ<iNPz){
		  adBxz[iCPt+iNPt*iCZ] += dCVal;
          adByz[iCPp+iNPp*iCZ] += dCVal;
			}
        }
	}
	sprintf (sBrkts,"(%i,:,:)",iCEl+1); 
	merged=sMerge("Spec2Dxy",sBrkts);
	MlMatrOut (fMOut, merged, adBxy, iNPt, iNPp, 'd');
	free(merged);
	merged=sMerge("Spec2Dxz",sBrkts);
	MlMatrOut (fMOut, merged, adBxz, iNPt, iNPz, 'd');
	free(merged);
	merged=sMerge("Spec2Dyz",sBrkts);
	MlMatrOut (fMOut, merged, adByz, iNPp, iNPz, 'd');
	free(merged);
    free(adBxz);
	free(adBxy);
    free(adByz);

#ifdef PilotPSI
	// Output photon emission in r/z coordinates, so that it can be compared with Abel-inverted values.
	if(GenOpt->cfl_print_rz_emission==1){
		iNTheta = 300;  /* number of integration points along theta direction */
		printf("Writing emission in rz coordinates... This only works for linear machines!!\n");
		iNPr = (iNPp<iNPt ? iNPp/2 : iNPt/2) - 1;
		adBrz=(double*)calloc( iNPr * iNPz,sizeof(double));
		for (iCPz=0; iCPz<iNPz; iCPz++) {
			for (iCPr=0; iCPr<iNPr; iCPr++) { 
			   for (iCTheta=0; iCTheta<iNTheta; iCTheta++){
					iCPt = floor( 0.5*iNPt + iCPr * sin(iCTheta * 2 * M_PI / iNTheta ) );
					iCPp = floor( 0.5*iNPp + iCPr * cos(iCTheta * 2 * M_PI / iNTheta ) );
					dCVal = v_all->b[iCPt][iCPp].vcell[iCPz].dAtomEmiss[iCEl] / iNTheta;
					adBrz[iCPr+iNPr*iCPz] += dCVal;
				}
			}
		}
		merged=sMerge("Spec2Drz",sBrkts);
		MlMatrOut (fMOut, merged, adBrz, iNPr, iNPz, 'd');
		free(adBrz);
		free(merged);
	}
#endif

  }
    /* --------- 2D averaging for Spectroscopy (end) --------*/


    /* -------------- 3D output ---------------- */
	if (GenOpt->cflMatl3DOut==1)
    for (iCPz=0; iCPz<iNPz; iCPz++){
	  fprintf(fMOut, "PDens(%i,%i,:,:) = [\n", iCEl+1,iCPz+1);
	  for (iCPp=0; iCPp<iNPp; iCPp++){
        for (iCPt=0; iCPt<iNPt; iCPt++) 	
			fprintf(fMOut," %.2e", 
			v_all->b[iCPt][iCPp].vcell[iCPz].atomps[iCEl]
			 * dNDens3 /* [1/cm^3] */);		                    
        fprintf(fMOut,"\n");
	  }
      fprintf(fMOut, "];\n%%\n");
	}

#ifdef DIV_DOPPLERSPEC
	// TAM, Toni Makkonen, sep 2012
	if (GenOpt->cflMatl3DOut==1) {
	  // 3D "spec dens" (emission)
	  for (iCPz=0; iCPz<iNPz; iCPz++){
	    fprintf(fMOut, "SpecDens(%i,%i,:,:) = [\n", iCEl+1,iCPz+1);
	    for (iCPp=0; iCPp<iNPp; iCPp++){
	      for (iCPt=0; iCPt<iNPt; iCPt++) fprintf(fMOut," %.2e", v_all->b[iCPt][iCPp].vcell[iCPz].dAtomEmiss[iCEl]);
	      fprintf(fMOut,"\n");
	    }
	    fprintf(fMOut, "];\n%%\n");
	  }
	  // X
	  for (iCPz=0; iCPz<iNPz; iCPz++){
	    fprintf(fMOut, "AvVelX(%i,%i,:,:) = [\n", iCEl+1,iCPz+1);
	    for (iCPp=0; iCPp<iNPp; iCPp++){
	      for (iCPt=0; iCPt<iNPt; iCPt++) fprintf(fMOut," %.2e", v_all->b[iCPt][iCPp].vcell[iCPz].AvVelX[iCEl]);
	      fprintf(fMOut,"\n");
	    }
	    fprintf(fMOut, "];\n%%\n");
	  }
	  // Y
	  for (iCPz=0; iCPz<iNPz; iCPz++){
	    fprintf(fMOut, "AvVelY(%i,%i,:,:) = [\n", iCEl+1,iCPz+1);
	    for (iCPp=0; iCPp<iNPp; iCPp++){
	      for (iCPt=0; iCPt<iNPt; iCPt++) fprintf(fMOut," %.2e", v_all->b[iCPt][iCPp].vcell[iCPz].AvVelY[iCEl]);
	      fprintf(fMOut,"\n");
	    }
	    fprintf(fMOut, "];\n%%\n");
	  }
	  // Z
	  for (iCPz=0; iCPz<iNPz; iCPz++){
	    fprintf(fMOut, "AvVelZ(%i,%i,:,:) = [\n", iCEl+1,iCPz+1);
	    for (iCPp=0; iCPp<iNPp; iCPp++){
	      for (iCPt=0; iCPt<iNPt; iCPt++) fprintf(fMOut," %.2e", v_all->b[iCPt][iCPp].vcell[iCPz].AvVelZ[iCEl]);
	      fprintf(fMOut,"\n");
	    }
	    fprintf(fMOut, "];\n%%\n");
	  }
	  // XX
	  for (iCPz=0; iCPz<iNPz; iCPz++){
	    fprintf(fMOut, "AvVelXX(%i,%i,:,:) = [\n", iCEl+1,iCPz+1);
	    for (iCPp=0; iCPp<iNPp; iCPp++){
	      for (iCPt=0; iCPt<iNPt; iCPt++) fprintf(fMOut," %.2e", v_all->b[iCPt][iCPp].vcell[iCPz].AvVelXX[iCEl]);
	      fprintf(fMOut,"\n");
	    }
	    fprintf(fMOut, "];\n%%\n");
	  }
	  // YY
	  for (iCPz=0; iCPz<iNPz; iCPz++){
	    fprintf(fMOut, "AvVelYY(%i,%i,:,:) = [\n", iCEl+1,iCPz+1);
	    for (iCPp=0; iCPp<iNPp; iCPp++){
	      for (iCPt=0; iCPt<iNPt; iCPt++) fprintf(fMOut," %.2e", v_all->b[iCPt][iCPp].vcell[iCPz].AvVelYY[iCEl]);
	      fprintf(fMOut,"\n");
	    }
	    fprintf(fMOut, "];\n%%\n");
	  }
	  // ZZ
	  for (iCPz=0; iCPz<iNPz; iCPz++){
	    fprintf(fMOut, "AvVelZZ(%i,%i,:,:) = [\n", iCEl+1,iCPz+1);
	    for (iCPp=0; iCPp<iNPp; iCPp++){
	      for (iCPt=0; iCPt<iNPt; iCPt++) fprintf(fMOut," %.2e", v_all->b[iCPt][iCPp].vcell[iCPz].AvVelZZ[iCEl]);
	      fprintf(fMOut,"\n");
	    }
	    fprintf(fMOut, "];\n%%\n");
	  }
	  // XY
	  for (iCPz=0; iCPz<iNPz; iCPz++){
	    fprintf(fMOut, "AvVelXY(%i,%i,:,:) = [\n", iCEl+1,iCPz+1);
	    for (iCPp=0; iCPp<iNPp; iCPp++){
	      for (iCPt=0; iCPt<iNPt; iCPt++) fprintf(fMOut," %.2e", v_all->b[iCPt][iCPp].vcell[iCPz].AvVelXY[iCEl]);
	      fprintf(fMOut,"\n");
	    }
	    fprintf(fMOut, "];\n%%\n");
	  }
	  // XZ
	  for (iCPz=0; iCPz<iNPz; iCPz++){
	    fprintf(fMOut, "AvVelXZ(%i,%i,:,:) = [\n", iCEl+1,iCPz+1);
	    for (iCPp=0; iCPp<iNPp; iCPp++){
	      for (iCPt=0; iCPt<iNPt; iCPt++) fprintf(fMOut," %.2e", v_all->b[iCPt][iCPp].vcell[iCPz].AvVelXZ[iCEl]);
	      fprintf(fMOut,"\n");
	    }
	    fprintf(fMOut, "];\n%%\n");
	  }
	  // YZ
	  for (iCPz=0; iCPz<iNPz; iCPz++){
	    fprintf(fMOut, "AvVelYZ(%i,%i,:,:) = [\n", iCEl+1,iCPz+1);
	    for (iCPp=0; iCPp<iNPp; iCPp++){
	      for (iCPt=0; iCPt<iNPt; iCPt++) fprintf(fMOut," %.2e", v_all->b[iCPt][iCPp].vcell[iCPz].AvVelYZ[iCEl]);
	      fprintf(fMOut,"\n");
	    }
	    fprintf(fMOut, "];\n%%\n");
	  }

	}
#endif

	/* -------------- 3D output (end) ---------- */
    
	printf ("Complete!\n");
  } /* for by Element */
  

  free(Buf);
  fclose (fMOut);
#ifdef JETLim
  fclose (fJET);
#endif  /* JETLim */
}

/* -------------------------------------------------------
   Writes out a Matlab file containing limiter information
   ------------------------------------------------------- */
void MatlLimiterOut (int                  surf, 
				     struct Ort           *pOrt,
				     int                  *iRand,
				     double		          *dBound,
					 struct SurfaceCell   **net,
					 struct Path	      *pfade,
#ifdef ITER_BM
					 struct strBMPar *pBMP,
#endif
					 struct strGenOptions *GenOpt)
{
   char sFN_NTrv[DEF_STR_LEN];
   FILE *fMOut;

   /* necessary for cells loops */
   int                x_cell, y_cell;         /* laufvariablen                      */
   struct SurfaceCell *y_snet;			      /* zugehoerige laufpointer */
   int iNCx, iNCy;
   double *MdthX, *MdthY, *MdthZ;

#ifdef ITER_BM
#define PRINT_OUT_NORMAL
#endif

#ifdef PRINT_OUT_NORMAL
   double *Norm0, *Norm1, *Norm2;
   double *Br0, *Br1, *Br2;
   struct Magnetic m_ptr;   /* BCmt: for B-field  . . . */
   struct BField   bfield;  /* BCmt: for B-field  . . . */
#endif 

   char sVarN[DEF_STR_LEN];

#ifdef ITER_BM
   double dRangeBM, dxBM;
#endif


   sprintf (sFN_NTrv,"limiter");
   fMOut= MFOpen (sFN_NTrv, pfade->RESULT, -1);
   
   fprintf (fMOut, "iSType = %i;\n", surf);
   fprintf (fMOut, "%%\n");

   fprintf(fMOut,"%% Volume taken by limiter (x1, x2, y1, y2,  dx, dy)\n");
   MlMatrOut (fMOut, "dBound", dBound, N_GRID_BOUND_ELEM, 1, 'd');
   fprintf(fMOut,"%% Number of cells (Ny, Nx, N elements, N marks)\n");
   MlMatrOut (fMOut, "iRand", iRand, 4, 1, 'i');
   fprintf (fMOut, "iNx   = %8i;\n", iRand[1]);
   fprintf (fMOut, "iNy   = %8i;\n", iRand[0]);
   fprintf (fMOut, "%%\n%%\n");


   fprintf(fMOut,"%% Data contained in structure 'Ort'.\n");
   fprintf(fMOut,"%% ----------------------------------\n");
#ifndef UGEOMETRY
   fprintf (fMOut, "RadX  = %.3f;\n", pOrt->t_rad);
   fprintf (fMOut, "LenX  = %.3f;\n", pOrt->t_len);
   fprintf (fMOut, "RadY  = %.3f;\n", pOrt->p_rad);
   fprintf (fMOut, "LenY  = %.3f;\n", pOrt->p_len);
   fprintf (fMOut, "RadZ  = %8i;\n", pOrt->z_rad);
   fprintf (fMOut, "LenZ  = %8i;\n", pOrt->z_len);
   fprintf(fMOut,"%% Angle of main limiter plane.\n");
   fprintf (fMOut, "Ang   = %8.3f;  %% RAD\n", (double)pOrt->angle);
   fprintf (fMOut, "AngG  = %8.3f;  %% Grad\n", (double)pOrt->angle/(RAD));
#endif
   fprintf(fMOut,"%% Gas puffing source data:\n");
   fprintf (fMOut, "SrcX  = %8.2f;\n", pOrt->s_tor);
   fprintf (fMOut, "SrcY  = %8.2f;\n", pOrt->s_pol);
   fprintf (fMOut, "SrcZ  = %8.2f;\n", pOrt->s_rad);
   fprintf (fMOut, "STAng = %8.3f;  %% RAD\n", (double)pOrt->dTorAng);
   fprintf (fMOut, "STAngG= %8.3f;  %% Grad\n", (double)pOrt->dTorAng/(RAD));
   fprintf (fMOut, "SPAng = %8.3f;  %% RAD\n", (double)pOrt->dPolAng);
   fprintf (fMOut, "SPAngG= %8.3f;  %% Grad\n", (double)pOrt->dPolAng/(RAD));
   fprintf (fMOut, "%%\n%%\n");
   /* struct Ort contains some other parameters! */

#ifdef ITER_BM
   fprintf (fMOut, "%% ----------------- ITER_BM -----------------\n");
   fprintf (fMOut, "BM_rho_plasma_tor  = %.4f;\n", pBMP->rho_plasma_tor);
   fprintf (fMOut, "BM_rho_pp  = %.4f;\n", pBMP->rho_pp);
   fprintf (fMOut, "BM_halfdim_tor  = %.4f;\n", pBMP->halfdim_tor);
   fprintf (fMOut, "BM_beta  = %.4f;\n", pBMP->beta);
   fprintf (fMOut, "BM_Zero  = %.2f;  %%[mm]\n", pBMP->dZeroMM);
   fprintf (fMOut, "BM_Bt_x  = %.2f;  %%[mm]\n", pBMP->Bt_x);
   fprintf (fMOut, "BM_Bp_y  = %.2f;  %%[mm]\n", pBMP->Bp_y);
   fprintf (fMOut, "BM_Br_z  = %.2f;  %%[mm]\n", pBMP->Br_z);
#ifndef JETLim
   fprintf (fMOut, "BM_Alp   = %.2f;  %%[Rad] =%.2fGrad\n", pBMP->BM_Alp, pBMP->BM_Alp*GRAD);
   fprintf (fMOut, "BM_treent  = %.4f;\n", pBMP->treent);
   fprintf (fMOut, "BM_slot  = %.4f;\n", pBMP->slot);
   fprintf (fMOut, "BM_RRidge  = %.2f;  %%[mm]\n", pBMP->dRRidge);
#endif

   assert((MdthX = (double*)calloc(1*201,sizeof(double)))!=NULL);
   assert((MdthY = (double*)calloc(1*201,sizeof(double)))!=NULL);
   dRangeBM = pBMP->halfdim_tor*1.1*1000; /* m->mm */
   dxBM = dRangeBM/100;
   for (iNCx=0; iNCx<201; iNCx++){
     MdthX[iNCx] = iNCx*dxBM - dRangeBM;
#ifdef JETLim
     MdthY[iNCx] = BM_dJETLimShape(MdthX[iNCx], 0., pBMP);
#else
	 MdthY[iNCx] = BM_dTorShape(MdthX[iNCx], pBMP);
#endif
   }
   MlMatrOut (fMOut, "BM_X", MdthX, 201, 1, 'd');
   MlMatrOut (fMOut, "BM_Y", MdthY, 201, 1, 'd');
   fprintf (fMOut, "%%\n");
#ifdef JETLim
   dRangeBM = pBMP->dim_pol/2*1.15*1000; /* m->mm */
   dxBM = dRangeBM/100;
   for (iNCx=0; iNCx<201; iNCx++){
     MdthX[iNCx] = iNCx*dxBM - dRangeBM;
     MdthY[iNCx] = BM_dJETLimShape(0., MdthX[iNCx], pBMP);
   }
   MlMatrOut (fMOut, "JETPolCS_X", MdthX, 201, 1, 'd');
   MlMatrOut (fMOut, "JETPolCS_Y", MdthY, 201, 1, 'd');
#endif
   free(MdthX);
   free(MdthY);
   fprintf (fMOut, "%%\n%%\n");
#endif

   /* ------------------- CELLs --------------------------*/
   fprintf(fMOut,"%% Surface cells output.\n");
   fprintf(fMOut,"%% ----------------------------------\n");
   
   iNCx=iRand[1]; iNCy=iRand[0];
   MdthX = (double*)calloc(iNCx*iNCy,sizeof(double));
   MdthY = (double*)calloc(iNCx*iNCy,sizeof(double));
   MdthZ = (double*)calloc(iNCx*iNCy,sizeof(double));
#ifdef PRINT_OUT_NORMAL 
   Norm0 = (double*)calloc(iNCx*iNCy,sizeof(double));
   Norm1 = (double*)calloc(iNCx*iNCy,sizeof(double));
   Norm2 = (double*)calloc(iNCx*iNCy,sizeof(double));
   Br0 = (double*)calloc(iNCx*iNCy,sizeof(double));
   Br1 = (double*)calloc(iNCx*iNCy,sizeof(double));
   Br2 = (double*)calloc(iNCx*iNCy,sizeof(double));
#endif /* PRINT_OUT_NORMAL */
  
   for(x_cell=0;x_cell<iNCx;x_cell++){

	 for(y_cell=0;y_cell<iNCy;y_cell++){
	   	 y_snet = &(net[x_cell][y_cell]);/*die richtige x-reihe des oberflachennetzes einsetzen*/
	   MdthX[x_cell + y_cell*iNCx] =	y_snet->midth[0];
	   MdthY[x_cell + y_cell*iNCx] =	y_snet->midth[1];
	   MdthZ[x_cell + y_cell*iNCx] =	y_snet->midth[2];

#ifdef PRINT_OUT_NORMAL 
	   Norm0[x_cell + y_cell*iNCx] =	y_snet->normal[0];
	   Norm1[x_cell + y_cell*iNCx] =	y_snet->normal[1];
	   Norm2[x_cell + y_cell*iNCx] =	y_snet->normal[2];

	   BM_B_field(&m_ptr, &bfield,y_snet->midth[0],y_snet->midth[1],y_snet->midth[2]);
	   Br0[x_cell + y_cell*iNCx] =	bfield.Br[0];
	   Br1[x_cell + y_cell*iNCx] =	bfield.Br[1];
	   Br2[x_cell + y_cell*iNCx] =	bfield.Br[2];
#endif /* PRINT_OUT_NORMAL */

	 }/*y_cell*/
   
   }/*x_cell*/
   sprintf(sVarN,"Cell.mdX");
   MlMatrOut (fMOut, sVarN, MdthX, iNCx, iNCy, 'd');
   sprintf(sVarN,"Cell.mdY");
   MlMatrOut (fMOut, sVarN, MdthY, iNCx, iNCy, 'd');
   sprintf(sVarN,"Cell.mdZ");
   MlMatrOut (fMOut, sVarN, MdthZ, iNCx, iNCy, 'd');

   free(MdthX);
   free(MdthY);
   free(MdthZ);

#ifdef PRINT_OUT_NORMAL 
   sprintf(sVarN,"Cell.Nrm0");
   MlMatrOut (fMOut, sVarN, Norm0, iNCx, iNCy, 'd');
   sprintf(sVarN,"Cell.Nrm1");
   MlMatrOut (fMOut, sVarN, Norm1, iNCx, iNCy, 'd');
   sprintf(sVarN,"Cell.Nrm2");
   MlMatrOut (fMOut, sVarN, Norm2, iNCx, iNCy, 'd');
   free(Norm0);
   free(Norm1);
   free(Norm2);

   sprintf(sVarN,"Cell.Br0");
   MlMatrOut (fMOut, sVarN, Br0, iNCx, iNCy, 'd');
   sprintf(sVarN,"Cell.Br1");
   MlMatrOut (fMOut, sVarN, Br1, iNCx, iNCy, 'd');
   sprintf(sVarN,"Cell.Br2");
   MlMatrOut (fMOut, sVarN, Br2, iNCx, iNCy, 'd');
   free(Br0);
   free(Br1);
   free(Br2);
#endif /* PRINT_OUT_NORMAL */
   /* ------------------- CELLs (end)----------------------*/


   fclose (fMOut);
   printf ("Limiter data output complete.\n");
}

#ifdef JETLim
/*--BCmt----------------------------------------------------------------------------------------*/
/*                        Integrated values inside JET Be limiter 'spot' (PFMC-2013)            */
/*--BCmt----------------------------------------------------------------------------------------*/
void JETSpoISOut( struct SurfaceCell     **netz,   /*heisst im move_ctrl() "net"                */
		     struct Path            *pfades,  /*heisst im move_ctrl() "pfade"              */
		     struct TimeCell         *timer,  /*heisst im move_ctrl() "time_ptr"           */                
		     int                      *rand,  /*heisst im move_ctrl() "grid_rand"          */
		     int                       step,  /*heisst im move_ctrl() "time_g"             */
			 struct PhErData        *strPhErData
)
{
  double dY_BeD=0., dY_BeBe=0.;
  double dTeA=0., dNeA=0., dHFlA=0., dArea=0., dBAngA=0., dPER=0., dSPER=0., dPERPl=0;
JETSpotVAR
  int                x_cell;         /* laufvariablen                      */
  int                y_cell;         /* laufvariablen                      */
  struct SurfaceCell *y_snet;        /* zugehoerige laufpointer im sputter */
  int counter;
  char sFN_NTrv[DEF_STR_LEN];		   /* Filename part */   
  struct SurfaceCell  strSCell;
  int i, iBB, iCN =0;

	/* BCMt: ----------------- Integrate values inside JET Be limiter 'spot' (PFMC-2013) --------------- */

iBB = 0;    /* flag fro Be self-sputter */
for (i=0; i<strPhErData->iNPIComb; i++)
	if (strcmp(strPhErData->sComb[i],"bebe")==0)
	{iBB = 1; break;}

SURF_START(rand[0],rand[1],y_snet,netz);
        dX = y_snet->midth[0];
	    dY = y_snet->midth[1];
        dZ = y_snet->midth[2];
            /* Making a projection on perpendicular cylinder crossscection . . . */
        dX0_ = dX0 - (dZ-dZ0)*dTanA;
	    dX = (dX - dX0_ )*dCosA;
	    if ((dX*dX + dY*dY)<dRSp2){                    /* Cylindric spot with R=60mm */
			   iCN ++;
			   dArea += y_snet->area;
               dY_BeD  +=   PlasmaSputterLogic(  "_d", "be",  strPhErData,       
						3*y_snet->Te,     /* BCmt: Ignored inside! */
						y_snet
					  )*y_snet->area;
			 
			   if (iBB>0)
               dY_BeBe +=   PlasmaSputterLogic(  "be", "be",  strPhErData,       
						3*y_snet->Te,    /* BCmt: Ignored inside! */
						y_snet
					  )*y_snet->area;
			   dHFlA   +=  y_snet->H_fluence;
			   dNeA    +=  y_snet->ne*y_snet->area;
			   dTeA    +=  y_snet->Te*y_snet->area;
			   dBAngA  +=  acos(y_snet->b_angle)*y_snet->area;
			   dPER    +=  y_snet->energetic[0];
			   dSPER   +=  y_snet->content[0].ngself;
			   dPERPl  +=  y_snet->content[0].dPEr_Plasma;
			   if (strPhErData->iShadowMod==0)
			        y_snet->iShadow *= -1;
               /* BCmt: misusing shadow to light out spot position . . . */ 
		}
SURF_ENDE ;

  sprintf (sFN_NTrv,"SurfJET");
  fJET = MFOpen (sFN_NTrv, pfades->RESULT, step);
	fprintf (fJET, "%%\n%% -------------------- Element '%s' -------------------\n",y_snet->content[0].element);
    fprintf (fJET, "iNC = %i; \n", iCN);
    fprintf (fJET, "SpotArea = %.3f;   %%[mm^2]\n", dArea);
    fprintf (fJET, "dY_Be_D  = %.3e;\n", dY_BeD/dArea);
	fprintf (fJET, "dY_BeBe  = %.3e;\n", dY_BeBe/dArea);
	fprintf (fJET, "dHFlux   = %.3e;\n", dHFlA/dArea);
	fprintf (fJET, "Av_Ne    = %.3e;\n", dNeA/dArea);	
	fprintf (fJET, "Av_Te    = %.3f;\n", dTeA/dArea);	
	fprintf (fJET, "dBAng    = %.3e;\n", dBAngA/dArea);	
	fprintf (fJET, "dBAngG   = %.2f;\n", dBAngA/dArea*GRAD);	
	fprintf (fJET, "dPER     = %.3e;\n", dPER/dArea);	
	fprintf (fJET, "dSPER    = %.3e;\n", dSPER/dArea);
    fprintf (fJET, "dPER_Pl  = %.3e;\n", dPERPl/dArea);
	strSCell.b_angle = cos(dBAngA/dArea);
	strSCell.Te = dTeA/dArea;
	fprintf (fJET, "Yld_Be_D_Aver    = %.3e;\n", 
		         PlasmaSputterLogic(  "_d", "be",  strPhErData,       
						3*y_snet->Te,     /* BCmt: Ignored inside! */
						&strSCell));
	strSCell.b_angle = 1.0;
	fprintf (fJET, "Yld_Be_D_NI      = %.3e;\n", 
		         PlasmaSputterLogic(  "_d", "be",  strPhErData,       
						3*y_snet->Te,     /* BCmt: Ignored inside! */
						&strSCell));
	fprintf (fJET, "%% -----------------------------------------------------------\n");
  fclose (fJET);
}
#endif  /* JETLim */

/*--SDcmt----------------------------------------------------------------------------------------*/
/*                                                                                               */
/* Write the available surface information into one matlab output file                           */
/*                                                                                               */
/*--SDcmt----------------------------------------------------------------------------------------*/
void MatlSurfaceOut( struct SurfaceCell     **netz,   /*heisst im move_ctrl() "net"                */
		     struct Path            *pfades,  /*heisst im move_ctrl() "pfade"              */
		     struct TimeCell         *timer,  /*heisst im move_ctrl() "time_ptr"           */
		     struct LostParticle    *looser,  /*heisst im move_ctrl() "lost"               */                 
		     struct SurfaceGridData *grid_d,  /*heisst im move_ctrl() "sgrid" (aus main.c) */   
		     int                      *rand,  /*heisst im move_ctrl() "grid_rand"          */
		     double                  *bound,  /*heisst im move_ctrl() "grid_bound"         */
		     int                       step,  /*heisst im move_ctrl() "time_g"             */
		     int                       surf   /* Surface type */,

		   		     int iMain_seed   )
{
  /* ----------------------------------------------------------------------
     Variablendefinition
     -------------------------------------------------------------------------*/
  /* beginning of each file             */
  int                elemente;       /* oldy but goody - how many elements */
  int                elem_no;        /* laufvariablen                      */
  int                x_cell;         /* laufvariablen                      */
  int                y_cell;         /* laufvariablen                      */

  struct SurfaceCell *y_snet;        /* zugehoerige laufpointer im sputter */
                                     /* zyklus                             */

#ifdef HardSpSeq
  const int iActNElem=rand[2];       /* Actual number of lements traced    */
  int       *piHES;                  /* Hard element sequence */
  int       iCEl;                    /* Loop var for element */
  int       iCL;                     /* Loop for layer */
  double    *dBuf;
#endif

#ifdef BE_CARBYDE
  double *dBufBC;
#endif

  double             nada;           /* hilfsvariable bei ausgabe          */	

  char sFN_NTrv[DEF_STR_LEN];		   /* Filename part */   
  char sElementName[] ="____";	   /* Fixed length Element name */

  /* Number of different particle "statuses": */	

#define PARTSTATNO 18

  char sParticleStatus[PARTSTATNO][8];/* Explanation of the particle status*/	
  int iStatusC;					    /* Particle status counter */

  int iCP;
  int iNCx, iNCy;
  double *MdthX, *MdthY;
#ifdef UGEOMETRY
  double *MdthZ;
  int *nVertices;
  double **xVert;
  double **yVert;
  double **zVert;
  int iv;
#endif
  int *iClX, *iClY;
  double *SurfT, *Cellsize,*Area, *ShPot, *H_flc;
  double *pTi, *pTe, *pNe;
#ifdef ITER_BM
  int *Shadow;
#endif
#ifdef JETLim
  double *CLen;
#endif
#ifdef BE_CARBYDE
  double *pdCrb;
#endif
  
  
  char sVarN[DEF_STR_LEN];
  int counter;


  FILE *fMOut, *f_FILMTMP;
  double dLostT[MAX_N_ELEM];
	
  int mem_ctrl,merke;
#ifndef HardSpSeq
  int k;
#endif
  struct Spezies CijSpez;
  /* ----------------------------------------------------------------------
     endof Variablendefinition
     -------------------------------------------------------------------------*/


  iNCx=rand[1]; iNCy=rand[0];

  elemente = rand[2];
	
  /*SDcmt: Initialize array of status strings:*/
  strcpy(sParticleStatus[0] ,"NG     ");
  strcpy(sParticleStatus[1] ,"NC     ");
  strcpy(sParticleStatus[2] ,"_T     ");
  strcpy(sParticleStatus[3] ,"_R     ");
  strcpy(sParticleStatus[4] ,"ME     ");
  strcpy(sParticleStatus[5] ,"MQ     ");
  strcpy(sParticleStatus[6] ,"SPER   ");
  strcpy(sParticleStatus[7] ,"UC     ");
  strcpy(sParticleStatus[8] ,"PQ     ");
  strcpy(sParticleStatus[9] ,"PD     ");
  strcpy(sParticleStatus[10],"PD_CHEM");
  strcpy(sParticleStatus[11],"BG     ");	
  strcpy(sParticleStatus[12],"SH     ");
  strcpy(sParticleStatus[13],"FLC    ");
  strcpy(sParticleStatus[14],"CER    ");
  strcpy(sParticleStatus[15],"PER    ");
  /* strcpy(sParticleStatus[16],"PER6   "); BCmt: doubling SZ . . . */
  strcpy(sParticleStatus[16],"PER(D) ");
  strcpy(sParticleStatus[17],"BAng   ");

  /*SDcmt: Open File for writing (w+)*/
  sprintf (sFN_NTrv,"Surface_step%i",step);
  fMOut= MFOpen (sFN_NTrv, pfades->RESULT, -1);
	
  /*SDcmt: Write infos valid for all elements and particle statuses */
  fprintf(fMOut, "%%\n");
  fprintf(fMOut, "%% Number of cell features\n");
  fprintf(fMOut, "iNClFtr=%i;\n",PARTSTATNO);
  fprintf(fMOut, "%% Number of elements followed\n");
  fprintf(fMOut, "iNElem=%i;\n",elemente);
  fprintf(fMOut, "%%\n");
  fprintf(fMOut, "%% -----------------------------------------------------:\n");
  fprintf(fMOut, "%% Informations valid for all elements and particle statuses:\n");
  fprintf(fMOut, "iSType=%d;\n",surf);
  if (timer!=NULL){  
    fprintf(fMOut, "iTimeStepNo=%d;\n",step);
    fprintf(fMOut, "dStepTime=%f;\n",timer->step); 
	fprintf(fMOut, "dRealTime=%f;\n",timer->dRealTm); 
    fprintf(fMOut, "dTotalTime=%f;\n",timer->total);}

	
  fprintf(fMOut, "dXmin= %.1f;\n",bound[0]);
  fprintf(fMOut, "dXmax= %.1f;\n",bound[1]);
  fprintf(fMOut, "dYmin= %.1f;\n",bound[2]);
  fprintf(fMOut, "dYmax= %.1f;\n",bound[3]);
#ifdef UGEOMETRY
  fprintf(fMOut, "dZmin= %.1f;\n",bound[4]);
  fprintf(fMOut, "dZmax= %.1f;\n",bound[5]);
  fprintf(fMOut, "dDeltaX=%.1f;\n",bound[6]);
  fprintf(fMOut, "dDeltaY=%.1f;\n",bound[7]);
  fprintf(fMOut, "dDeltaZ=%.1f;\n",bound[8]);
#else
  fprintf(fMOut, "dDeltaX=%.1f;\n",bound[4]);
  fprintf(fMOut, "dDeltaY=%.1f;\n",bound[5]);
#endif

  fprintf(fMOut, "iNX=%d;\n",iNCx);
  fprintf(fMOut, "iNY=%d;\n",iNCy);

  fprintf(fMOut, "Random_Seed=%d;\n",iMain_seed);

  nada     = 0.0;
  elem_no  = 0;

	
  /* -------------- Cells indexes and positions --------------------*/
  fprintf(fMOut, "%%\n");

  fprintf (fMOut, "%%--------------------------------------------LOST------------------------------------------------\n");
  fprintf(fMOut,"%% %5s  %12s %12s %12s %12s %12s %12s   %12s %12s %12s %12s\n", "Elem",
	  "XMin", "XMax", "YMin", "YMax", "ZMin", "ZMax", "Total", "Ionized", "Reflected","Time out");
  /* for by Element */
  for(elem_no=0;elem_no<elemente;elem_no++){
	  dLostT[elem_no] = (looser->x_min[elem_no] + looser->x_max[elem_no] + looser->y_min[elem_no] + 
		looser->y_max[elem_no] + looser->z_min[elem_no] + looser->z_max[elem_no]);
	fprintf(fMOut,"%% %5s  %12.2e %12.2e %12.2e %12.2e %12.2e %12.2e   %12.2e %12.2e %12i %12.2e\n", 
		looser->name[elem_no], looser->x_min[elem_no], looser->x_max[elem_no], 
		looser->y_min[elem_no], looser->y_max[elem_no], looser->z_min[elem_no], looser->z_max[elem_no],
		dLostT[elem_no], looser->ionized[elem_no],looser->refl_cnt[elem_no], looser->time_out[elem_no]);
  }
#ifdef WPLATE
    y_snet = &(netz[0][0]); /* H.E.S. is implied . . .*/ 
	fprintf(fMOut,"%%   Ch. Elem  %12s %12s %12s %12s %12s %12s \n",  
		          "PlateIn", "PlateInS", "PlateInT", "PlateOut", "ShieldIn", "ShieldOut");
	for (iCL = 0; iCL<MAXBeCH; iCL++){
	  fprintf(fMOut,"%% %2i <%s> %12.2e %12.2e %12.2e %12.2e %12.2e %12.2e \n",  iCL,  y_snet->content[0].element,
	            looser->WP.PlateIn[iCL], 
				looser->WP.PlateInS[iCL], looser->WP.PlateInT[iCL],
				looser->WP.PlateOut[iCL], 
				looser->WP.ShieldIn[iCL], looser->WP.ShieldOut[iCL]);
	}
	fprintf(fMOut,"%% Total WPLATE[%i]: %12.2e \n" ,elem_no, looser->WPlate[0]);
    fprintf(fMOut,"WPlateD = [ \n");
    for (iCL = 0; iCL<MAXBeCH; iCL++)
	fprintf(fMOut," %12.2e         %12.2e %12.2e %12.2e %12.2e %12.2e %12.2e  \n",  
	            looser->WPlate[elem_no],
	            looser->WP.PlateIn[iCL], 
				looser->WP.PlateInS[iCL], looser->WP.PlateInT[iCL],
				looser->WP.PlateOut[iCL], 
				looser->WP.ShieldIn[iCL], looser->WP.ShieldOut[iCL]);
    fprintf(fMOut,"]; \n");
#endif
  fprintf (fMOut, "%%--------------------------------------------LOST------------------------------------------------\n");
  fprintf (fMOut, "%%\n");  


/* -------------- Cells indexes and positions --------------------*/
  fprintf (fMOut, "\n");    
  fprintf (fMOut, "%%-------------------------------Element independent cell properties------------------------------\n");
  fprintf (fMOut, "%%\n");  

  MdthX = (double*)calloc(iNCx*iNCy,sizeof(double));
  MdthY = (double*)calloc(iNCx*iNCy,sizeof(double));	
#ifdef UGEOMETRY
  MdthZ = (double*)calloc(iNCx*iNCy,sizeof(double));
  nVertices = (int*)calloc(iNCx*iNCy,sizeof(int));
  xVert = (double**)calloc(5,sizeof(double));	// reserved up to 5 vertices
  yVert = (double**)calloc(5,sizeof(double));	// reserved up to 5 vertices
  zVert = (double**)calloc(5,sizeof(double));	// reserved up to 5 vertices
  for(iv = 0; iv<5; iv++)
  { 
	  xVert[iv] = (double*)calloc(iNCx*iNCy,sizeof(double));
	  yVert[iv] = (double*)calloc(iNCx*iNCy,sizeof(double));
	  zVert[iv] = (double*)calloc(iNCx*iNCy,sizeof(double));	
  }
#endif
  iClX = (int*)calloc(iNCx*iNCy,sizeof(int));	
  iClY = (int*)calloc(iNCx*iNCy,sizeof(int));	
  Area = (double*)calloc(iNCx*iNCy,sizeof(double));	
  Cellsize = (double*)calloc(iNCx*iNCy,sizeof(double));
  SurfT = (double*)calloc(iNCx*iNCy,sizeof(double));	
  ShPot = (double*)calloc(iNCx*iNCy,sizeof(double));   
  pTe = (double*)calloc(iNCx*iNCy,sizeof(double));
  pTi = (double*)calloc(iNCx*iNCy,sizeof(double));
  pNe = (double*)calloc(iNCx*iNCy,sizeof(double));	
  H_flc = (double*)calloc(iNCx*iNCy,sizeof(double));

#ifdef BE_CARBYDE
  pdCrb = (double*)calloc(iNCx*iNCy,sizeof(double));
#endif
#ifdef ITER_BM
  Shadow = (int*)calloc(iNCx*iNCy,sizeof(int));
#ifdef JETLim
  CLen = (double*)calloc(iNCx*iNCy,sizeof(double));
#endif
#endif


  for(x_cell=0;x_cell<iNCx;x_cell++){	 
    for(y_cell=0;y_cell<iNCy;y_cell++){
      y_snet = &(netz[x_cell][y_cell]);/*die richtige x-reihe des oberflachennetzes einsetzen*/
      MdthX[x_cell + y_cell*iNCx] =	y_snet->midth[0];
      MdthY[x_cell + y_cell*iNCx] =	y_snet->midth[1];
	#ifdef UGEOMETRY
	  MdthZ[x_cell + y_cell*iNCx] =	y_snet->midth[2];
	  nVertices[x_cell + y_cell*iNCx] =	y_snet->nVertices;
	  for (iv = 0; iv < 5; iv++)
	  {
		  if (iv < y_snet->nVertices)
		  {
			xVert[iv][x_cell + y_cell*iNCx] =	y_snet->xvertex[iv];
			yVert[iv][x_cell + y_cell*iNCx] =	y_snet->yvertex[iv];
			zVert[iv][x_cell + y_cell*iNCx] =	y_snet->zvertex[iv];
		  }
		  else
		  {
			xVert[iv][x_cell + y_cell*iNCx] =	0.0;
			yVert[iv][x_cell + y_cell*iNCx] =	0.0;
			zVert[iv][x_cell + y_cell*iNCx] =	0.0;
		  }
	  }
    #endif
	  iClX[x_cell + y_cell*iNCx] =	x_cell + 1;
      iClY[x_cell + y_cell*iNCx] =	y_cell + 1;        
      /*SDcmt: local surface temperature*/
      SurfT[x_cell + y_cell*iNCx]=y_snet->AT_surf;
      /*SDcmt: cell size*/
      Cellsize[x_cell + y_cell*iNCx]=y_snet->cellsize;								
      /*SDcmt: area of the cell [mm^2](AREA)*/										
      if(y_snet->area !=0.){ 
	Area[x_cell + y_cell*iNCx]=y_snet->area;								
      }else
	Area[x_cell + y_cell*iNCx]=nada;

	  ShPot[x_cell + y_cell*iNCx]=y_snet->dBiasShPot;
	  pTe[x_cell + y_cell*iNCx]=y_snet->Te;
	  pTi[x_cell + y_cell*iNCx]=y_snet->Ti;
	  pNe[x_cell + y_cell*iNCx]=y_snet->ne;
	  H_flc[x_cell + y_cell*iNCx]=y_snet->H_fluence/y_snet->area;
#ifdef BE_CARBYDE
	  pdCrb[x_cell + y_cell*iNCx]=y_snet->BeCarb[BeC_All];
#endif
#ifdef ITER_BM
      Shadow[x_cell + y_cell*iNCx] = y_snet->iShadow;
#ifdef JETLim
      CLen[x_cell + y_cell*iNCx] = y_snet->dCL;
#endif
#endif
    }/*y_cell*/   
  }/*x_cell*/

  iCP = 0;

  iCP++;
  fprintf(fMOut,"CellProp{%i}.name='dMdlX';\n",iCP);
  sprintf(sVarN,"CellProp{%i}.data",iCP);
  MlMatrOut (fMOut, sVarN, MdthX, iNCx, iNCy, 'd');
  fprintf(fMOut,"dMdlX=CellProp{%i}.data;\n\n",iCP);

  iCP++;
  fprintf(fMOut,"CellProp{%i}.name='dMdlY';\n",iCP);
  sprintf(sVarN,"CellProp{%i}.data",iCP);
  MlMatrOut (fMOut, sVarN, MdthY, iNCx, iNCy, 'd');
  fprintf(fMOut,"dMdlY=CellProp{%i}.data;\n\n",iCP);

#ifdef UGEOMETRY
  iCP++;
  fprintf(fMOut,"CellProp{%i}.name='dMdlZ';\n",iCP);
  sprintf(sVarN,"CellProp{%i}.data",iCP);
  MlMatrOut (fMOut, sVarN, MdthZ, iNCx, iNCy, 'd');
  fprintf(fMOut,"dMdlZ=CellProp{%i}.data;\n\n",iCP);

  iCP++;
  fprintf(fMOut,"CellProp{%i}.name='iNVert';\n",iCP);
  sprintf(sVarN,"CellProp{%i}.data",iCP);
  MlMatrOut (fMOut, sVarN, nVertices, iNCx, iNCy, 'i');
  fprintf(fMOut,"iNVert=CellProp{%i}.data;\n\n",iCP);

  for (iv = 0; iv < 5; iv++)
  {
	  iCP++;
	  fprintf(fMOut,"CellProp{%i}.name='dXVert%i';\n",iCP,iv);
	  sprintf(sVarN,"CellProp{%i}.data",iCP);
	  MlMatrOut (fMOut, sVarN, xVert[iv], iNCx, iNCy, 'd');
	  fprintf(fMOut,"dXVert%i=CellProp{%i}.data;\n\n",iv,iCP);

	  iCP++;
	  fprintf(fMOut,"CellProp{%i}.name='dYVert%i';\n",iCP,iv);
	  sprintf(sVarN,"CellProp{%i}.data",iCP);
	  MlMatrOut (fMOut, sVarN, yVert[iv], iNCx, iNCy, 'd');
	  fprintf(fMOut,"dYVert%i=CellProp{%i}.data;\n\n",iv,iCP);

	  iCP++;
	  fprintf(fMOut,"CellProp{%i}.name='dZVert%i';\n",iCP,iv);
	  sprintf(sVarN,"CellProp{%i}.data",iCP);
	  MlMatrOut (fMOut, sVarN, zVert[iv], iNCx, iNCy, 'd');
	  fprintf(fMOut,"dZVert%i=CellProp{%i}.data;\n\n",iv,iCP);
  }

#endif

  iCP++;
  fprintf(fMOut,"CellProp{%i}.name='iClX';\n",iCP);
  sprintf(sVarN,"CellProp{%i}.data",iCP);
  MlMatrOut (fMOut, sVarN, iClX, iNCx, iNCy, 'i');
  fprintf(fMOut,"iClX=CellProp{%i}.data;\n\n",iCP);

  iCP++;
  fprintf(fMOut,"CellProp{%i}.name='iClY';\n",iCP);
  sprintf(sVarN,"CellProp{%i}.data",iCP);
  MlMatrOut (fMOut, sVarN, iClY, iNCx, iNCy, 'i');
  fprintf(fMOut,"iClY=CellProp{%i}.data;\n\n",iCP);
  
  iCP++;
  fprintf(fMOut,"CellProp{%i}.name='dSurfT';\n",iCP);
  sprintf(sVarN,"CellProp{%i}.data",iCP);
  MlMatrOut (fMOut, sVarN, SurfT, iNCx, iNCy, 'd');

  iCP++;
  fprintf(fMOut,"CellProp{%i}.name='dArea';\n",iCP);
  sprintf(sVarN,"CellProp{%i}.data",iCP);
  MlMatrOut (fMOut, sVarN, Area, iNCx, iNCy, 'd');

  iCP++;
  fprintf(fMOut,"CellProp{%i}.name='dCellSize';\n",iCP);
  sprintf(sVarN,"CellProp{%i}.data",iCP);
  MlMatrOut (fMOut, sVarN, Cellsize, iNCx, iNCy, 'd');

  iCP++;
  fprintf(fMOut,"CellProp{%i}.name='dSheathPot';\n",iCP);
  sprintf(sVarN,"CellProp{%i}.data",iCP);
  MlMatrOut (fMOut, sVarN, ShPot, iNCx, iNCy, 'd');

  iCP++;
  fprintf(fMOut,"CellProp{%i}.name='dH_fluence';\n", iCP);
  sprintf(sVarN,"CellProp{%i}.data",iCP);
  MlMatrOut (fMOut, sVarN, H_flc, iNCx, iNCy, 'd');

  iCP++;
  fprintf(fMOut,"CellProp{%i}.name='Te';\n",iCP);
  sprintf(sVarN,"CellProp{%i}.data",iCP);
  MlMatrOut (fMOut, sVarN, pTe, iNCx, iNCy, 'd');

  iCP++;
  fprintf(fMOut,"CellProp{%i}.name='Ti';\n",iCP);
  sprintf(sVarN,"CellProp{%i}.data",iCP);
  MlMatrOut (fMOut, sVarN, pTi, iNCx, iNCy, 'd');

  iCP++;
  fprintf(fMOut,"CellProp{%i}.name='ne';\n",iCP);
  sprintf(sVarN,"CellProp{%i}.data",iCP);
  MlMatrOut (fMOut, sVarN, pNe, iNCx, iNCy, 'd');

#ifdef BE_CARBYDE
  iCP++;
  fprintf(fMOut,"CellProp{%i}.name='Be2C';\n",iCP);
  sprintf(sVarN,"CellProp{%i}.data",iCP);
  MlMatrOut (fMOut, sVarN, pdCrb, iNCx, iNCy, 'd');
#endif

#ifdef ITER_BM
  iCP++;
  fprintf(fMOut,"CellProp{%i}.name='Shadow';\n",iCP);
  sprintf(sVarN,"CellProp{%i}.data",iCP);
  MlMatrOut (fMOut, sVarN, Shadow, iNCx, iNCy, 'i');
#ifdef JETLim
  iCP++;
  fprintf(fMOut,"CellProp{%i}.name='CLen';\n",iCP);
  sprintf(sVarN,"CellProp{%i}.data",iCP);
  MlMatrOut (fMOut, sVarN, CLen, iNCx, iNCy, 'd');
#endif
#endif

  free(iClX);
  free(iClY);
  free(MdthX);
  free(MdthY);
#ifdef UGEOMETRY
  free(MdthZ);
  free(nVertices);
  free(xVert);
  free(yVert);
  free(zVert);
#endif
  free(Area);
  free(Cellsize);
  free(SurfT);
  free(ShPot);
  free(pTe);
  free(pTi);
  free(pNe);
  free(H_flc);
#ifdef BE_CARBYDE
  free(pdCrb);
#endif
#ifdef ITER_BM
  free(Shadow);
#ifdef JETLim
  free(CLen);
#endif
#endif
  /*SDcmt: All elements:*/
  for(elem_no=0;elem_no<elemente;elem_no++)
    {	       
      y_snet = &(netz[0][0]);
       
      strncpy (sElementName, y_snet->content[elem_no].element, strlen (y_snet->content[elem_no].element));		
      fprintf (fMOut, "%%\n%% -------------------- Element '%s' -------------------\n",sElementName);
      fprintf (fMOut, "ElName(%i,:) = '%s';\n%%\n",elem_no+1,sElementName);
      
	  strncpy(CijSpez.name,y_snet->content[elem_no].element,2);
	  CijSpez.name[2]='\0';
	  ElemNo(&CijSpez);


      for (iStatusC=0;iStatusC<PARTSTATNO;iStatusC++)
	{			
	   
	   
	  fprintf (fMOut, "%%\n%% -------------------- '%.7s' -------------------\n",sParticleStatus[iStatusC]);
	  fprintf (fMOut, "StatusName(%i,:) = '%.7s';\n%%\n",iStatusC+1,sParticleStatus[iStatusC]);
	   
	   
	  /*-----------------------------------------------------
	    zweifache for()- Schleife:
	    1) * 2)  ==alle SurfaceCells 
	    # in: s_rand[1,0]
	    ------------------------------------------------------*/
	  fprintf(fMOut,"%%!Transposed Array!\n");
	  fprintf(fMOut, "SurfCell(%i,%i,:,:) = [\n",elem_no+1,iStatusC+1);
			

	  for(x_cell=0;x_cell<iNCx;x_cell++)
	    {				
	      for(y_cell=0;y_cell<iNCy;y_cell++)
		{
		  y_snet = &(netz[x_cell][y_cell]);

#ifndef HardSpSeq
          /* Searching for an element compartmnent  . . . */ 
		  mem_ctrl=FALSE;
				for(k=0;k<rand[2];k++){
					if(y_snet->which[y_snet->upmost][k] == CijSpez.elem){
						mem_ctrl = TRUE;
						merke = k;
						break;
					}
				}
#else /* HardSpSeq */
				merke = elem_no;
				mem_ctrl = TRUE;
#endif

		  switch (iStatusC)
		    {
		    case 0: /*SDcmt: abtransportierte Teilchen */							
		      fprintf(fMOut," %.4e",y_snet->content[elem_no].newgo/y_snet->area);
		      break;
		    case 1: /*SDcmt: neugelandete Teilchen */
		      fprintf(fMOut," %.4e",y_snet->content[elem_no].newcome/y_snet->area);
		      break;
		    case 2: /*SDcmt: in allen bisherigen!! Schritten transportierte Teilchen */
		      fprintf(fMOut," %.4e",y_snet->content[elem_no].transport/y_snet->area);
		      break;
		    case 3:	/*SDcmt: verbliebene/entstandene Resource */
		      fprintf(fMOut," %.4e",y_snet->content[elem_no].resource/y_snet->area);
		      break;
		    case 4: /*SDcmt: mittlere Energie der gelandeten Teilchen*/
		      if(y_snet->content[elem_no].newcome!=0.)
			{ 								
			  fprintf(fMOut," %.4e",
				  y_snet->content[elem_no].energy/y_snet->content[elem_no].newcome);								
			}
		      else
			fprintf(fMOut," %.4e",nada);							
		      break;
		    case 5:	/*SDcmt: mittlere Ladung der gelandeten Teilchen*/							
		      if(y_snet->content[elem_no].newcome!=0.)
			{ 
			  fprintf(fMOut," %.4e",
				  y_snet->content[elem_no].charge/y_snet->content[elem_no].newcome);								
			}
		      else
			fprintf(fMOut," %.4e",nada);
		      break;
		    case 6: /*SDcmt: in diesem Schritt abtransportierte Teilchen (SelbstZerstaeubung)*/					
		      fprintf(fMOut," %.4e",y_snet->content[elem_no].ngself/y_snet->area);
		      break;
		    case 7: /*SDcmt: Konzentration zu Beginn dieses Schrittes*/		
				if (mem_ctrl==FALSE){
					fprintf(fMOut," %.4e",0.);
				}else
					fprintf(fMOut," %.4e",y_snet->cij[merke]);
				break;
		    case 8:	/*SDcmt: mittlere Ladung der prompt gelandeten Teilchen*/						
		      if(y_snet->content[elem_no].prompt!=0.)
			{ 
			  fprintf(fMOut," %.4e",
				  y_snet->content[elem_no].prochar/y_snet->content[elem_no].prompt);								
			}
		      else
			fprintf(fMOut," %.4e",nada);
		      break;
		    case 9:  /*SDcmt: in diesem Schritt prompt deponierte Teilchen*/
		      fprintf(fMOut," %.4e",y_snet->content[elem_no].prompt/y_snet->area);
		      break;
		    case 10: /*SDcmt: in diesem Schritt prompt deponierte Teilchen (chemisch)*/
		      fprintf(fMOut," %.4e",y_snet->content[elem_no].prompt_chem/y_snet->area);
		      break;
		    case 11: /*SDcmt: in diesem Schritt aus background deponierte C-Teilchen*/										
		      fprintf(fMOut," %.4e",y_snet->content[elem_no].c_back/y_snet->area);
		      break;
		    case 12: /*BCmt: shortage of particle in the layer) */										
		      fprintf(fMOut," %.4e",y_snet->dAtShortage[elem_no]);																			
		      break;
		    case 13:
		      fprintf(fMOut," %.4e",y_snet->content[elem_no].fluence/y_snet->area);
		      break;
		    case 14: /* Amount of chemically eroded particles . . . */		
				if (mem_ctrl==FALSE){
					fprintf(fMOut," %.4e",0.);
				}else
					fprintf(fMOut," %.4e",y_snet->thermic[merke]/y_snet->area);
				break;
		    case 15: /* Amount of physically eroded particles . . . */		
				if (mem_ctrl==FALSE){
					fprintf(fMOut," %.4e",0.);
				}else
					fprintf(fMOut," %.4e",y_snet->energetic[merke]/y_snet->area);
				break;
		    case 16: /* Amount of physically eroded particles . . . */		
				if (mem_ctrl==FALSE){
					fprintf(fMOut," %.4e",0.);
				}else
					fprintf(fMOut," %.4e",y_snet->content[elem_no].dPEr_Plasma/y_snet->area);
				break;
			case 17: /* B-field angle with surfface normal */
				/* BCmt: absolete PER6 Amount of physically eroded particles . . . */		
				if (mem_ctrl==FALSE){
					fprintf(fMOut," %.4e",0.);
				}else
					/* fprintf(fMOut," %.4e",y_snet->newenergetic[merke]/y_snet->area); BCmt: 'PER6' - absolete */
					fprintf(fMOut," %.4e",y_snet->b_angle);  
				break;
				
		    } /*SDcmt: Endof switch iStatusC */										
				
		}/*x_cell*/
	      fprintf(fMOut,"\n");
			
	    }/*y_cell*/
	  fprintf(fMOut, "]';\n%%\n");
	} /* SDcmt: for (iStatusC=0;iStatusC<iPartStatNo;iStatus++)*/
    } /* SDcmt endof for(elem_no=0;elem_no<elemente;elem_no++)   */							  
  fclose (fMOut);


/*--------------------------------------------------------------------------
        STARTOF OUTPUT of concentration file "xxxxfilm_st??.m"
--------------------------------------------------------------------------*/
  f_FILMTMP= MFOpen ("Filmdata", pfades->RESULT, step);
  /*------------------------------------------------------------------
  zweifache for()- Schleife ersetzt durch Makro
  -------------------------------------------------------------------*/
   
#ifdef HardSpSeq
  fprintf (f_FILMTMP, "%% ==== Hard element sequence ====\n");
  fprintf (f_FILMTMP, "iNTrEl = %i;\n", iActNElem);
  piHES = aiGetHES();
  MlMatrOut (f_FILMTMP, "HES_Ind", piHES, iActNElem, 1, 'i');
  fprintf (f_FILMTMP, "HES_Nm = [");
  for (iCEl =0; iCEl<iActNElem; iCEl++) 
	  fprintf (f_FILMTMP, " '%s' ", sHSeqNm(iCEl));
  fprintf (f_FILMTMP, "];\n");
  fprintf (f_FILMTMP, "%% ===============================\n%%\n");
#endif /* HardSpSeq */

  SURF_START(rand[0],rand[1],y_snet,netz);
  /*-------------------------------------------------------------------------
  writing the y_snet->LAYERs & y_snet->WHICHes
  --------------------------------------------------------------------------*/
  fprintf (f_FILMTMP, "%%\n%%------------------------------------------------------\n"); 
  fprintf (f_FILMTMP, "%% Cell (x=%i, y=%i), MidPoint XYZ(%.2e, %.2e, %.2e)\n",
      x_cell, y_cell, y_snet->midth[0],y_snet->midth[1], y_snet->midth[2]);
  sprintf (sVarN, "Layer(%i,%i,:,:)",x_cell+1, y_cell+1);
#ifndef HardSpSeq
  MlMatrOut (f_FILMTMP, sVarN, y_snet->layer, (MAX_N_ELEM), 50, 'd');
  sprintf (sVarN, "Which(%i,%i,:,:)",x_cell+1, y_cell+1);
  MlMatrOut (f_FILMTMP, sVarN, y_snet->which, (MAX_N_ELEM), 50, 'i');
#else
  dBuf = (double *) calloc (50*iActNElem, sizeof(double));
  for(iCL=0; iCL<50; iCL++)
   for(iCEl=0; iCEl<iActNElem; iCEl++)
     dBuf[iCEl+iActNElem*iCL] = y_snet->layer[iCL][iCEl];
  MlMatrOut (f_FILMTMP, sVarN, dBuf, iActNElem, 50, 'd');
  free(dBuf);
#endif

  SURF_ENDE ;

#ifdef BE_CARBYDE
  dBufBC = (double *) calloc (rand[0]*rand[1], sizeof(double));
  SURF_START(rand[0],rand[1],y_snet,netz);
     dBufBC[counter++]=y_snet->BeCarb[BeC_All];
  SURF_ENDE;
  MlMatrOut (f_FILMTMP, "Be2C", dBufBC, rand[0], rand[1], 'd');
  free(dBufBC);
#endif

  fclose (f_FILMTMP);
   
/*--------------------------------------------------------------------------
        ENDOF OUTPUT of concentration file "xxxxfilm_st??.m"
--------------------------------------------------------------------------*/

				
}/*end of the routine!!!*/

#ifdef TRIM
/*--SDcmt----------------------------------------------------------------------------------------*/
/*                                                                                               */
/* Write the available surface information from SDTRimSP into one matlab output file                           */
/*                                                                                               */
/*--SDcmt----------------------------------------------------------------------------------------*/
void  matlab_trim_out( struct SurfaceCell      **netz,  /*heisst im move_ctrl() "net"                */
		       struct Path            *pfades,  /*heisst im move_ctrl() "pfade"              */
		       struct TimeCell         *timer,  /*heisst im move_ctrl() "time_ptr"           */
		       struct LostParticle    *looser,  /*heisst im move_ctrl() "lost"               */                 
		       int                      *rand,  /*heisst im move_ctrl() "grid_rand"          */
		       double                  *bound,  /*heisst im move_ctrl() "grid_bound"         */
		       int                       step,  /*heisst im move_ctrl() "time_g"             */
		       int                       surf,  /* Surface type */
		       int                 iMain_seed,
			   struct SurfaceGridData  *sgrid			   
		       )
{
  /* ----------------------------------------------------------------------
     Variablendefinition
     -------------------------------------------------------------------------*/
  /* beginning of each file             */

  int                elem_no,slab_no;/* laufvariablen                      */
  int                bgtr;       /* laufvariablen                      */
  int                x_cell;         /* laufvariablen                      */
  int                y_cell;         /* laufvariablen                      */

  struct SurfaceCell *y_snet;        /* zugehoerige laufpointer im sputter */
  /* zyklus                             */	

  double             nada;           /* hilfsvariable bei ausgabe          */	

  char sFN_NTrv[N_OUTPUT_FILES][DEF_STR_LEN];		   /* Filename part */   
  char sElementName[] ="____";	   /* Fixed length Element name */

  /* Number of different trim particle "statuses": */	
#define TRIMPARTSTATNO 15
  char sParticleStatus[TRIMPARTSTATNO][15];/* Explanation of the particle status*/	
  int iStatusC,iFiles;					    /* Particle status counter, filecounter */
  char sBgtr[2][3];

  int iNCx, iNCy;
  double *MdthX, *MdthY, *DeltaDepth, *TotDeltaDepth, *SurfT, *Cellsize,*Area, *ShPot;
  int *iClX, *iClY;
  char sVarN[DEF_STR_LEN]; 

  double *qux_buf;
  
  int distrib_no;

  FILE *fMout[N_OUTPUT_FILES];
#ifdef BE_CARBYDE
  double Be2C_concentration[3][MAX_N_SLABS];
#endif

  /* ----------------------------------------------------------------------
     endof Variablendefinition
     -------------------------------------------------------------------------*/
  iNCx=rand[1]; iNCy=rand[0];

  /*SDcmt: Initialize array of status strings:
   *       Changing s.th. with the particle status means, that s.th. 
   *       might need to be changed in matlab_trim_in and in MERO.
   */
  strcpy(sParticleStatus[0] ,"fluence       ");
  strcpy(sParticleStatus[1] ,"backscattered ");
  strcpy(sParticleStatus[SC_SPUTTERED] ,"sputtered     "); /*see matlab_trim_in*/
  strcpy(sParticleStatus[3] ,"reemitted     ");
  strcpy(sParticleStatus[4] ,"deposited     ");
  strcpy(sParticleStatus[5] ,"transmitted   ");
  strcpy(sParticleStatus[6] ,"stopped recoil");
  strcpy(sParticleStatus[7] ,"resource      ");
  strcpy(sParticleStatus[8] ,"chem sputtered");
  strcpy(sParticleStatus[SC_TOTAL_RESOURCE] ,"total resource");
  strcpy(sParticleStatus[10],"chem reflected");
  strcpy(sParticleStatus[SC_DELETED_COMP],"deleted compon");
  strcpy(sParticleStatus[12],"stopped proj  ");
  strcpy(sParticleStatus[13],"av penetration");
  strcpy(sParticleStatus[14],"mx penetration");
  strcpy(sBgtr[0],"TOT"); /*from Background*/
  strcpy(sBgtr[1],"FO"); /*from Followed particles*/
		
  /*SDcmt: Initialize buffer for qux*/
  qux_buf = (double*)calloc(sgrid->IN_ti.ncp*sgrid->IN_ti.nqx,sizeof(double));


  /*SDcmt: Open File for writing (w+)*/
  sprintf (sFN_NTrv[0],"Trim_step%i",step);
  sprintf (sFN_NTrv[1],"Trim3D_step%i",step);
  sprintf (sFN_NTrv[2],"TrimDistr_step%i",step);

  for(iFiles=0;iFiles<N_OUTPUT_FILES;iFiles++){
	fMout[iFiles]= MFOpen (sFN_NTrv[iFiles], pfades->RESULT, -1);
	
  /*SDcmt: Write infos valid for all elements (& files)*/
  fprintf(fMout[iFiles], "%%\n");
  fprintf(fMout[iFiles], "%% Number of elements:\n");
  fprintf(fMout[iFiles], "iNElem=%i;\n",sgrid->IN_ti.ncp);
  fprintf(fMout[iFiles], "%%\n");
  fprintf(fMout[iFiles], "%% Number of cell features\n");
  switch (iFiles){
			   case 0:
				   fprintf(fMout[iFiles], "iNClFtr=%i;\n",2*TRIMPARTSTATNO);
				   break;
			   case 1:
				   fprintf(fMout[iFiles], "iNClFtr=%i;\n",0);
				   break;
			   case 2:fprintf(fMout[iFiles], "iNClFtr=%i;\n",0);
				   break;	
			   default:
				   fprintf(fMout[iFiles], "iNClFtr=%i;\n",-1); /*ERROR?!*/
  }

  fprintf(fMout[iFiles], "%%\n");
  fprintf(fMout[iFiles], "%% -----------------------------------------------------:\n");
  fprintf(fMout[iFiles], "%% Informations valid for all elements and cells:\n");
  fprintf(fMout[iFiles], "iSType=%d;\n",surf);
  fprintf(fMout[iFiles], "iTimeStepNo=%d;\n",step);
  fprintf(fMout[iFiles], "dStepTime=%f;\n",timer->step); 
  fprintf(fMout[iFiles], "dRealTime=%f;\n",timer->dRealTm); 
  fprintf(fMout[iFiles], "dTotalTime=%f;\n",timer->total);
	
  fprintf(fMout[iFiles], "dXmin= %.1f;\n",bound[0]);
  fprintf(fMout[iFiles], "dXmax= %.1f;\n",bound[1]);
  fprintf(fMout[iFiles], "dYmin= %.1f;\n",bound[2]);
  fprintf(fMout[iFiles], "dYmax= %.1f;\n",bound[3]);
  fprintf(fMout[iFiles], "dDeltaX=%.1f;\n",bound[4]);
  fprintf(fMout[iFiles], "dDeltaY=%.1f;\n",bound[5]);
  fprintf(fMout[iFiles], "iNX=%d;\n",iNCx);
  fprintf(fMout[iFiles], "iNY=%d;\n",iNCy);

  fprintf(fMout[iFiles], "Random_Seed=%d;\n",iMain_seed);

  fprintf(fMout[iFiles], "idrel=%i;\n",sgrid->IN_ti.idrel); 
  fprintf(fMout[iFiles], "isbv=%i;\n",sgrid->IN_ti.isbv); 
  fprintf(fMout[iFiles], "nqx=%i;\n",sgrid->IN_ti.nqx); 
  fprintf(fMout[iFiles], "ttarget=%f;\n",sgrid->IN_ti.ttarget); 

  fprintf(fMout[iFiles], "%% -----------------------------------------------------:\n");
  fprintf(fMout[iFiles], "%% Informations of the elements:\n");

  for (elem_no=0;elem_no<sgrid->IN_ti.ncp;elem_no++){	
    strncpy (sElementName, sgrid->speciesTRIM[elem_no].name, strlen (sgrid->speciesTRIM[elem_no].name));	
    fprintf (fMout[iFiles], "%%\n%% -------------------- Element '%s' -------------------\n",sElementName);
    fprintf (fMout[iFiles], "ElName(%i,:) = '%s';\n%%\n",elem_no+1,sElementName);
    fprintf (fMout[iFiles], "qubeam(%i) = %e;\n%%\n"  ,elem_no+1, sgrid->IN_PART_cond[elem_no].qubeam);
    fprintf (fMout[iFiles], "qu(%i) = %e;\n%%\n"      ,elem_no+1, sgrid->IN_PART_cond[elem_no].qu);
    fprintf (fMout[iFiles], "e_cutoff(%i) = %e;\n%%\n",elem_no+1, sgrid->IN_PART_cond[elem_no].e_cutoff);
  }
  }

  nada     = 0.0;
  elem_no  = 0;

	
  /* -------------- Cells indexes and positions --------------------*/
  fprintf(fMout[0], "%%\n");
  MdthX = (double*)calloc(iNCx*iNCy,sizeof(double));
  MdthY = (double*)calloc(iNCx*iNCy,sizeof(double));	
  iClX = (int*)calloc(iNCx*iNCy,sizeof(int));	
  iClY = (int*)calloc(iNCx*iNCy,sizeof(int));	
  DeltaDepth = (double*)calloc(iNCx*iNCy,sizeof(double));	
  TotDeltaDepth = (double*)calloc(iNCx*iNCy,sizeof(double));
  Area = (double*)calloc(iNCx*iNCy,sizeof(double));	
  Cellsize = (double*)calloc(iNCx*iNCy,sizeof(double));
  SurfT = (double*)calloc(iNCx*iNCy,sizeof(double));	
  ShPot = (double*)calloc(iNCx*iNCy,sizeof(double));	

  for(x_cell=0;x_cell<iNCx;x_cell++){	 
    for(y_cell=0;y_cell<iNCy;y_cell++){
      y_snet = &(netz[x_cell][y_cell]);/*die richtige x-reihe des oberflachennetzes einsetzen*/
      MdthX[x_cell + y_cell*iNCx] =	y_snet->midth[0];
      MdthY[x_cell + y_cell*iNCx] =	y_snet->midth[1];
      iClX[x_cell + y_cell*iNCx] =	x_cell + 1;
      iClY[x_cell + y_cell*iNCx] =	y_cell + 1;        
      TotDeltaDepth[x_cell + y_cell*iNCx]=y_snet->dDeltaDepth+y_snet->dTotDeltaDepth;
      DeltaDepth[x_cell + y_cell*iNCx]=y_snet->dDeltaDepth;
      /*SDcmt: local surface temperature*/
      SurfT[x_cell + y_cell*iNCx]=y_snet->AT_surf;
      /*SDcmt: cell size*/
      Cellsize[x_cell + y_cell*iNCx]=y_snet->cellsize;								
      /*SDcmt: area of the cell [mm^2](AREA)*/										
	  if(y_snet->area !=0.){ 
		  Area[x_cell + y_cell*iNCx]=y_snet->area;								
	  }else
		  Area[x_cell + y_cell*iNCx]=nada;

	  ShPot[x_cell + y_cell*iNCx]=y_snet->dBiasShPot;
    }/*y_cell*/   
  }/*x_cell*/

  for(iFiles=0;iFiles<N_OUTPUT_FILES;iFiles++){
	  fprintf(fMout[iFiles],"CellProp{%i}.name='dMdlX';\n",CP_DMDLX);
	  sprintf(sVarN,"CellProp{%i}.data",CP_DMDLX);
	  MlMatrOut (fMout[iFiles], sVarN, MdthX, iNCx, iNCy, 'd');
	  fprintf(fMout[iFiles],"dMdlX=CellProp{%i}.data;\n\n",CP_DMDLX);

	  fprintf(fMout[iFiles],"CellProp{%i}.name='dMdlY';\n",CP_DMDLY);
	  sprintf(sVarN,"CellProp{%i}.data",CP_DMDLY);
	  MlMatrOut (fMout[iFiles], sVarN, MdthY, iNCx, iNCy, 'd');
	  fprintf(fMout[iFiles],"dMdlY=CellProp{%i}.data;\n\n",CP_DMDLY);

	  fprintf(fMout[iFiles],"CellProp{%i}.name='iClX';\n",CP_ICLX);
	  sprintf(sVarN,"CellProp{%i}.data",CP_ICLX);
	  MlMatrOut (fMout[iFiles], sVarN, iClX, iNCx, iNCy, 'i');
	  fprintf(fMout[iFiles],"iClX=CellProp{%i}.data;\n\n",CP_ICLX);

	  fprintf(fMout[iFiles],"CellProp{%i}.name='iClY';\n",CP_ICLY);
	  sprintf(sVarN,"CellProp{%i}.data",CP_ICLY);
	  MlMatrOut (fMout[iFiles], sVarN, iClY, iNCx, iNCy, 'i');
	  fprintf(fMout[iFiles],"iClY=CellProp{%i}.data;\n\n",CP_ICLY);
  }
  
  fprintf(fMout[0],"CellProp{%i}.name='dSurfT';\n",CP_DSURFT);
  sprintf(sVarN,"CellProp{%i}.data",CP_DSURFT);
  MlMatrOut (fMout[0], sVarN, SurfT, iNCx, iNCy, 'd');

  fprintf(fMout[0],"CellProp{%i}.name='dArea';\n",CP_DAREA);
  sprintf(sVarN,"CellProp{%i}.data",CP_DAREA);
  MlMatrOut (fMout[0], sVarN, Area, iNCx, iNCy, 'd');

  fprintf(fMout[0],"CellProp{%i}.name='dCellSize';\n",CP_DCELLSIZE);
  sprintf(sVarN,"CellProp{%i}.data",CP_DCELLSIZE);
  MlMatrOut (fMout[0], sVarN, Cellsize, iNCx, iNCy, 'd');

  fprintf(fMout[0],"CellProp{%i}.name='dSheathPot';\n",CP_DSHPOT);
  sprintf(sVarN,"CellProp{%i}.data",CP_DSHPOT);
  MlMatrOut (fMout[0], sVarN, ShPot, iNCx, iNCy, 'd');

  fprintf(fMout[0],"CellProp{%i}.name='dDeltaDepth';\n",CP_DDELTADEPTH);
  sprintf(sVarN,"CellProp{%i}.data",CP_DDELTADEPTH);
  MlMatrOut (fMout[0], sVarN, DeltaDepth, iNCx, iNCy, 'd');

  fprintf(fMout[0],"CellProp{%i}.name='dTotDeltaDepth';\n",CP_DTOTDELTADEPTH);
  sprintf(sVarN,"CellProp{%i}.data",CP_DTOTDELTADEPTH);
  MlMatrOut (fMout[0], sVarN, TotDeltaDepth, iNCx, iNCy, 'd');

 

  free(iClX);
  free(iClY);
  free(MdthX);
  free(MdthY);
  free(DeltaDepth);
  free(TotDeltaDepth);
  free(Area);
  free(Cellsize);
  free(SurfT);
  free(ShPot);

 fprintf(fMout[0], "%% -----------------------------------------------------:\n");
 fprintf(fMout[0], "%% Informations of the cells:\n");

  /*SDcmt: All elements:*/
 for(elem_no=0;elem_no<sgrid->IN_ti.ncp;elem_no++){     
      y_snet = &(netz[0][0]);
      strncpy (sElementName, sgrid->speciesTRIM[elem_no].name, strlen (sgrid->speciesTRIM[elem_no].name));		
      fprintf (fMout[0], "%%\n%% -------------------- Element '%s' -------------------\n",sElementName);
      for (bgtr=0;bgtr<2;bgtr++){
	for (iStatusC=0;iStatusC<TRIMPARTSTATNO;iStatusC++){				   	   
	  fprintf (fMout[0], "%%\n%% -------------------- '%.14s' -------------------\n",sParticleStatus[iStatusC]);
	  fprintf (fMout[0], "StatusName(%i,:) = '%.2s:%.14s';\n%%\n",bgtr*TRIMPARTSTATNO+(iStatusC+1),sBgtr[bgtr],sParticleStatus[iStatusC]);
	   	   
	  /*-----------------------------------------------------
	    zweifache for()- Schleife:
	    1) * 2)  ==alle SurfaceCells 
	    # in: s_rand[1,0]
	    ------------------------------------------------------*/
	  fprintf(fMout[0],"%%!Transposed Array!\n");
	  fprintf(fMout[0], "SurfCell(%i,%i,:,:) = [\n",elem_no+1,bgtr*TRIMPARTSTATNO+(iStatusC+1));
			
	  for(x_cell=0;x_cell<iNCx;x_cell++)
	    {				
	      for(y_cell=0;y_cell<iNCy;y_cell++)
		{
		  y_snet = &(netz[x_cell][y_cell]);
		  if (bgtr==0 && iStatusC==0 && sgrid->surface_model!=3 /*Once for every element and cell!*/){ /*Really, really total resource*/
		    y_snet->PART_quant[0][elem_no].dTotResource+=y_snet->PART_quant[1][elem_no].dTotResource;
		    y_snet->PART_quant[1][elem_no].dTotResource=0.;
		    if (sgrid->add_FO_to_BG && elem_no==0 /*Once for every cell!*/)
		      Add_PART_quant(&(y_snet->PART_quant[0][0]),&(y_snet->PART_quant[0][0]),&(y_snet->PART_quant[1][0]),sgrid->IN_ti.ncp); 
		  }
		  
		  switch (iStatusC)
		    {     
		    case 0: /*SDcmt:  fluence */							
		      fprintf(fMout[0]," %.4e",y_snet->PART_quant[bgtr][elem_no].dFlc/y_snet->area);
		      break;
		    case 1: /*SDcmt:  reflected */
		      fprintf(fMout[0]," %.4e",y_snet->PART_quant[bgtr][elem_no].dReflected/y_snet->area);
		      break;
		    case SC_SPUTTERED: /*SDcmt:  sputtered */
		      fprintf(fMout[0]," %.4e",y_snet->PART_quant[bgtr][elem_no].dSputtered/y_snet->area); 
		      break;
		    case 3: /*SDcmt:  re-emitted */
		      fprintf(fMout[0]," %.4e",y_snet->PART_quant[bgtr][elem_no].dReemitted/y_snet->area);
		      break;
		    case 4: /*SDcmt:  deposited */
		      fprintf(fMout[0]," %.4e",y_snet->PART_quant[bgtr][elem_no].dDeposited/y_snet->area);
		      break;
		    case 5: /*SDcmt:  transmitted */
		      fprintf(fMout[0]," %.4e",y_snet->PART_quant[bgtr][elem_no].dTransmitted/y_snet->area);
		      break;		  							
		    case 6: /*SDcmt:  stopped recoils */
		      fprintf(fMout[0]," %.4e",y_snet->PART_quant[bgtr][elem_no].dStoppedRecoil/y_snet->area);
		      break;
		    case 7: /*SDcmt:  resource */
		      fprintf(fMout[0]," %.4e",y_snet->PART_quant[bgtr][elem_no].dResource/y_snet->area);
		      break;		  							
		    case 8: /*SDcmt:  chemical sputtered */
		      fprintf(fMout[0]," %.4e",y_snet->PART_quant[bgtr][elem_no].dSputtered_chem/y_snet->area);
		      break;		  							
		    case SC_TOTAL_RESOURCE: /*SDcmt:  total resource*/
		      fprintf(fMout[0]," %.4e",y_snet->PART_quant[bgtr][elem_no].dTotResource/y_snet->area);
		      break;		  							
		    case 10: /*SDcmt: chemical reflected*/
		      fprintf(fMout[0]," %.4e",y_snet->PART_quant[bgtr][elem_no].dReflected_chem/y_snet->area);
		      break;		  		
			case SC_DELETED_COMP: /*SDcmt: from completely deleted components*/
				fprintf(fMout[0]," %.4e",y_snet->PART_quant[bgtr][elem_no].dDeletedComp/y_snet->area);
		      break;
		    case 12: /*SDcmt: stopped projectiles*/
		      fprintf(fMout[0]," %.4e",y_snet->PART_quant[bgtr][elem_no].dStoppedProj/y_snet->area);
		      break;
		    case 13: /*SDcmt: average penetration depth*/
		      fprintf(fMout[0]," %.4e",y_snet->PART_quant[bgtr][elem_no].dAvPenetration);
		      break;		  												
 		    case 14: /*SDcmt: average penetration depth*/
		      fprintf(fMout[0]," %.4e",y_snet->PART_quant[bgtr][elem_no].dMaxPenetration);
		      break;		  												

		    } /*SDcmt: Endof switch iStatusC */										
		  
		}/*x_cell*/
	      fprintf(fMout[0],"\n");
			
	    }/*y_cell*/
	  fprintf(fMout[0], "]';\n%%\n");
	} /* SDcmt: for (iStatusC=0;iStatusC<iPartStatNo;iStatus++)*/
      } /* SDcmt: for (bgtr=0...) */
    } /* SDcmt endof for(elem_no=0;elem_no<elemente;elem_no++)   */				
			  
 fprintf (fMout[1], "%%\n%% -----------------------------------------------------------------\n");
 fprintf (fMout[1], "%%\n%% ----------------- START OF DEPTH DEPENDENT DATA -----------------\n");
 fprintf(fMout[1],"%%Depth dependent data: array or matrix for each surface cell:\n");
 fprintf(fMout[1],"CellProp{%i}.name='qux'; %%atomic fraction\n",CP_QUX);
 fprintf(fMout[1],"CellProp{%i}.name='xno'; %%depth [angstrom]\n",CP_XNO);
 fprintf(fMout[1],"CellProp{%i}.name='dns'; %%average atomic density [atoms angstrom^(-3)]\n",CP_DNS);																	
 fprintf(fMout[1],"CellProp{%i}.name='BG:pdepth'; %%penetration depth of projectiles from background [atoms]\n",CP_PDEPTHBG);																	
 fprintf(fMout[1],"CellProp{%i}.name='FO:pdepth'; %%penetration depth of projectiles from followed particles [atoms]\n",CP_PDEPTHFO);																	
#ifdef BE_CARBYDE
 fprintf(fMout[1],"CellProp{%i}.name='Be2C'; %%concentration of Be2C \n",CP_BE2C_TRIM);
#endif
 
 for(x_cell=0;x_cell<iNCx;x_cell++){		
   for(y_cell=0;y_cell<iNCy;y_cell++){ 	
     y_snet = &(netz[x_cell][y_cell]);
     /* SDcmt: Write the qux buffer for MlMatrOut,*/
     /* Everything in mm^2*/
     y_snet = &(netz[x_cell][y_cell]);
     for(elem_no=0;elem_no<sgrid->IN_ti.ncp;elem_no++){
       for (slab_no=0;slab_no<sgrid->IN_ti.nqx;slab_no++){
	     qux_buf[slab_no+elem_no*sgrid->IN_ti.nqx]=y_snet->qux[elem_no][slab_no];
       }
     }

     fprintf (fMout[1], "%%\n%% -------------------- DEPTH DEPENDENT DATA OF CELL (%i,%i) -------------------\n",x_cell+1,y_cell+1);

     fprintf(fMout[1],"%% qux (atomic fraction)\n");
     sprintf(sVarN,"CellProp{%i}.data(%i,%i,:,:)",CP_QUX,x_cell+1,y_cell+1);
     MlMatrOut (fMout[1], sVarN, qux_buf, sgrid->IN_ti.nqx, sgrid->IN_ti.ncp, 'd');     
    
     fprintf(fMout[1],"%% xno (depth [angstrom])\n");
     sprintf(sVarN,"CellProp{%i}.data(%i,%i,:)",CP_XNO,x_cell+1,y_cell+1);
     MlMatrOut (fMout[1], sVarN, y_snet->xno, sgrid->IN_ti.nqx+1, 1, 'd');

     fprintf(fMout[1],"%% dns (average atomic density [atoms angstrom^(-3)])\n");
     sprintf(sVarN,"CellProp{%i}.data(%i,%i,:)",CP_DNS,x_cell+1,y_cell+1);
     MlMatrOut (fMout[1], sVarN, y_snet->dns, sgrid->IN_ti.nqx, 1, 'd');           

	 /* SDcmt: Write the qux buffer for MlMatrOut,*/
	 for(elem_no=0;elem_no<sgrid->IN_ti.ncp;elem_no++){
       for (slab_no=0;slab_no<sgrid->IN_ti.nqx;slab_no++){
		   qux_buf[slab_no+elem_no*sgrid->IN_ti.nqx]=y_snet->pdepthBG[elem_no][slab_no];
       }
     }
	 fprintf(fMout[1],"%% penetration depth background\n");
     sprintf(sVarN,"CellProp{%i}.data(%i,%i,:,:)",CP_PDEPTHBG,x_cell+1,y_cell+1);
     MlMatrOut (fMout[1], sVarN, qux_buf, sgrid->IN_ti.nqx, sgrid->IN_ti.ncp, 'd');     

	 /* SDcmt: Write the qux buffer for MlMatrOut,*/
	 for(elem_no=0;elem_no<sgrid->IN_ti.ncp;elem_no++){
       for (slab_no=0;slab_no<sgrid->IN_ti.nqx;slab_no++){
		   qux_buf[slab_no+elem_no*sgrid->IN_ti.nqx]=y_snet->pdepthFO[elem_no][slab_no];
       }
     }
	 fprintf(fMout[1],"%% penetration depth followed\n");
     sprintf(sVarN,"CellProp{%i}.data(%i,%i,:,:)",CP_PDEPTHFO,x_cell+1,y_cell+1);
     MlMatrOut (fMout[1], sVarN, qux_buf, sgrid->IN_ti.nqx, sgrid->IN_ti.ncp, 'd');     

#ifdef BE_CARBYDE
	 /* SDcmt: Write the qux buffer for MlMatrOut,*/

	CFORT_SORTCARBIDETRIM_SIMPLE(&(y_snet->qux[0][0]),&(Be2C_concentration[0][0]));
	 for(elem_no=0;elem_no<sgrid->IN_ti.ncp;elem_no++){ /*SDcmt: REDUNDANT data for every Element this is BeC_All!!*/
       for (slab_no=0;slab_no<sgrid->IN_ti.nqx;slab_no++){
		   qux_buf[slab_no+elem_no*sgrid->IN_ti.nqx]=Be2C_concentration[BeC_All][slab_no];
       }
     }
	 fprintf(fMout[1],"%% Carbide concentration \n");
     sprintf(sVarN,"CellProp{%i}.data(%i,%i,:,:)",CP_BE2C_TRIM,x_cell+1,y_cell+1);
     MlMatrOut (fMout[1], sVarN, qux_buf, sgrid->IN_ti.nqx, sgrid->IN_ti.ncp, 'd');     
#endif /*BE_CARBYDE*/
   }
 }

 if (sgrid->surface_model==3){
 fprintf (fMout[2], "%%\n%% -----------------------------------------------------------------------------\n");
 fprintf (fMout[2], "%%\n%% ----------------- START OF FOLLOWED PARTICLES DISTRIBUTIONS -----------------\n");
 
 sprintf(sVarN,"e0_bins");
 MlMatrOut (fMout[2], sVarN, sgrid->e0_bins,N_E0_BINS,1,'d');

 sprintf(sVarN,"alpha0_bins");
 MlMatrOut (fMout[2], sVarN, sgrid->alpha0_bins,N_ALPHA0_BINS,1,'d');

 for(x_cell=0;x_cell<iNCx;x_cell++){		
	 for(y_cell=0;y_cell<iNCy;y_cell++){ 	
		 y_snet = &(netz[x_cell][y_cell]);
		 /* SDcmt: Write the distrib buffer for MlMatrOut,*/
		 fprintf (fMout[2], "%%\n%% ---------------- FOLLOWED PARTICLES DISTRIBUTION DATA OF CELL (%i,%i) ---------------\n",x_cell+1,y_cell+1);		
		 
		 for(elem_no=0;elem_no<sgrid->IN_ti.ncp;elem_no++){		
			 fprintf(fMout[2],"CellProp{%i}.extra_data(%i,%i,%i).number_of_distribs=%i;\n",CP_MAX+1,x_cell+1,y_cell+1,elem_no+1,y_snet->number_of_distribs[elem_no]);
			 for(distrib_no=0;distrib_no<y_snet->number_of_distribs[elem_no];distrib_no++){
				 fprintf(fMout[2],"CellProp{%i}.name='dDistrib_%i';\n",CP_MAX+1+distrib_no,1+distrib_no);
		
				 fprintf(fMout[2],"CellProp{%i}.extra_data(%i,%i,%i).suppress_bsc=%i;\n",CP_MAX+1+distrib_no,x_cell+1,y_cell+1,elem_no+1,y_snet->IN_followed_part[distrib_no*(MAX_N_ELEM)+elem_no].su_bsc);		
				 fprintf(fMout[2],"CellProp{%i}.extra_data(%i,%i,%i).suppress_dep=%i;\n",CP_MAX+1+distrib_no,x_cell+1,y_cell+1,elem_no+1,y_snet->IN_followed_part[distrib_no*(MAX_N_ELEM)+elem_no].su_dep);		
				 fprintf(fMout[2],"CellProp{%i}.extra_data(%i,%i,%i).suppress_sput=%i;\n",CP_MAX+1+distrib_no,x_cell+1,y_cell+1,elem_no+1,y_snet->IN_followed_part[distrib_no*(MAX_N_ELEM)+elem_no].su_sput);						

				 sprintf(sVarN,"CellProp{%i}.data(%i,%i,%i,:,:)",CP_MAX+1+distrib_no,x_cell+1,y_cell+1,elem_no+1);
				 MlMatrOut (fMout[2], sVarN, y_snet->IN_followed_part[distrib_no*(MAX_N_ELEM)+elem_no].fluence, N_ALPHA0_BINS, N_E0_BINS, 'd');     						 	 
			 }
		 }
	 }
 }
 }
 for (iFiles=0;iFiles<N_OUTPUT_FILES;iFiles++)
	fclose (fMout[iFiles]);
  free(qux_buf);			
}
/*end of the routine  matlab_trim_out !!!*/


void Add_PART_quant(struct RES_particle_quantities *PQsum,struct RES_particle_quantities *PQ1,struct RES_particle_quantities *PQ2,int iActNElemTRIM){
  int i;

  for (i=0;i<iActNElemTRIM;i++){
    PQsum[i].dFlc=PQ1[i].dFlc + PQ2[i].dFlc ;
    PQsum[i].dReflected=PQ1[i].dReflected + PQ2[i].dReflected ;
    PQsum[i].dReemitted=PQ1[i].dReemitted + PQ2[i].dReemitted ;
	PQsum[i].dDeletedComp=PQ1[i].dDeletedComp + PQ2[i].dDeletedComp ;
    PQsum[i].dSputtered=PQ1[i].dSputtered + PQ2[i].dSputtered ;
    PQsum[i].dSputtered_chem=PQ1[i].dSputtered_chem + PQ2[i].dSputtered_chem ;
    PQsum[i].dDeposited=PQ1[i].dDeposited + PQ2[i].dDeposited ; 
    PQsum[i].dTransmitted=PQ1[i].dTransmitted + PQ2[i].dTransmitted ;
    PQsum[i].dStoppedRecoil=PQ1[i].dStoppedRecoil + PQ2[i].dStoppedRecoil ;
    PQsum[i].dStoppedProj=PQ1[i].dStoppedProj + PQ2[i].dStoppedProj ;
    PQsum[i].dResource=PQ2[i].dResource ;
    PQsum[i].dReflected_chem=PQ1[i].dReflected_chem + PQ2[i].dReflected_chem;
  }
}

#endif /*ifdef TRIM*/
#ifdef TRIM
/*--SDcmt----------------------------------------------------------------------------------------*/
/*                     function matlab_trim_in(...)                                              */
/* Get the available surface information for SDTrimSP                                            */
/*                                                                                               */
/*--SDcmt----------------------------------------------------------------------------------------*/
void  matlab_trim_in( struct SurfaceCell       **s_net,  /*heisst im move_ctrl() "s_net"                */
		      struct Path               *pfade,  /*heisst im move_ctrl() "pfade"              */
		      struct TimeCell         *timer,  /*heisst im move_ctrl() "time_ptr"           */
		      struct LostParticle    *looser,  /*heisst im move_ctrl() "lost"               */                 
		      struct SurfaceGridData *grid_d,  /*heisst im move_ctrl() "sgrid" (aus main.c) */   
		      int                      *rand,  /*heisst im move_ctrl() "grid_rand"          */
		      double                  *bound,  /*heisst im move_ctrl() "grid_bound"         */
		      int                 time_g   /*heisst im move_ctrl() "time_g"             */
){
  /*---------------------------------------------------------------------------
    Variablendefinition
  ----------------------------------------------------------------------------*/
  int i,x_cell,y_cell,iNElm,elem_no,distrib_no,iFiles;
  double **qux_buf;
  double *dBuf,*xno_buf;
  struct SurfaceCell       *y_snet;
  char sFN_NTrv[N_OUTPUT_FILES][DEF_STR_LEN];
  char sKWrd[DEF_STR_LEN];
  FILE *filehandler[N_OUTPUT_FILES]; /* for fileopen if file doesn't exists */
  long lqux=0,ldns=0,lTotDepth=0,lNewEn=0,lTotR=0, lDistrib=0,ldummy=-1;
  char *merged;
  /*---------------------------------------------------------------------------
    Ende von Variablendefinition
  ----------------------------------------------------------------------------*/
  sprintf (sFN_NTrv[0],"%sTrim_step%i",pfade->RESULT,time_g-1);
  sprintf (sFN_NTrv[1],"%sTrim3D_step%i",pfade->RESULT,time_g-1);
  sprintf (sFN_NTrv[2],"%sTrimDistr_step%i",pfade->RESULT,time_g-1);

   for(iFiles=0;iFiles<N_OUTPUT_FILES;iFiles++)
     {
       merged=sMerge(sFN_NTrv[iFiles],".m");
       filehandler[iFiles]=fopenPlus(merged);
     }



  /*-----------------------------------------------------
    zweifache for()- Schleife:
  ------------------------------------------------------*/
  for(x_cell=0;x_cell<rand[1];x_cell++) {
    for(y_cell=0;y_cell<rand[0];y_cell++){
      y_snet = &(s_net[x_cell][y_cell]);
      
      /*Read qux from matlab style ascii file*/
      sprintf(sKWrd,"CellProp{%i}.data(%i,%i,:,:)",CP_QUX,x_cell+1,y_cell+1);
      if ((iNElm=iFBlockRead (filehandler[1], sKWrd, CT_MATR, &qux_buf,&lqux))!=grid_d->IN_ti.nqx*grid_d->IN_ti.ncp){
	printf ("Error in prozesses: inappropriate number of elements (%i) in matrix '%s'.\n"
		,iNElm,sKWrd);
	exit(1);
      }
      /*Re-organise*/
      for(elem_no=0;elem_no<grid_d->IN_ti.ncp;elem_no++){
	memcpy (y_snet->qux[elem_no],qux_buf[elem_no],grid_d->IN_ti.nqx *sizeof(double));
      }
      free(qux_buf);

      /*Read xno from matlab style ascii file*/
      if (x_cell==0 && y_cell==0){ /*read only once*/
	sprintf(sKWrd,"CellProp{%i}.data(%i,%i,:)",CP_XNO,x_cell+1,y_cell+1);
	if ((iNElm=iFBlockRead (filehandler[1], sKWrd, CT_VEC, &xno_buf,&ldummy))!=(grid_d->IN_ti.nqx+1)){
	  printf ("Error in prozesses: inappropriate number of elements (%i) in vector '%s'.\n"
		  ,iNElm,sKWrd);
	  exit(1);
	}
      }
      /*Re-organise*/
      memcpy (y_snet->xno,xno_buf,(grid_d->IN_ti.nqx+1) *sizeof(double));
      
      /*Calculate middle of slabs xxx from beginning of slabs xno:*/
      for (i=0;i<grid_d->IN_ti.nqx;i++){
	y_snet->xxx[i]=y_snet->xno[i]+(y_snet->xno[i+1]-y_snet->xno[i])/2;
      }

     /*Read dns from matlab style ascii file*/
      sprintf(sKWrd,"CellProp{%i}.data(%i,%i,:)",CP_DNS,x_cell+1,y_cell+1);
      if ((iNElm=iFBlockRead (filehandler[1], sKWrd, CT_VEC, &dBuf,&ldns))!=(grid_d->IN_ti.nqx)){
	printf ("Error in prozesses: inappropriate number of elements (%i) in vector '%s'.\n"
		,iNElm,sKWrd);
	exit(1);
      }
      /*Re-organise*/
      memcpy (y_snet->dns,dBuf,(grid_d->IN_ti.nqx) *sizeof(double));
      free(dBuf);
    }/*y_cell*/
  }/*x_cell*/
  free(xno_buf);

   for(elem_no=0;elem_no<grid_d->IN_ti.ncp;elem_no++){
  /*
   * Read sputtered from followed particles (FO) for newenergetics. 
   *                               FO* TRIMPARTSTATNO + sputtered + matlab index starts @ 1
   *   --> SurfCell(element_number,1 * TRIMPARTSTATNO + 2         + 1)
   */
     if (grid_d->surface_model!=3){/*new energetics automatically included in SM3...*/
     sprintf(sKWrd,"SurfCell(%i,%i,:,:)",elem_no+1,1*TRIMPARTSTATNO+SC_SPUTTERED+1);   
     if ((iNElm=iFBlockRead (filehandler[0], sKWrd, CT_MATR, &qux_buf,&lNewEn))!=rand[0]*rand[1]){
       printf ("Error in prozesses: inappropriate number of elements (%i) in matrix '%s'.\n"
	       ,iNElm,sKWrd);
       exit(1);
     }
	 }
     /*Re-organise*/
	 for(x_cell=0;x_cell<rand[1];x_cell++) {
		 for(y_cell=0;y_cell<rand[0];y_cell++){
			 y_snet = &(s_net[x_cell][y_cell]);
			 if (grid_d->surface_model!=3){
				 y_snet->newenergetic[elem_no]=qux_buf[y_cell][x_cell]*y_snet->area;
			 }else if (grid_d->surface_model==3){
				 y_snet->newenergetic[elem_no]=0.;
			 }
		 }
	 }
  /*
   * Read total resource from particles (TOT)
   *                               TOT* TRIMPARTSTATNO + Resource  + matlab index starts @ 1
   *   --> SurfCell(element_number,0 *  TRIMPARTSTATNO + 9         + 1)
   */
     sprintf(sKWrd,"SurfCell(%i,%i,:,:)",elem_no+1,SC_TOTAL_RESOURCE+1);   
     if ((iNElm=iFBlockRead (filehandler[0], sKWrd, CT_MATR, &qux_buf,&lTotR))!=rand[0]*rand[1]){
       printf ("Error in prozesses: inappropriate number of elements (%i) in matrix '%s'.\n"
	       ,iNElm,sKWrd);
       exit(1);
     }
     /*Re-organise*/
     for(x_cell=0;x_cell<rand[1];x_cell++) {
       for(y_cell=0;y_cell<rand[0];y_cell++){
	 y_snet = &(s_net[x_cell][y_cell]);
	 y_snet->PART_quant[0][elem_no].dTotResource=qux_buf[y_cell][x_cell]*y_snet->area;
       }
     }

     free(qux_buf);
   } /*for (elemno...)*/

   /*Get the total erosion/deposition*/
   sprintf(sKWrd,"CellProp{%i}.data",CP_DTOTDELTADEPTH);   
   if ((iNElm=iFBlockRead (filehandler[0], sKWrd, CT_MATR, &qux_buf,&lTotDepth))!=rand[0]*rand[1]){
     printf ("Error in prozesses: inappropriate number of elements (%i) in matrix '%s'.\n"
	     ,iNElm,sKWrd);
     exit(1);
   }
   /*Re-organise*/
   for(x_cell=0;x_cell<rand[1];x_cell++) {
     for(y_cell=0;y_cell<rand[0];y_cell++){
       y_snet = &(s_net[x_cell][y_cell]);
       y_snet->dTotDeltaDepth=qux_buf[y_cell][x_cell];/*Transpose array!*/
     }
   }
   free(qux_buf);

   /*Get the distributions*/
   if (grid_d->surface_model==3){
	   sprintf(sKWrd,"e0_bins");
	   if ((iNElm=iFBlockRead (filehandler[2], sKWrd, CT_VEC_NO_ALLOC, grid_d->e0_bins,&lDistrib))!=(N_E0_BINS)){
		   printf ("Error in prozesses: inappropriate number of elements (%i) in matrix '%s'.\n"
			   ,iNElm,sKWrd);
		   exit(1);
	   }
	   sprintf(sKWrd,"alpha0_bins");
	   if ((iNElm=iFBlockRead (filehandler[2], sKWrd, CT_VEC_NO_ALLOC, grid_d->alpha0_bins,&lDistrib))!=(N_ALPHA0_BINS)){
		   printf ("Error in prozesses: inappropriate number of elements (%i) in matrix '%s'.\n"
			   ,iNElm,sKWrd);
		   exit(1);
	   }

     for(x_cell=0;x_cell<rand[1];x_cell++){		
	   for(y_cell=0;y_cell<rand[0];y_cell++){ 	
		   y_snet = &(s_net[x_cell][y_cell]);
		   /* SDcmt: Write the distrib buffer for MlMatrOut,*/

		   for(elem_no=0;elem_no<grid_d->IN_ti.ncp;elem_no++){	
			   sprintf(sKWrd,"CellProp{%i}.extra_data(%i,%i,%i).number_of_distribs",CP_MAX+1,x_cell+1,y_cell+1,elem_no+1);   
			   if ((iNElm=iFBlockRead (filehandler[2], sKWrd, CT_PARI, &(y_snet->number_of_distribs[elem_no]),&lDistrib))!=1){
				   printf ("Error in prozesses: inappropriate number of elements (%i) in matrix '%s'.\n"
					   ,iNElm,sKWrd);
				   exit(1);
			   }
			   
			   for(distrib_no=0;distrib_no<y_snet->number_of_distribs[elem_no];distrib_no++){
				   sprintf(sKWrd,"CellProp{%i}.extra_data(%i,%i,%i).suppress_bsc",CP_MAX+1+distrib_no,x_cell+1,y_cell+1,elem_no+1);   
				   if ((iNElm=iFBlockRead (filehandler[2], sKWrd, CT_PARI, &(y_snet->IN_followed_part[distrib_no*(MAX_N_ELEM)+elem_no].su_bsc),&lDistrib))!=1){
					   printf ("Error in prozesses: inappropriate number of elements (%i) in matrix '%s'.\n"
						   ,iNElm,sKWrd);
					   exit(1);
				   }
				   sprintf(sKWrd,"CellProp{%i}.extra_data(%i,%i,%i).suppress_dep",CP_MAX+1+distrib_no,x_cell+1,y_cell+1,elem_no+1);   
				   if ((iNElm=iFBlockRead (filehandler[2], sKWrd, CT_PARI, &(y_snet->IN_followed_part[distrib_no*(MAX_N_ELEM)+elem_no].su_dep),&lDistrib))!=1){
					   printf ("Error in prozesses: inappropriate number of elements (%i) in matrix '%s'.\n"
						   ,iNElm,sKWrd);
					   exit(1);
				   }
				   sprintf(sKWrd,"CellProp{%i}.extra_data(%i,%i,%i).suppress_sput",CP_MAX+1+distrib_no,x_cell+1,y_cell+1,elem_no+1);   
				   if ((iNElm=iFBlockRead (filehandler[2], sKWrd, CT_PARI, &(y_snet->IN_followed_part[distrib_no*(MAX_N_ELEM)+elem_no].su_sput),&lDistrib))!=1){
					   printf ("Error in prozesses: inappropriate number of elements (%i) in matrix '%s'.\n"
						   ,iNElm,sKWrd);
					   exit(1);
				   }
				   sprintf(sKWrd,"CellProp{%i}.data(%i,%i,%i,:,:)",CP_MAX+1+distrib_no,x_cell+1,y_cell+1,elem_no+1);   
				   if ((iNElm=iFBlockRead (filehandler[2], sKWrd, CT_MATR_NO_ALLOC, (y_snet->IN_followed_part[distrib_no*(MAX_N_ELEM)+elem_no].fluence),&lDistrib))!=(N_ALPHA0_BINS)*(N_E0_BINS)){
					   printf ("Error in prozesses: inappropriate number of elements (%i) in matrix '%s'.\n"
						   ,iNElm,sKWrd);
					   exit(1);
				   }				   				  
			   }/*for distrib_no*/
		   }/*for elem_no*/
	   }/*for y_cell*/
     }/*for x_cell*/
   }/*endif surface_model==3*/
   for(iFiles=0;iFiles<N_OUTPUT_FILES;iFiles++)
	fclose (filehandler[iFiles]);
   printf ("Data for %s  was succesfully read\n\n",sFN_NTrv);
}
/*   ENDOF  function matlab_trim_in(...)                                                          */
#endif /*TRIM*/

/*----------------------------------------------------------------
Prints out E_field in Matlab format
-----------------------------------------------------------------*/
void MatlEFOut( 
	     struct Path          *pfade,
		 double		          *dBound,
	     struct  Numeric           *num,
	     struct Ort                *ort,  /* contains limiter info                        */
	     struct Electric       *e_param,  /* therefore we are here                        */
	     struct BField               *b,  
	     struct Magnetic      *magnetic,
	     struct Plasma           *p_ptr,
             struct nTdata          *nT_ptr,
	     struct PConst      *plas_const,
#ifdef PISCES
	     struct SurfaceCell       **net_ptr,  
		 int                         *rand,  
#endif /* PISCES */
	     int                     surf
	   )
{
#ifdef PISCES
   const int iNP=21;
#else
   const int iNP=41;
#endif
#ifndef UGEOMETRY
   const double dRS=1.2;
#else
   const double dRS=1.0;
   /*struct EroSurfaceOutputData eroData;*/
#endif

   double dx,dy,dz, dxL,dyL,dzL, dX1, dY1, dZ1;
   double dXL, dXR, dYL, dYR;
   int ix,iy,iz;
   double *dBufDz, *dBufEx, *dBufEy, *dBufEz, *dX, *dY, *dZ;
   double dZmin=0, dZmax;

   struct ParticleStruct   w;
   struct EField           e;
  
   FILE *fMOut;
   char sVarN[DEF_STR_LEN];

   double progress;

   printf("Writing E_field . . .\n");
     
#ifdef PISCES
   dZmax=0.5;

   dXL = -50;
   dXR = 50;
   dYL = -50;
   dYR = 50;
#else /* PISCES */
   
#ifndef UGEOMETRY
   if ((ort->t_rad>0)||(ort->p_rad>0))
     dZmax=maximum(ort->t_rad, ort->p_rad);
   else
     dZmax=ort->angle*ort->t_len;
#else
   dZmin = dBound[4];
   dZmax = dBound[5];
#endif

#ifdef JETLim
   dZmax = 200.;
#endif

   dXL = dBound[0];
   dXR = dBound[1];
   dYL = dBound[2];
   dYR = dBound[3];
#endif /* PISCES */

#ifdef JETLim
   dXL = -200;/*dBound[0];*/
   dXR = +200;/*dBound[1];*/
   dYL = -170;/*dBound[2];*/
   dYR = +170;/*dBound[3];*/
#endif

   dxL=dXR-dXL;
   dx=dxL/(iNP-1)*dRS;
   dX1=dXL - dxL*(dRS-1)/2;
   dyL=dYR - dYL;
   dy=dyL/(iNP-1)*dRS;
   dY1=dYL - dyL*(dRS-1)/2;
   
   dzL=dZmax-dZmin;
   dz=dzL/(iNP-1)*dRS;
   dZ1=dZmin - dzL*(dRS-1)/2;

   fMOut= MFOpen ("EField", pfade->RESULT, -1);

   fprintf(fMOut,"iNP = %i;\n%%\n",iNP);
   fprintf(fMOut,"iSType = %i;\n%%\n",surf);

   fprintf(fMOut,"dx = %.2e;\n",dx);
   fprintf(fMOut,"dy = %.2e;\n",dy);
   fprintf(fMOut,"dz = %.2e;\n",dz);

   dBufDz=(double *) calloc (iNP*iNP,sizeof(double));
   dBufEx=(double *) calloc (iNP*iNP,sizeof(double));
   dBufEy=(double *) calloc (iNP*iNP,sizeof(double));
   dBufEz=(double *) calloc (iNP*iNP,sizeof(double));
   dX=(double *) calloc (iNP,sizeof(double));
   dY=(double *) calloc (iNP,sizeof(double));
   dZ=(double *) calloc (iNP,sizeof(double));

   for (ix=0; ix<iNP; ix++){
	 dX[ix]=dX1+ix*dx;
	 dY[ix]=dY1+ix*dy;
	 dZ[ix]=dZ1+ix*dz;
   }
   MlMatrOut (fMOut, "dX", dX, iNP, 1, 'd');
   MlMatrOut (fMOut, "dY", dY, iNP, 1, 'd');
   MlMatrOut (fMOut, "dZ", dZ, iNP, 1, 'd');

   printf("output progress . . . 00%%");

   for (iz=0; iz<iNP; iz++){

	 progress = (double)iz/(iNP-1)*100;
	 printf("\r\r\routput progress . . . %03d%%", (int)progress);

     for (ix=0; ix<iNP; ix++)
	 for (iy=0; iy<iNP; iy++){

       w.location[0]=dX[ix];
	   w.location[1]=dY[iy];
       w.location[2]=dZ[iz];

	   #ifdef UGEOMETRY
	   w.eroSData.voxel = -1;
	   w.eroSData.distance = NOVALUE;
	   w.eroSData.normal[0] = -1./sqrt(3);
	   w.eroSData.normal[1] = -1./sqrt(3);
	   w.eroSData.normal[2] = -1./sqrt(3);
	   w.eroSData.normalInterp[0] = -1./sqrt(3);
	   w.eroSData.normalInterp[1] = -1./sqrt(3);
	   w.eroSData.normalInterp[2] = -1./sqrt(3);
	   w.eroSData.proxLim = num->proxLim;
	   w.eroSData.projection[0] = dBound[0] - 100.; /* to guarantee position outside the box */
	   w.eroSData.projection[1] = dBound[2] - 100.;
	   w.eroSData.projection[2] = dBound[4] - 100.;
	   w.eroSData.poly = NULL;
	   #endif

       E_field(&w,  /* Arbeitszeiger                                */
	     num,
	     ort,  /* contains limiter info                        */
	     e_param,  /* therefore we are here                        */
	     b,  
	     magnetic,
	     &e,
	     p_ptr,
         nT_ptr,
	     plas_const,
#ifdef PISCES
	     net_ptr,
		 rand,  
         dBound,
#endif /* PISCES */
#ifdef UGEOMETRY
	     0
#else
		 surf
#endif
	    );

	   dBufDz[ix+iNP*iy]=e.dz;
	   dBufEx[ix+iNP*iy]=e.E[0];
	   dBufEy[ix+iNP*iy]=e.E[1];
	   dBufEz[ix+iNP*iy]=e.E[2];
	 }

	 sprintf(sVarN,"dE_dz(%i,:,:)",iz+1);
	 MlMatrOut (fMOut, sVarN, dBufDz, iNP, iNP, 'd');
	 sprintf(sVarN,"dE_x(%i,:,:)",iz+1);
	 MlMatrOut (fMOut, sVarN, dBufEx, iNP, iNP, 'd');
	 sprintf(sVarN,"dE_y(%i,:,:)",iz+1);
	 MlMatrOut (fMOut, sVarN, dBufEy, iNP, iNP, 'd');
	 sprintf(sVarN,"dE_z(%i,:,:)",iz+1);
	 MlMatrOut (fMOut, sVarN, dBufEz, iNP, iNP, 'd');
   }

#ifdef JETLim
  for (ix=0; ix<iNP; ix++){
	 dZ[ix]=0+ 0.001*pow(1.3,ix);
	 w.location[0] = 0;
	 w.location[1] = 0;
	 w.location[2] = BM_dJETLimShape(0., 0., NULL)+ dZ[ix];
	 E_field(&w,  /* Arbeitszeiger                                */
	     num,
	     ort,  /* contains limiter info                        */
	     e_param,  /* therefore we are here                        */
	     b,  
	     magnetic,
	     &e,
	     p_ptr,
         nT_ptr,
	     plas_const,
#ifdef PISCES
	     net_ptr,
		 rand,  
         dBound,
#endif /* PISCES */
#ifdef UGEOMETRY
	     0
#else
		 surf
#endif
	    );
	 dBufEz[ix]=e.E[2];
	 dBufDz[ix]=e.dz;
   }
   MlMatrOut (fMOut, "dZ_00z", dZ, iNP, 1, 'd');
   MlMatrOut (fMOut, "d_dz_00z", dBufDz, iNP, 1, 'd');
   MlMatrOut (fMOut, "dE_00z", dBufEz, iNP, 1, 'd');
#endif

   free(dBufDz);
   free(dBufEx);
   free(dBufEy);
   free(dBufEz);
   free(dX);
   free(dY);
   free(dZ);
   fclose (fMOut);
   printf("\nWriting of E_field complete.\n");
}

/*----------------------------------------------------------------
Prints out Plasma parmaeters (ne, Te, Ti) in Matlab format
-----------------------------------------------------------------*/
void MatlPlParOut( 
	     struct Path          *pfade,
		 double		          *dBound,
	     struct  Numeric           *num,
	     struct Ort                *ort,  /* contains limiter info                        */
	     struct Electric       *e_param,  /* therefore we are here                        */
	     struct BField               *b,  
	     struct Magnetic      *magnetic,
	     struct Plasma           *p_ptr,
             struct nTdata          *nT_ptr,
	     struct PConst      *plas_const,
	     int                     surf,
		 char  cflRM,
		 struct VolumeAll  *v_all,
		 int    iTStep
	   )
{
   const int iNP=21;
#ifndef UGEOMETRY
   const double dRS=1.5;
#else
   const double dRS=1.0;
   struct EroBfieldOutputData eroBData;
#endif
   double dx,dy,dz, dxL,dyL,dzL, dX1, dY1, dZ1;
   int iNPx, iNPy, iNPz;
   int ix,iy,iz;
   double *dBufNe, *dBufTe, *dBufTi, *dX, *dY, *dZ;
   double dZmin=0, dZmax;
   /*AL ORNL 2016 print grad Te, grad Ti 
   double *dBufgradTe_x, *dBufgradTe_y, *dBufgradTe_z;
   double *dBufgradTi_x, *dBufgradTi_y, *dBufgradTi_z;*/
   
   FILE *fMOut;
   char sVarN[DEF_STR_LEN];

   printf("Writing Plasma Parameters . . .\n");

#ifndef UGEOMETRY
   if ((ort->t_rad>0)||(ort->p_rad>0))
     dZmax=maximum(ort->t_rad, ort->p_rad);
   else
     dZmax=ort->angle*ort->t_len;
#else
   dZmin = dBound[4];
   dZmax = dBound[5];
#endif

  if  (cflRM == 1){ 
   iNPx=iNP; iNPy=iNP; iNPz=iNP;
   dxL=dBound[1]-dBound[0];
   dx=dxL/(iNP-1)*dRS;
   dX1=dBound[0]- dxL*(dRS-1)/2;
   dyL=dBound[3]-dBound[2];
   dy=dyL/(iNP-1)*dRS;
   dY1=dBound[2]- dyL*(dRS-1)/2;
   
   dzL=dZmax-dZmin;
   dz=dzL/(iNP-1)*dRS;
   dZ1=dZmin - dzL*(dRS-1)/2;
  }
  else{
    dx          = v_all->diff[0];
    dy          = v_all->diff[1];
    dz          = v_all->diff[2];
    dX1       = v_all->adRange[0];
    dY1       = v_all->adRange[1]; 
	dZ1       = 0;
	iNPx      = v_all->anz[0];    /* points in toroidal dir. */
    iNPy      = v_all->anz[1];    /* points in poloidal dir. */
    iNPz      = v_all->anz[2];
  }
  
  if  (cflRM == 1){ 
     fMOut= MFOpen ("PlParRM", pfade->RESULT, iTStep);
     fprintf(fMOut,"iNP = %i;\n%%\n",iNP);
  }else{
	 fMOut= MFOpen ("PlPar", pfade->RESULT, iTStep);
	 fprintf(fMOut,"iNPx = %i;\n"    ,iNPx);
	 fprintf(fMOut,"iNPy = %i;\n"    ,iNPy);
	 fprintf(fMOut,"iNPz = %i;\n%%\n",iNPz);
  }

   fprintf(fMOut,"iSType = %i;\n%%\n",surf);

   fprintf(fMOut,"dx = %.2e;\n",dx);
   fprintf(fMOut,"dy = %.2e;\n",dy);
   fprintf(fMOut,"dz = %.2e;\n",dz);

   dBufNe=(double *) calloc (iNPx*iNPy,sizeof(double));
   dBufTe=(double *) calloc (iNPx*iNPy,sizeof(double));
   dBufTi=(double *) calloc (iNPx*iNPy,sizeof(double));
   /*AL ORNL 2016 print grad Te, grad Ti 
   dBufgradTe_x=(double *) calloc (iNPx*iNPy,sizeof(double));
   dBufgradTe_y=(double *) calloc (iNPx*iNPy,sizeof(double));
   dBufgradTe_z=(double *) calloc (iNPx*iNPy,sizeof(double));
   dBufgradTi_x=(double *) calloc (iNPx*iNPy,sizeof(double));
   dBufgradTi_y=(double *) calloc (iNPx*iNPy,sizeof(double));
   dBufgradTi_z=(double *) calloc (iNPx*iNPy,sizeof(double));
   */
   
   dX=(double *) calloc (iNPx,sizeof(double));
   dY=(double *) calloc (iNPy,sizeof(double));
   dZ=(double *) calloc (iNPz,sizeof(double));

   for (ix=0; ix<iNPx; ix++)
	 dX[ix]=dX1+ix*dx;
   for (ix=0; ix<iNPy; ix++)
	 dY[ix]=dY1+ix*dy;
   for (ix=0; ix<iNPz; ix++)
	 dZ[ix]=dZ1+ix*dz;
  
   MlMatrOut (fMOut, "dX", dX, iNPx, 1, 'd');
   MlMatrOut (fMOut, "dY", dY, iNPy, 1, 'd');
   MlMatrOut (fMOut, "dZ", dZ, iNPz, 1, 'd');

   for (iz=0; iz<iNPz; iz++){
     for (ix=0; ix<iNPx; ix++)
	 for (iy=0; iy<iNPy; iy++){

#ifdef DIV_NEWGEOM
           /* Eetu Ahonen 120811: added struct Plasma p_ptr to arguments */
	   dBufTe[ix+iNPx*iy]=
	     Te(dX[ix], dY[iy], dZ[iz], nT_ptr, num, magnetic, ort, b, plas_const, p_ptr, surf);
	   dBufTi[ix+iNPx*iy]=
		  Ti(dX[ix], dY[iy], dZ[iz], nT_ptr, num, magnetic, ort, b, plas_const, p_ptr, surf, dBufTe[ix+iNPx*iy]);
#else
	   dBufTe[ix+iNPx*iy]=
	     Te(dX[ix], dY[iy], dZ[iz], nT_ptr, num, magnetic, ort, b, plas_const, surf);
	   dBufTi[ix+iNPx*iy]=
		  Ti(dX[ix], dY[iy], dZ[iz], nT_ptr, num, magnetic, ort, b, plas_const, surf, dBufTe[ix+iNPx*iy]);
#endif
#ifdef UGEOMETRY
	   eroBData.SOL = -1;
	   eroBData.SHADOW = -1;
	   dBufNe[ix+iNPx*iy]=
		  ne(dX[ix], dY[iy], dZ[iz], nT_ptr, num, magnetic, ort, b, plas_const, p_ptr, surf, &eroBData);
#else
	   dBufNe[ix+iNPx*iy]=
		  ne(dX[ix], dY[iy], dZ[iz], nT_ptr, num, magnetic, ort, b, plas_const, p_ptr, surf);
#endif

	   /*AL ORNL 2016, print grad Te, grad Ti, as used in thermal force calculation 

	     grad_Te(dX[ix], dY[iy], dZ[iz], nT_ptr, num, p_ptr, magnetic, ort, b, plas_const, surf);
	     grad_Ti(dX[ix], dY[iy], dZ[iz], nT_ptr, num, p_ptr, magnetic, ort, b, plas_const, surf

#ifdef UGEOMETRY
		     , &eroBData
#endif
		     );
	     

	     dBufgradTe_x[ix+iNPx*iy]=nT_ptr->grad_Te[0];
	     dBufgradTe_y[ix+iNPx*iy]=nT_ptr->grad_Te[1];
	     dBufgradTe_z[ix+iNPx*iy]=nT_ptr->grad_Te[2];

	     dBufgradTi_x[ix+iNPx*iy]=nT_ptr->grad_Ti[0];
	     dBufgradTi_y[ix+iNPx*iy]=nT_ptr->grad_Ti[1];
	     dBufgradTi_z[ix+iNPx*iy]=nT_ptr->grad_Ti[2];
	   */	     
	 }

	 sprintf(sVarN,"Ne(%i,:,:)",iz+1);
	 MlMatrOut (fMOut, sVarN, dBufNe, iNPx, iNPy, 'd');
	 sprintf(sVarN,"Te(%i,:,:)",iz+1);
	 MlMatrOut (fMOut, sVarN, dBufTe, iNPx, iNPy, 'd');
	 sprintf(sVarN,"Ti(%i,:,:)",iz+1);
	 MlMatrOut (fMOut, sVarN, dBufTi, iNPx, iNPy, 'd');
	 /*AL ORNL 2016 print grad Te, grad Ti 
	 sprintf(sVarN,"gradTe_x(%i,:,:)",iz+1);
	 MlMatrOut (fMOut, sVarN, dBufgradTe_x, iNPx, iNPy, 'd');
	 sprintf(sVarN,"gradTe_y(%i,:,:)",iz+1);
	 MlMatrOut (fMOut, sVarN, dBufgradTe_y, iNPx, iNPy, 'd');
	 sprintf(sVarN,"gradTe_z(%i,:,:)",iz+1);
	 MlMatrOut (fMOut, sVarN, dBufgradTe_z, iNPx, iNPy, 'd');
	 sprintf(sVarN,"gradTi_x(%i,:,:)",iz+1);
	 MlMatrOut (fMOut, sVarN, dBufgradTi_x, iNPx, iNPy, 'd');
	 sprintf(sVarN,"gradTi_y(%i,:,:)",iz+1);
	 MlMatrOut (fMOut, sVarN, dBufgradTi_y, iNPx, iNPy, 'd');
	 sprintf(sVarN,"gradTi_z(%i,:,:)",iz+1);
	 MlMatrOut (fMOut, sVarN, dBufgradTi_z, iNPx, iNPy, 'd');*/
   }

   free(dBufNe);
   free(dBufTe);
   free(dBufTi);
   /*AL ORNL 2016 print grad Te, grad Ti 
   free(dBufgradTe_x);
   free(dBufgradTe_y);
   free(dBufgradTe_z);
   free(dBufgradTi_x);
   free(dBufgradTi_y);
   free(dBufgradTi_z);*/
   
   free(dX);
   free(dY);
   free(dZ);

   fclose (fMOut);
   printf("Writing of plasma parameter complete.\n");
}

/*AL ORNL 2016; export ERO's volumetric input information for GITR*/

void Mtlb_GITRin(struct Path          *pfade,
		 double                   *dBound,
		 struct  Numeric           *num,
		 struct Ort                *ort,  /* contains limiter info                        */
		 struct Electric       *e_param,  /* therefore we are here                        */
		 struct BField               *b,
		 struct Magnetic      *magnetic,
		 struct Force             *force,
		 /*		 struct ParticleStruct   *w_ptr,*/
		 struct Plasma           *p_ptr,
		 struct nTdata          *nT_ptr,
		 struct PConst      *plas_const,
		 int                     surf,
		 struct VolumeAll  *v_all,
		 int    iTStep,
		 struct Randomstruct                     *idum
		 )
{
  const int iNP=21;
  const double dRS=1.0;
  
  double dx,dy,dz, dxL,dyL,dzL, dX1, dY1, dZ1;
  int iNPx, iNPy, iNPz;
  int ix,iy,iz;
  double *dBufNe, *dBufTe, *dBufTi, *dX, *dY, *dZ;

  double *dBufgradTe_x, *dBufgradTe_y, *dBufgradTe_z;
  double *dBufgradTi_x, *dBufgradTi_y, *dBufgradTi_z;
  double *dBufthermf_x, *dBufthermf_y, *dBufthermf_z;
  
  double *dBufEx, *dBufEy, *dBufEz;

  double *dBufcs;

  double loc_Te, loc_Ti;
  
  double dZmin=0, dZmax, zsurf;
  double eps=0.001; /*add small number when calling E_field, to avoid zone borderline issues*/
  double z_eff, mu, alpha, beta;
  double thermalf[3];
  
  int SOL;

  
  struct ParticleStruct   w;
  struct EField           e;
  
  FILE *fMOut;
  char sVarN[DEF_STR_LEN];

  w.mass=183.85;
  w.charge=1;
  
  printf("\nExporting plasma parameters for GITR . . .\n"); 
  printf("Includes:  ne, Te, Ti, grad_ne, grad_Te, grad_Ti, flow velocity, \n"); /*plus flow v when implemented*/ 
  printf("           E_field and thermal force (for mass %g charge %d) \n", w.mass, w.charge);
  
  if ((ort->t_rad>0)||(ort->p_rad>0))
    dZmax=maximum(ort->t_rad, ort->p_rad);
  else
    dZmax=ort->angle*ort->t_len;
  
  
  dx          = v_all->diff[0];
  dy          = v_all->diff[1];
  dz          = v_all->diff[2];
  dX1       = v_all->adRange[0];
  dY1       = v_all->adRange[1];
  dZ1       = 0;
  iNPx      = v_all->anz[0];    /* points in toroidal dir. */
  iNPy      = v_all->anz[1];    /* points in poloidal dir. */
  iNPz      = v_all->anz[2];
  zsurf=0; /*initialize*/
  
  dzL=dZmax-dZmin;
  /*dz=dzL/(iNP-1)*dRS;*/
  dZ1=dZmin - dzL*(dRS-1)/2;

  
  fMOut= MFOpen ("GITRinput", pfade->RESULT, iTStep);
  fprintf(fMOut,"iNPx = %i;\n"    ,iNPx);
  fprintf(fMOut,"iNPy = %i;\n"    ,iNPy);
  fprintf(fMOut,"iNPz = %i;\n%%\n",iNPz);
  
  
  fprintf(fMOut,"iSType = %i;\n%%\n",surf);
  
  fprintf(fMOut,"dx = %.2e;\n",dx);
  fprintf(fMOut,"dy = %.2e;\n",dy);
  fprintf(fMOut,"dz = %.2e;\n",dz);
  
  dBufNe=(double *) calloc (iNPx*iNPy,sizeof(double));
  dBufTe=(double *) calloc (iNPx*iNPy,sizeof(double));
  dBufTi=(double *) calloc (iNPx*iNPy,sizeof(double));
  
  dBufgradTe_x=(double *) calloc (iNPx*iNPy,sizeof(double));
  dBufgradTe_y=(double *) calloc (iNPx*iNPy,sizeof(double));
  dBufgradTe_z=(double *) calloc (iNPx*iNPy,sizeof(double));
  
  dBufgradTi_x=(double *) calloc (iNPx*iNPy,sizeof(double));
  dBufgradTi_y=(double *) calloc (iNPx*iNPy,sizeof(double));
  dBufgradTi_z=(double *) calloc (iNPx*iNPy,sizeof(double));

  dBufthermf_x=(double *) calloc (iNPx*iNPy,sizeof(double));
  dBufthermf_y=(double *) calloc (iNPx*iNPy,sizeof(double));
  dBufthermf_z=(double *) calloc (iNPx*iNPy,sizeof(double));

  
  dBufEx=(double *) calloc (iNPx*iNPy,sizeof(double));
  dBufEy=(double *) calloc (iNPx*iNPy,sizeof(double));
  dBufEz=(double *) calloc (iNPx*iNPy,sizeof(double));
  
  dBufcs=(double *) calloc (iNPx*iNPy,sizeof(double));
  
  dX=(double *) calloc (iNPx,sizeof(double));
  dY=(double *) calloc (iNPy,sizeof(double));
  dZ=(double *) calloc (iNPz,sizeof(double));

  
  for (ix=0; ix<iNPx; ix++)
    dX[ix]=dX1+ix*dx;
  for (ix=0; ix<iNPy; ix++)
    dY[ix]=dY1+ix*dy;
  for (ix=0; ix<iNPz; ix++)
    dZ[ix]=dZ1+ix*dz;
  

  MlMatrOut (fMOut, "dX", dX, iNPx, 1, 'd');
  MlMatrOut (fMOut, "dY", dY, iNPy, 1, 'd');
  MlMatrOut (fMOut, "dZ", dZ, iNPz, 1, 'd');

  for (iz=0; iz<iNPz; iz++){
    for (ix=0; ix<iNPx; ix++)
      for (iy=0; iy<iNPy; iy++){

	/* BG PLASMA: Ne, Te, Ti ; grad_Te, grad_Ti */
	
	dBufNe[ix+iNPx*iy]=
	  ne(dX[ix], dY[iy], dZ[iz], nT_ptr, num, magnetic, ort, b, plas_const, p_ptr, surf);
 
	
	loc_Te=Te(dX[ix], dY[iy], dZ[iz], nT_ptr, num, magnetic, ort, b, plas_const, surf);
	loc_Ti=Ti(dX[ix], dY[iy], dZ[iz], nT_ptr, num, magnetic, ort, b, plas_const, surf, dBufTe[ix+iNPx*iy]);
	dBufTe[ix+iNPx*iy]= loc_Te;
	dBufTi[ix+iNPx*iy]= loc_Ti;

	grad_Te(dX[ix], dY[iy], dZ[iz], nT_ptr, num, p_ptr, magnetic, ort, b, plas_const, surf);
	grad_Ti(dX[ix], dY[iy], dZ[iz], nT_ptr, num, p_ptr, magnetic, ort, b, plas_const, surf);

	dBufgradTe_x[ix+iNPx*iy]=nT_ptr->grad_Te[0];
	dBufgradTe_y[ix+iNPx*iy]=nT_ptr->grad_Te[1];
	dBufgradTe_z[ix+iNPx*iy]=nT_ptr->grad_Te[2];
	
	dBufgradTi_x[ix+iNPx*iy]=nT_ptr->grad_Ti[0];
	dBufgradTi_y[ix+iNPx*iy]=nT_ptr->grad_Ti[1];
	dBufgradTi_z[ix+iNPx*iy]=nT_ptr->grad_Ti[2];

	/* THERMAL FORCE -- BY grad_Te, grad_Ti */

	z_eff=z_effective(dX[ix], dY[iy], dZ[iz], p_ptr, nT_ptr);
	mu = w.mass /(w.mass + z_eff);

	alpha  = 0.71 * pow(w.charge,2.);
	beta = -3. *(  1 - mu - 5*M_SQRT2 * pow(w.charge,2.) * ( 1.1*pow(mu,2.5) - 0.35 *pow(mu,1.5) )
				/ (2.6 - 2. * mu + 5.4 *mu *mu ));

        thermalf[0] = 9.84425e12/w.mass * ( alpha * nT_ptr->grad_Te[0]  + beta * nT_ptr->grad_Ti[0] );
	thermalf[1] = 9.84425e12/w.mass * ( alpha * nT_ptr->grad_Te[1]  + beta * nT_ptr->grad_Ti[1] );
	thermalf[2] = 9.84425e12/w.mass * ( alpha * nT_ptr->grad_Te[2]  + beta * nT_ptr->grad_Ti[2] );
	
	
	dBufthermf_x[ix+iNPx*iy]=thermalf[0]*num->coll_t;
	dBufthermf_y[ix+iNPx*iy]=thermalf[1]*num->coll_t;
	dBufthermf_z[ix+iNPx*iy]=thermalf[2]*num->coll_t;



	/*FLOW VELOCITY */

	/*can't call local_flow_velocity(...) for the surface location at the surface zone limit -> add eps*/
	w.location[0]=dX[ix]+eps;
	w.location[1]=dY[iy]+eps;
	w.location[2]=dZ[iz]+eps;
	
	
	/*set to 1 cause it's needed in local_flow_Velocity routine; though used after computing cs;*/
	w.v[1]=0.2;
	w.v[2]=1.0;
	w.v[3]=1.0;

	
	SOL=local_flow_velocity( &w, num, p_ptr, nT_ptr, magnetic,  ort,  b , force,  plas_const, surf, num->coll_t, idum  );
	/*	plas_const->cs0=1.0;
		plas_const->cs=1.0;*/
	  
	dBufcs[ix+iNPx*iy]=plas_const->cs;


	/*E_FIELD*/
	
	/*can't call E_Field(...) for the surface location at the surface zone limit -> add eps*/

	w.location[0]=dX[ix]+eps;
	w.location[1]=dY[iy]+eps;
	w.location[2]=dZ[iz]+eps;
	
	E_field(&w,  
		num,
		ort, 
		e_param,
		b,
		magnetic,
		&e,
		p_ptr,
		nT_ptr,
		plas_const,
		surf
		);
	
	dBufEx[ix+iNPx*iy]=e.E[0];
	dBufEy[ix+iNPx*iy]=e.E[1];
	dBufEz[ix+iNPx*iy]=e.E[2];
	
      }
    
    
    sprintf(sVarN,"Ne(%i,:,:)",iz+1);
    MlMatrOut (fMOut, sVarN, dBufNe, iNPx, iNPy, 'd');
    sprintf(sVarN,"Te(%i,:,:)",iz+1);
    MlMatrOut (fMOut, sVarN, dBufTe, iNPx, iNPy, 'd');
    sprintf(sVarN,"Ti(%i,:,:)",iz+1);
    MlMatrOut (fMOut, sVarN, dBufTi, iNPx, iNPy, 'd');

    sprintf(sVarN,"gradTe_x(%i,:,:)",iz+1);
    MlMatrOut (fMOut, sVarN, dBufgradTe_x, iNPx, iNPy, 'd');
    sprintf(sVarN,"gradTe_y(%i,:,:)",iz+1);
    MlMatrOut (fMOut, sVarN, dBufgradTe_y, iNPx, iNPy, 'd');
    sprintf(sVarN,"gradTe_z(%i,:,:)",iz+1);
    MlMatrOut (fMOut, sVarN, dBufgradTe_z, iNPx, iNPy, 'd');

    sprintf(sVarN,"gradTi_x(%i,:,:)",iz+1);
    MlMatrOut (fMOut, sVarN, dBufgradTi_x, iNPx, iNPy, 'd');
    sprintf(sVarN,"gradTi_y(%i,:,:)",iz+1);
    MlMatrOut (fMOut, sVarN, dBufgradTi_y, iNPx, iNPy, 'd');
    sprintf(sVarN,"gradTi_z(%i,:,:)",iz+1);
    MlMatrOut (fMOut, sVarN, dBufgradTi_z, iNPx, iNPy, 'd');

    sprintf(sVarN,"thermF_x(%i,:,:)",iz+1);
    MlMatrOut (fMOut, sVarN, dBufthermf_x, iNPx, iNPy, 'd');
    sprintf(sVarN,"thermF_y(%i,:,:)",iz+1);
    MlMatrOut (fMOut, sVarN, dBufthermf_y, iNPx, iNPy, 'd');
    sprintf(sVarN,"thermF_z(%i,:,:)",iz+1);
    MlMatrOut (fMOut, sVarN, dBufthermf_z, iNPx, iNPy, 'd');


    sprintf(sVarN,"v_flow(%i,:,:)",iz+1);
    MlMatrOut (fMOut, sVarN, dBufcs, iNPx, iNPy, 'd');
    
    sprintf(sVarN,"dE_x(%i,:,:)",iz+1);
    MlMatrOut (fMOut, sVarN, dBufEx, iNPx, iNPy, 'd');
    sprintf(sVarN,"dE_y(%i,:,:)",iz+1);
    MlMatrOut (fMOut, sVarN, dBufEy, iNPx, iNPy, 'd');
    sprintf(sVarN,"dE_z(%i,:,:)",iz+1);
    MlMatrOut (fMOut, sVarN, dBufEz, iNPx, iNPy, 'd');
    
  }

  free(dBufNe);
  free(dBufTe);
  free(dBufTi);

  free(dBufgradTe_x);
  free(dBufgradTe_y);
  free(dBufgradTe_z);
  
  free(dBufgradTi_x);
  free(dBufgradTi_y);
  free(dBufgradTi_z);

  free(dBufthermf_x);
  free(dBufthermf_y);
  free(dBufthermf_z);

  free(dBufcs);
  
  free(dBufEx);
  free(dBufEy);
  free(dBufEz);
  
  free(dX);
  free(dY);
  free(dZ);

  fclose (fMOut);
  printf("Exporting plasma parameters for GITR complete\n");
}

/* AL ORNL 2016; end of export ERO's volumetric info for GITR*/


void matlab_particle_out(struct ParticleStruct *particles,
			 int loops,/*number of particles*/
			 char *pfades,  /*heisst im move_ctrl() "pfade"              */
			 struct TimeCell*timer,  /*heisst im move_ctrl() "time_ptr"           */
			 int step,
			 int xcell,
			 int ycell,
			 int serial_number){
/* ----------------------------------------------------------------------
   Variablendefinition
-------------------------------------------------------------------------*/
  char sFN[DEF_STR_LEN];		   /* Filename */   
  char sFN_NTrv[DEF_STR_LEN];		   /* Filename part */   
  FILE *fMOut;
  char sElementName[] ="    ";	   /* Fixed length Element name */

  int i;
/* ----------------------------------------------------------------------
   endof Variablendefinition
-------------------------------------------------------------------------*/

             
  /*SDcmt: Open File for writing (w+)*/
  sprintf (sFN_NTrv,"particle_infos%i",step);
  if(pfades==NULL)
    sprintf (sFN,"%s%s.m", "./",sFN_NTrv);
  else
    sprintf (sFN,"%s%s.m", pfades,sFN_NTrv);

  if ( (fMOut=fopen(sFN,"r")) == NULL )
    if ( (fMOut= MFOpen (sFN_NTrv, pfades, -1)) == NULL){
      printf("matlab_particle_out(...): ERROR -  Could not open file!\n");
      return;
    }	 
    else{
      fprintf(fMOut,"%%WRITING PARTICLE INFORMATIONS: \n\n");
      fprintf(fMOut,"%%location (...,:)  %% Ortskoordinaten x,y,z\n");
      fprintf(fMOut,"%%azimuth  (...  )  %% 0...2pi, location[3] == phi - azimut\n");  
      fprintf(fMOut,"%%polar    (...  )  %% 0...pi/2, location[4] == theta - polarwinkel\n");  
      fprintf(fMOut,"%%dz       (...  )  %% normal distance from surface \n");  
      fprintf(fMOut,"%%time     (...  )  %% laufzeit \n");  
      fprintf(fMOut,"%%v_abs    (...  )  %% v[0] - Betrag der geschwindigkeit in cm/s \n");  
      fprintf(fMOut,"%%velocity (...,:)  %% v[1]...[3] - x,y,z komponenten der geschw.\n");
      fprintf(fMOut,"%%v_azimuth(...  )  %% 0...2pi, ,v[5] - winkel relativ zur ...richtung\n");  
      fprintf(fMOut,"%%v_polar  (...  )  %% 0...pi/2, v[4] - winkel relativ zur ...richtung\n");  
      fprintf(fMOut,"%%energy   (...  )  %% Energy of particle \n");  
      fprintf(fMOut,"%%Ti       (...  )  %% Temperatur des Teilchens \n");      
      fprintf(fMOut,"%%charge   (...  )  %% Charge of particle\n");  
      fprintf(fMOut,"%%mass     (...  )  %% masse des Teilchens\n");  
      fprintf(fMOut,"%%z        (...  )  %% Ordnungszahl\n");     
      fprintf(fMOut,"%%elem_name     (...,:) %% Name des Elementes(SI,C,HE,H,D,T,...)\n");  
      fprintf(fMOut,"%%elem_type_phys(...,:) %% Type of species - redeposited or substrate or puffed\n");  
      fprintf(fMOut,"%%atoms    (...  )  %% # der atome in quasiteilchen\n");
    }
  else{ 
    fclose(fMOut);
    if ((fMOut=fopen(sFN,"at")) == NULL){
      printf("matlab_particle_out(...): ERROR -  Could not re-open file!\n");
      return;      
    }
  }	      
  
  fprintf(fMOut,"\n\n%% **** PARTICLE INFORMATIONS FOR CELL (%d,%d):\n\n",xcell+1,ycell+1);
  fprintf(fMOut,"\n\n%% **** Serial number of particles: %d\n\n",serial_number+1);
  for (i=0;i<loops;i++){
    fprintf(fMOut,"\n%% ** PARTICLE NUMBER %d (of %d)\n",i+1,loops);
    fprintf(fMOut,"location (%d,%d,%d,%i,:) = [ %.2e, %.2e, %.2e];\n",xcell+1,ycell+1,serial_number+1,i+1,particles[i].location[0],particles[i].location[1],particles[i].location[2]);
    fprintf(fMOut,"azimuth  (%d,%d,%d,%i  ) = %.2e;\n",xcell+1,ycell+1,serial_number+1,i+1,particles[i].location[3]);  
    fprintf(fMOut,"polar    (%d,%d,%d,%i  ) = %.2e;\n",xcell+1,ycell+1,serial_number+1,i+1,particles[i].location[4]);  
    fprintf(fMOut,"dz       (%d,%d,%d,%i  ) = %.2e;\n",xcell+1,ycell+1,serial_number+1,i+1,particles[i].dz);  
    fprintf(fMOut,"time     (%d,%d,%d,%i  ) = %.2e;\n",xcell+1,ycell+1,serial_number+1,i+1,particles[i].time);  
    fprintf(fMOut,"v_abs    (%d,%d,%d,%i  ) = %.2e; \n",xcell+1,ycell+1,serial_number+1,i+1,particles[i].v[0]);  
    fprintf(fMOut,"velocity (%d,%d,%d,%i,:) = [ %.2e, %.2e, %.2e];\n",xcell+1,ycell+1,serial_number+1,i+1,particles[i].v[1],particles[i].v[2],particles[i].v[3]);
    fprintf(fMOut,"v_azimuth(%d,%d,%d,%i  ) = %.2e;\n",xcell+1,ycell+1,serial_number+1,i+1,particles[i].v[5]);  
    fprintf(fMOut,"v_polar  (%d,%d,%d,%i  ) = %.2e;\n",xcell+1,ycell+1,serial_number+1,i+1,particles[i].v[4]);  
    fprintf(fMOut,"energy   (%d,%d,%d,%i  ) = %.2e;\n",xcell+1,ycell+1,serial_number+1,i+1,particles[i].energy);  
    fprintf(fMOut,"Ti       (%d,%d,%d,%i  ) = %.2e;\n",xcell+1,ycell+1,serial_number+1,i+1,particles[i].Ti);      
    fprintf(fMOut,"charge   (%d,%d,%d,%i  ) = %d;\n",xcell+1,ycell+1,serial_number+1,i+1,particles[i].charge);  
    fprintf(fMOut,"mass     (%d,%d,%d,%i  ) = %.2e;\n",xcell+1,ycell+1,serial_number+1,i+1,particles[i].mass);  
    fprintf(fMOut,"z        (%d,%d,%d,%i  ) = %d;\n",xcell+1,ycell+1,serial_number+1,i+1,particles[i].z);     
    strncpy (sElementName,"     ",4);strncpy (sElementName,particles[i].elem_name,strlen (particles[i].elem_name));		
    fprintf(fMOut,"elem_name     (%d,%d,%d,%i,:) = '%s';\n",xcell+1,ycell+1,serial_number+1,i+1,sElementName);  
    strncpy (sElementName,"     ",4);strncpy (sElementName,particles[i].elem_type_phys,strlen (particles[i].elem_type_phys));		
    fprintf(fMOut,"elem_type_phys(%d,%d,%d,%i,:) = '%s';\n",xcell+1,ycell+1,serial_number+1,i+1,sElementName);  
    fprintf(fMOut,"atoms    (%d,%d,%d,%i  ) = %.2e;\n",xcell+1,ycell+1,serial_number+1,i+1,particles[i].atoms);
  }
  fclose(fMOut);
} 


void CylProf(char *sName,
	     char cAx, 
		 double x1, double y1, double z1, 
		 double dR,
		 double dL,
		 int iNP,
		 struct Path          *pfade,
		 struct VolumeAll  *v_all,
		 int               iFMC,
	     struct V_TraceData  *VTr,
		 int    iTStep
	   )
{
   /* A*x + B*y + z = C    Profile axis */
   /* BCmt: later . . .           double A, B, C; */

   double dx        = v_all->diff[0];
   double dy        = v_all->diff[1];
   double dz        = v_all->diff[2];
   double dR2= dR*dR;
   double dX1       = v_all->adRange[0];
   double dY1       = v_all->adRange[1]; 
   double dZ1       = 0;
   int iNPx         = v_all->anz[0];    /* points in toroidal dir. */
   int iNPy         = v_all->anz[1];    /* points in poloidal dir. */
   int iNPz         = v_all->anz[2];
   int iNEl         = v_all->anz[3];    /* number of traced species */
   double dDens3D   = v_all->dDensFac;
                /* 1000. / (dx * dy * dz * time_c->step) */

   int ix,iy,iz;
   double *dBufNe, *dBufTe, *dLp;
   int *piNCl;
   double dxL;
   int indL, iCEl;
   double *dBufNp[MXN_TRACED_EL], *dBufSp[MXN_TRACED_EL];
   double dRC, dXC, dYC, dZC;
   
   char sVarN [DEF_STR_LEN], sSpNm[DEF_STR_LEN];
   FILE *fMOut;

  /* ---- Opening file, writing trivial things . . . ---- */
  if (iFMC==0)
    sprintf (sVarN,"Prf_%s", sName);
  else
    sprintf (sVarN,"Prf_%s_fm%i", sName, iFMC);
  fMOut= MFOpen (sVarN, pfade->RESULT, iTStep);
  fprintf(fMOut,"PrfName = '%s';\n%%\n",sName);

   dLp   =(double *) calloc (iNP,sizeof(double));
   piNCl =(int    *) calloc (iNP,sizeof(int));
   dBufNe=(double *) calloc (iNP,sizeof(double));
   dBufTe=(double *) calloc (iNP,sizeof(double));
   for (iCEl = 0; iCEl < iNEl; iCEl++){
     dBufNp[iCEl]=(double *) calloc (iNP,sizeof(double));
	 dBufSp[iCEl]=(double *) calloc (iNP,sizeof(double));
   }

   if (cAx == 'z')  dxL = dL/iNP;                /* axial profile  */
   if (cAx == 'Z')  dxL = dR/iNP;                /* radial profile */

   for (ix=0; ix<iNP; ix++)
   dLp[ix] = dxL*ix;

   if ((cAx == 'z')||(cAx == 'Z')){

     for (ix=0; ix<iNPx; ix++)
	 for (iy=0; iy<iNPy; iy++){
	     dXC       = v_all->b[ix][iy].edge[0]     + dx/2.;
         dYC       = v_all->b[ix][iy].edge[1]     + dy/2.;
        
		 dRC = sqrt(dXC*dXC + dYC*dYC);
	   for (iz=0; iz<iNPz; iz++){
		 dZC       = v_all->b[ix][iy].vcell[iz].z + dz/2. - z1;

		 if ( (dRC < dR) && (dZC > 0) && (dZC < dL) ){
		   if (cAx == 'z')                  /* axial profile  */
		     indL = (int)floor(dZC / dxL);
		   if (cAx == 'Z')                  /* radial profile */
		     indL = (int)floor(dRC / dxL);

		   piNCl [indL]+=1;
	       dBufTe[indL]+= v_all->b[ix][iy].vcell[iz].dTe;
	       dBufNe[indL]+= v_all->b[ix][iy].vcell[iz].dNe;
           for (iCEl = 0; iCEl < iNEl; iCEl++){
             dBufNp[iCEl][indL]+= v_all->b[ix][iy].vcell[iz].atomps[iCEl]*dDens3D;
	         dBufSp[iCEl][indL]+= v_all->b[ix][iy].vcell[iz].dAtomEmiss[iCEl];
		   }
		 }
	   }
	 }	 

   }else{       /* Axis type */
	  printf ("Error (CylProf).\n");
	  exit(0);
   }

   /* Averaging . . . */
   for (ix=0; ix<iNP; ix++){
	 if (piNCl[ix]>0){
       dBufTe[ix]/=piNCl[ix];
	   dBufNe[ix]/=piNCl[ix];
	   for (iCEl = 0; iCEl < iNEl; iCEl++){
         dBufNp[iCEl][ix]/=piNCl[ix];
	     dBufSp[iCEl][ix]/=piNCl[ix];
	   }
	 }
   }
   
   fprintf(fMOut,"%%\n");
   MlMatrOut (fMOut, "RpX", dLp   , iNP, 1, 'd');
   MlMatrOut (fMOut, "NVC", piNCl , iNP, 1, 'i');
#ifdef PISCES
   if (cAx == 'z'){         /* axial profile  */
	 fprintf(fMOut,"XAxNm = 'Z';\n");
	 fprintf(fMOut,"XAxUn = 'mm';\n");
   };               
   if (cAx == 'Z'){         /* radial profile */
	 fprintf(fMOut,"XAxNm = 'R';\n");
	 fprintf(fMOut,"XAxUn = 'mm';\n");
   };
#endif
   fprintf(fMOut,"%%\n");

   /* sprintf(sVarN,"(%i,:,:)",iz+1);*/
   fprintf(fMOut,"D.N    = %i;\n%%\n"      , 2+iNEl*2);

   fprintf(fMOut,"D.D{1}.sLab = '%s';\n"    , "n_e");
   MlMatrOut (fMOut, "D.D{1}.PX", dLp   , iNP, 1, 'd');
   MlMatrOut (fMOut, "D.D{1}.PY", dBufNe, iNP, 1, 'd');
   fprintf (fMOut, "%%\n%%\n");

   fprintf(fMOut,"D.D{2}.sLab = '%s';\n"    , "T_e");
   MlMatrOut (fMOut, "D.D{2}.PX", dLp   , iNP, 1, 'd');
   MlMatrOut (fMOut, "D.D{2}.PY", dBufTe, iNP, 1, 'd');
   fprintf (fMOut, "%%\n%%\n");

   for (iCEl = 0; iCEl < iNEl; iCEl++){
     if (v_all->charge[iCEl]>0)
       sprintf (sSpNm, "%s^{+%i}", v_all->name[iCEl], v_all->charge[iCEl]); 
	 else
	   sprintf (sSpNm, "%s^0", v_all->name[iCEl]); 

	 /* Density of traced species . . . */
	 fprintf(fMOut,"D.D{%i}.sLab = 'n(%s)';\n", 3+2*iCEl, sSpNm);
	 sprintf(sVarN,"D.D{%i}.PX",3+2*iCEl);
	 MlMatrOut (fMOut, sVarN, dLp         , iNP, 1, 'd');
	 sprintf(sVarN,"D.D{%i}.PY",3+2*iCEl);
	 MlMatrOut (fMOut, sVarN, dBufNp[iCEl], iNP, 1, 'd');
	 fprintf (fMOut, "%%\n%%\n");

	 /* Emission of traced species . . . */
	 if (VTr->iDataTp[iCEl] > 2)
	   fprintf(fMOut,"D.D{%i}.WL = %i;\n", 4+2*iCEl, (int)floor(VTr->dWLA[iCEl]));
	 fprintf(fMOut,"D.D{%i}.sLab = 'I(%s)';\n", 4+2*iCEl, sSpNm);
     sprintf(sVarN,"D.D{%i}.PX",4+2*iCEl);
	 MlMatrOut (fMOut, sVarN, dLp         , iNP, 1, 'd');
	 sprintf(sVarN,"D.D{%i}.PY",4+2*iCEl);
	 MlMatrOut (fMOut, sVarN, dBufSp[iCEl], iNP, 1, 'd');
	 fprintf (fMOut, "%%\n%%\n");
   }

   fprintf (fMOut, "%%\n%% ------------- End of ERO print --------------- \n");

   fclose (fMOut);
   free (dLp);
   free (dBufNe);
   free (dBufTe);
   free(piNCl);
   for (iCEl = 0; iCEl < iNEl; iCEl++){
     free(dBufNp[iCEl]);
	 free(dBufSp[iCEl]);
   }
}

/*----------------------------------------------------------------
   Prints out Ionization rates
-----------------------------------------------------------------*/
void MatlIRatesOut( 
	     struct VolumeAll    *v_all,
#ifdef OLD_IZ_MODEL
	     struct Ion		     *aIon,           /* Ionization potentials, etc.   */
#else
	     struct IZ_AtData    *IzData,
		 struct IZ_AtData    *RecData,
#ifdef METAST
		 struct METASTDATA *pMSD,
#endif
#endif /* OLD_IZ_MODEL */
	     struct V_TraceData  *VTr,
	     struct Path         *pfade,
         struct nTdata       *nT_ptr,
		 int    iTStep,
		 int    iMxCh
	   )
{
   const int iNP_ = 10;    /* Number of points to +- directions from 
						      cental value at LCFS */
   const int iNP = iNP_*2+1;
   const double dRS=1.5;   /* Factor between point by Ne, Te */
   double dCRS =1.0;
   double *dT, *dN, *dBufIzRt, *dBufRecRt, *dBufERt;
   int iNEl = v_all->anz[3];  /* Number of spectr. traced species */
   int i, iCNT, iCNN, iCEl;   /* Counters for ?, temper., dens., element */
#ifdef METAST
   int iMS;
#endif
   
   struct ParticleStruct Part;

   struct Spezies  spez;    /* Just to find max. charge . . . */

   FILE *fMOut;
   char sVarN[DEF_STR_LEN];
   
   printf("Writing Ionization rates . . .\n");

   dT    =(double *) calloc (iNP,sizeof(double));
   dN    =(double *) calloc (iNP,sizeof(double));
   dBufIzRt=(double *) calloc (iNP*iNP,sizeof(double));
   dBufRecRt=(double *) calloc (iNP*iNP,sizeof(double));
   dBufERt=(double *) calloc (iNP*iNP,sizeof(double));

   dCRS=1.0;
   dT[iNP_]=nT_ptr->Te;
   dN[iNP_]=nT_ptr->ne;
   if (dN[iNP_]<0) dN[iNP_] = -dN[iNP_];
   for (i=0; i<iNP_; i++){
     dCRS*=dRS;
     dT[iNP_-i-1]=dT[iNP_]/dCRS;
	 dT[iNP_+i+1]=dT[iNP_]*dCRS;
	 dN[iNP_-i-1]=dN[iNP_]/dCRS;
	 dN[iNP_+i+1]=dN[iNP_]*dCRS;
   }

   fMOut= MFOpen ("IzRates", pfade->RESULT, iTStep);
   fprintf(fMOut,"iNP = %i;\n"    ,iNP);
   fprintf(fMOut,"dRS = %.3f;\n"  ,dRS);
   fprintf(fMOut,"iNEl = %i;\n%%\n",iNEl);

   MlMatrOut (fMOut, "dTe", dT, iNP, 1, 'd');
   MlMatrOut (fMOut, "dNe", dN, iNP, 1, 'd');

#ifndef OLD_IZ_MODEL
   /* --- ionization and recombination statistic ----------*/
   fprintf(fMOut,"\n\n%% =============== Ionization/recombination statisctics ============== \n");
   for (iCEl=0; iCEl<IzData->iNEl; iCEl++){
	 fprintf(fMOut,"\n\n%% ------------------------ Specie #%i '%s' -------------------------- \n",
		   iCEl+1, IzData->sElNm[iCEl]);  
	 fprintf(fMOut,"IRSt_ElNm{%i} = '%s';\n",iCEl+1, IzData->sElNm[iCEl]);
	 sprintf(sVarN,"IRSt_Iz(%i,:)",iCEl+1);
     MlMatrOut(fMOut, sVarN,  IzData->dHistory[iCEl], MAX_IZ_ST, 1, 'd');
	 sprintf(sVarN,"IRSt_Rec(%i,:)",iCEl+1);
     MlMatrOut(fMOut, sVarN, RecData->dHistory[iCEl], MAX_IZ_ST, 1, 'd');
   }
   fprintf(fMOut,   "%% =================================================================== \n%%\n%%\n");
   /* --- ionization and recombination statistic (end) ----*/
#endif /* OLD_IZ_MODEL */
 
   for (iCEl=0; iCEl<iNEl; iCEl++){ 
	 Part.charge = v_all->charge[iCEl];
	 fprintf(fMOut,"\n%%\n%% =============== Specie #%i '%s' ============== \n",
		   iCEl, v_all->name[iCEl]);  

	 fprintf(fMOut,"SpNm{%i} = '%s';\n", iCEl+1, v_all->name[iCEl]);

	 if (Part.charge>1)
	   fprintf(fMOut,"SpNmCh{%i} = '%s^{+%i}';\n", iCEl+1, v_all->name[iCEl], Part.charge);
	 else if (Part.charge==1)
       fprintf(fMOut,"SpNmCh{%i} = '%s^{+}';\n", iCEl+1, v_all->name[iCEl]);
	 else 
	   fprintf(fMOut,"SpNmCh{%i} = '%s^{0}';\n", iCEl+1, v_all->name[iCEl]);
         
	 fprintf(fMOut,"SpCh{%i} = %i;\n", iCEl+1, Part.charge);

#ifdef PISCES
	 if (strncmp (v_all->name[iCEl], "Be", 2) == 0) v_all->name[iCEl][0] = 'b';
#endif

#ifdef OLD_IZ_MODEL
	 strncpy(Part.elem_name, v_all->name[iCEl], 2);
     Part.elem_name[2]='\0';
	 l_ratecoeff(aIon,     /* structure where the coeff will be written      */
		     &Part,	       /* contains info on the name of the element            */
		     iMxCh,        /* max. ionisation state considered                    */
					       /* defined above as CODE specific def                  */
		     pfade); 
#endif /* OLD_IZ_MODEL */
     
	 strncpy(spez.name, v_all->name[iCEl], 2);
	 spez.name[2]='\0';
	 ElemNo (&spez);

     Part.elem_name[2]='\0';  
     for (iCNT=0; iCNT<iNP; iCNT++)
	   for (iCNN=0; iCNN<iNP; iCNN++){  

         if (spez.z <= v_all->charge[iCEl]){      /* All electrons gone . . . */
			 dBufIzRt[iCNT+iNP*iCNN]= -1.;
			 break;
		 }

#ifdef OLD_IZ_MODEL
	     dBufIzRt[iCNT+iNP*iCNN]= rates(aIon, &Part, dT[iCNT]);
#else
		 dBufIzRt[iCNT+iNP*iCNN]= dIonizRt(IzData, v_all->name[iCEl], v_all->charge[iCEl], dT[iCNT], dN[iCNN]);
         if (RecData->iNEl>0)
		    dBufRecRt[iCNT+iNP*iCNN]= dIonizRt(RecData, v_all->name[iCEl], v_all->charge[iCEl], dT[iCNT], dN[iCNN]);
		 else 
			dBufRecRt[iCNT+iNP*iCNN]=-1.;
#endif /* OLD_IZ_MODEL */
		 
		 dBufERt[iCNT+iNP*iCNN]= 4*M_PI*emission(VTr,dT[iCNT], dN[iCNN],iCEl);
 	   }
      
	 sprintf(sVarN,"RtIz(%i,:,:)",iCEl+1);
	 MlMatrOut (fMOut, sVarN, dBufIzRt, iNP, iNP, 'd');
	 sprintf(sVarN,"RtRec(%i,:,:)",iCEl+1);
	 MlMatrOut (fMOut, sVarN, dBufRecRt, iNP, iNP, 'd');
	 sprintf(sVarN,"RtEm(%i,:,:)",iCEl+1);
	 MlMatrOut (fMOut, sVarN, dBufERt, iNP, iNP, 'd');

#ifdef METAST
    if (v_all->iMS[iCEl]!=MSUnres)
    if (Part.charge == 0)
    if ( strcmp(v_all->name[iCEl], "be")==0 ){
       /* Metastable transitions . . . */
       for (iMS=0; iMS<pMSD->iNTr; iMS++){
       for (iCNT=0; iCNT<iNP; iCNT++)
	   for (iCNN=0; iCNN<iNP; iCNN++){  
		 dBufIzRt[iCNT+iNP*iCNN] 
		    = dMetaStRt(pMSD->TrD[iMS], dT[iCNT], dN[iCNN]);
 	   }
	   fprintf(fMOut,"MS_RtNm{%i,%i} = '%s';\n", iCEl+1, iMS+1, pMSD->TrNm[iMS]);
	   sprintf(sVarN,"MS_RT(%i,%i,:,:)", iCEl+1, iMS+1);
	   MlMatrOut (fMOut, sVarN, dBufIzRt, iNP, iNP, 'd');
	 }/* iMS */
       /* Metastable emission. . .  */
       for (iMS=0; iMS<pMSD->iNLn*N_MSLEV; iMS++){
       for (iCNT=0; iCNT<iNP; iCNT++)
	   for (iCNN=0; iCNN<iNP; iCNN++){  
		 dBufIzRt[iCNT+iNP*iCNN] 
		    = dMetaStRt(pMSD->SpecD[iMS], dT[iCNT], dN[iCNN]);
 	   }
	   fprintf(fMOut,"MS_WL(%i) = %i;\n", 
		             iMS/N_MSLEV+1, pMSD->iaWL[iMS/N_MSLEV]);
	   sprintf(sVarN,"MS_PEC(%i,%i,:,:)", iMS/N_MSLEV+1, iMS%N_MSLEV+1);
	   MlMatrOut (fMOut, sVarN, dBufIzRt, iNP, iNP, 'd');
	 }/* iMS */
    } /* if be */
#endif
   } /* iCEl */

   free(dN);
   free(dT);
   free(dBufIzRt);
   free(dBufRecRt);
   free(dBufERt);

   fclose (fMOut);

   printf("Writing of ionization rates complete.\n");
}

#ifdef UGEOMETRY
#ifdef PRINT_BFIELD
void MatlBFOut( 
	     struct Path          *pfade,
		 double		          *dBound,
	     struct  Numeric           *num,
	     struct Ort                *ort,  /* contains limiter info                        */
	     struct Electric       *e_param,  /* therefore we are here                        */
	     struct BField               *b,  
	     struct Magnetic      *magnetic,
	     struct Plasma           *p_ptr,
             struct nTdata          *nT_ptr,
	     struct PConst      *plas_const,
	     int                     surf
	   )
{
   const int iNPx=51;
   const int iNPy=16;
   const int iNPz=21;
   const double dRS=1;
   double dx,dy,dz, dxL,dyL,dzL, dX1, dY1, dZ1;
   double dXL, dXR, dYL, dYR;
   int ix,iy,iz;
   double *dBufDzpro, *dBufDzcontra, *dBufPro, *dBufContra, *dX, *dY, *dZ;
   double *pos;
   double dZmin=0, dZmax;

   FILE *fMOut;
   char sVarN[DEF_STR_LEN];

   double progress;

   struct EroBfieldInputList eroList;
   struct EroBfieldOutputData eroData;

   printf("Writing B_field . . .\n");

   dXL = dBound[0];
   dXR = dBound[1];
   dxL=dXR-dXL;
   dx=dxL/(iNPx-1)*dRS;
   dX1=dXL - dxL*(dRS-1)/2;

   dYL = dBound[2];
   dYR = dBound[3];
   dyL=dYR - dYL;
   dy=dyL/(iNPy-1)*dRS;
   dY1=dYL - dyL*(dRS-1)/2;
   
   dZmin = dBound[4];
   dZmax = dBound[5];
   dzL=dZmax-dZmin;
   dz=dzL/(iNPz-1)*dRS;
   dZ1=dZmin - dzL*(dRS-1)/2;

   fMOut= fopen("BField.m", "wt");

   fprintf(fMOut,"iNPx = %i;\n%%\n",iNPx);
   fprintf(fMOut,"iNPy = %i;\n%%\n",iNPy);
   fprintf(fMOut,"iNPz = %i;\n%%\n",iNPz);
   fprintf(fMOut,"iSType = %i;\n%%\n",surf);

   fprintf(fMOut,"dx = %.2e;\n",dx);
   fprintf(fMOut,"dy = %.2e;\n",dy);
   fprintf(fMOut,"dz = %.2e;\n",dz);

   pos = (double *) calloc(3, sizeof(double));
   dBufDzpro = (double *) calloc (iNPx*iNPy*iNPz,sizeof(double));
   dBufDzcontra = (double *) calloc (iNPx*iNPy*iNPz,sizeof(double));
   dX=(double *) calloc (iNPx,sizeof(double));
   dY=(double *) calloc (iNPy,sizeof(double));
   dZ=(double *) calloc (iNPz,sizeof(double));

   for (ix=0; ix<iNPx; ix++)
	 dX[ix]=dX1+ix*dx;
   for (iy=0; iy<iNPy; iy++)
	 dY[iy]=dY1+iy*dy;
   for (iz=0; iz<iNPz; iz++)
	 dZ[iz]=dZ1+iz*dz;

   MlMatrOut (fMOut, "dX", dX, iNPx, 1, 'd');
   MlMatrOut (fMOut, "dY", dY, iNPy, 1, 'd');
   MlMatrOut (fMOut, "dZ", dZ, iNPz, 1, 'd');

   fclose (fMOut);

   printf("output progress . . . 00%%");

   #ifdef PRINT_BTRACE
   fBfield_init("BfieldTrace.m");
   #endif

   eroList.minBound[0] = num->UGbbox[0];
   eroList.minBound[1] = num->UGbbox[2];
   eroList.minBound[2] = num->UGbbox[4];
   eroList.maxBound[0] = num->UGbbox[1];
   eroList.maxBound[1] = num->UGbbox[3];
   eroList.maxBound[2] = num->UGbbox[5];
   eroList.dxAuto = num->BdxAuto;
   eroList.dxOut = num->BdxOut;
   eroList.dxIn = num->BdxIn;
   #ifdef ITER_BM
   eroList.Bfunc = &BM_B_field;
   #else
   eroList.Bfunc = &B_field;
   #endif
   eroList.pB = b;
   eroList.pMagnetic = magnetic;
   eroList.printBTrace = true;	/* to output the trace of the field lines if #ifdef PRINT_BTRACE*/

   for (iz=0; iz<iNPz; iz++){

	 progress = (double)iz/(iNPz-1)*100;
	 printf("\r\r\routput progress . . . %03d%%", (int)progress);

     for (ix=0; ix<iNPx; ix++)
	 for (iy=0; iy<iNPy; iy++){

       pos[0]=dX[ix];
	   pos[1]=dY[iy];
       pos[2]=dZ[iz];

	   findBfieldDistance(pos, &eroList, &eroData);
     
	   if (eroData.fwdFlag)
		dBufDzpro[ix+iNPx*(iy+iNPy*iz)] = eroData.fwdDist;
	   else
		dBufDzpro[ix+iNPx*(iy+iNPy*iz)] = -1.;
	   if (eroData.bwdFlag)
		dBufDzcontra[ix+iNPx*(iy+iNPy*iz)] = eroData.bwdDist;
	   else
		dBufDzcontra[ix+iNPx*(iy+iNPy*iz)] = -1.;
	 }
   }

   #ifdef PRINT_BTRACE
   fBfield_close();
   #endif
   
   dBufPro = (double *) calloc (iNPx*iNPy,sizeof(double));
   dBufContra = (double *) calloc (iNPx*iNPy,sizeof(double));

   fMOut= fopen("BField.m", "a+t");
   for (iz=0; iz<iNPz; iz++){
	 for (ix=0; ix<iNPx; ix++)
	 for (iy=0; iy<iNPy; iy++){
		 dBufPro[ix+iNPx*iy] = dBufDzpro[ix+iNPx*(iy+iNPy*iz)];
		 dBufContra[ix+iNPx*iy] = dBufDzcontra[ix+iNPx*(iy+iNPy*iz)];
	 }
     sprintf(sVarN,"dB_dz_pro(%i,:,:)",iz+1);
	 MlMatrOut (fMOut, sVarN, dBufPro, iNPx, iNPy, 'd');
	 sprintf(sVarN,"dB_dz_contra(%i,:,:)",iz+1);
	 MlMatrOut (fMOut, sVarN, dBufContra, iNPx, iNPy, 'd');
   }
   fclose (fMOut);

   free(pos);
   free(dBufPro);
   free(dBufContra);
   free(dBufDzpro);
   free(dBufDzcontra);
   free(dX);
   free(dY);
   free(dZ);
   printf("\nWriting of B_field complete.\n");
}
#endif
#endif

