#ifndef _GEOM_
#define _GEOM_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "Boundary.h"
#include "Particles.h"
#include "Surfaces.h"
#include "surfaceModel.h"
#include <cmath>

CUDA_CALLABLE_MEMBER_DEVICE
float findT(float x0, float x1, float y0, float y1, float intersectionx) {

  float a, b, c, a1, a2, t=0, discriminant, realPart, imaginaryPart;
  a = (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0);
  b = 2.0 * x0 * (x1 - x0) + 2.0 * y0 * (y1 - y0);
  c = x0 * x0 + y0 * y0 - intersectionx * intersectionx;
  discriminant = b * b - 4 * a * c;

  if (discriminant > 0) {
    a1 = (-b + std::sqrt(discriminant)) / (2 * a);
    a2 = (-b - std::sqrt(discriminant)) / (2 * a);
    //std::cout << "Roots are real and different." << std::endl;
    //std::cout << "a1 = " << a1 << std::endl;
    //std::cout << "a2 = " << a2 << std::endl;
    t = std::min(std::abs(a1),std::abs(a2));
  }

  else if (discriminant == 0) {
    // cout << "Roots are real and same." << endl;
    a1 = (-b + std::sqrt(discriminant)) / (2 * a);
    // cout << "a1 = a2 =" << a1 << endl;
  }

  else {
    realPart = -b / (2 * a);
    imaginaryPart = std::sqrt(-discriminant) / (2 * a);
    // cout << "Roots are complex and different."  << endl;
    // cout << "a1 = " << realPart << "+" << imaginaryPart << "i" << endl;
    // cout << "a2 = " << realPart << "-" << imaginaryPart << "i" << endl;
  }

  return t;
}
struct geometry_check {
  Particles *particlesPointer;
  const int nLines;
  Boundary *boundaryVector;
  Surfaces *surfaces;
  float dt;
  // int& tt;
  int nHashes;
  int *nR_closeGeom;
  int *nY_closeGeom;
  int *nZ_closeGeom;
  int *n_closeGeomElements;
  float *closeGeomGridr;
  float *closeGeomGridy;
  float *closeGeomGridz;
  int *closeGeom;
  int nEdist;
  float E0dist;
  float Edist;
  int nAdist;
  float A0dist;
  float Adist;

  geometry_check(Particles *_particlesPointer, int _nLines,
                 Boundary *_boundaryVector, Surfaces *_surfaces, float _dt,
                 int _nHashes, int *_nR_closeGeom, int *_nY_closeGeom,
                 int *_nZ_closeGeom, int *_n_closeGeomElements,
                 float *_closeGeomGridr, float *_closeGeomGridy,
                 float *_closeGeomGridz, int *_closeGeom, int _nEdist,
                 float _E0dist, float _Edist, int _nAdist, float _A0dist,
                 float _Adist)
      :

        particlesPointer(_particlesPointer), nLines(_nLines),
        boundaryVector(_boundaryVector), surfaces(_surfaces), dt(_dt),
        nHashes(_nHashes), nR_closeGeom(_nR_closeGeom),
        nY_closeGeom(_nY_closeGeom), nZ_closeGeom(_nZ_closeGeom),
        n_closeGeomElements(_n_closeGeomElements),
        closeGeomGridr(_closeGeomGridr), closeGeomGridy(_closeGeomGridy),
        closeGeomGridz(_closeGeomGridz), closeGeom(_closeGeom), nEdist(_nEdist),
        E0dist(_E0dist), Edist(_Edist), nAdist(_nAdist), A0dist(_A0dist),
        Adist(_Adist) {}

  CUDA_CALLABLE_MEMBER_DEVICE
  void operator()(std::size_t indx) const {
    // std::cout << "geometry check particle x" << particlesPointer->x[indx] <<
    // particlesPointer->x[indx]previous <<std::endl; std::cout << "geometry
    // check particle y" << particlesPointer->y[indx] <<
    // particlesPointer->y[indx]previous <<std::endl; std::cout << "geometry
    // check particle z" << particlesPointer->z[indx] <<
    // particlesPointer->z[indx]previous <<std::endl; std::cout << "geometry
    // check particle hitwall" << p.hitWall <<std::endl;
    if (particlesPointer->hitWall[indx] == 0.0) {
      int hitSurface = 0;
      float x = particlesPointer->x[indx];
      float y = particlesPointer->y[indx];
      float z = particlesPointer->z[indx];
      float xprev = particlesPointer->xprevious[indx];
      float yprev = particlesPointer->yprevious[indx];
      float zprev = particlesPointer->zprevious[indx];
      float dpath =
          std::sqrt((x - xprev) * (x - xprev) + (y - yprev) * (y - yprev) +
                    (z - zprev) * (z - zprev));
#if FLUX_EA > 0
      float dEdist = (Edist - E0dist) / static_cast<float>(nEdist);
      float dAdist = (Adist - A0dist) / static_cast<float>(nAdist);
      int AdistInd = 0;
      int EdistInd = 0;
#endif
      float vxy[3] = {0.0f};
      float vtheta[3] = {0.0f};
#if USECYLSYMM > 0
      if (boundaryVector[nLines].periodic) // if periodic
      {
        float pi = 3.14159265;
        float theta =
            std::atan2(particlesPointer->y[indx], particlesPointer->x[indx]);
        float thetaPrev = std::atan2(particlesPointer->yprevious[indx],
                                 particlesPointer->xprevious[indx]);
        // float vtheta =
        // std::atan2(particlesPointer->vy[indx],particlesPointer->vx[indx]);
        float rprev = std::sqrt(particlesPointer->xprevious[indx] *
                               particlesPointer->xprevious[indx] +
                           particlesPointer->yprevious[indx] *
                               particlesPointer->yprevious[indx]);
        float r = std::sqrt(particlesPointer->x[indx] * particlesPointer->x[indx] +
                       particlesPointer->y[indx] * particlesPointer->y[indx]);
        float rHat[3] = {0.0f};
        float vr[3] = {0.0f};
        rHat[0] = particlesPointer->x[indx];
        rHat[1] = particlesPointer->y[indx];

        vectorNormalize(rHat, rHat);
        vxy[0] = particlesPointer->vx[indx];
        vxy[1] = particlesPointer->vy[indx];
        vectorScalarMult(vectorDotProduct(rHat, vxy), rHat, vr);
        float vrMag = vectorNorm(vr);
        vectorSubtract(vxy, vr, vtheta);
        float vthetaMag = vectorNorm(vtheta);
        float vx0 = 0.0;
        float vy0 = 0.0;
        if (theta <= boundaryVector[nLines].y1) {
          particlesPointer->xprevious[indx] =
              r * std::cos(boundaryVector[nLines].y2 + theta);
          particlesPointer->yprevious[indx] =
              r * std::sin(boundaryVector[nLines].y2 + theta);
          particlesPointer->x[indx] =
              rprev * std::cos(boundaryVector[nLines].y2 + theta);
          particlesPointer->y[indx] =
              rprev * std::sin(boundaryVector[nLines].y2 + theta);

          vx0 = vrMag * std::cos(boundaryVector[nLines].y2 + theta) -
                vthetaMag * std::sin(boundaryVector[nLines].y2 + theta);
          vy0 = vrMag * std::sin(boundaryVector[nLines].y2 + theta) +
                vthetaMag * std::cos(boundaryVector[nLines].y2 + theta);
          particlesPointer->vx[indx] = vx0;
          particlesPointer->vy[indx] = vy0;
        } else if (theta >= boundaryVector[nLines].y2) {
          particlesPointer->xprevious[indx] =
              rprev * std::cos(thetaPrev - boundaryVector[nLines].y2);
          particlesPointer->yprevious[indx] =
              rprev * std::sin(thetaPrev - boundaryVector[nLines].y2);
          particlesPointer->x[indx] =
              r * std::cos(theta - boundaryVector[nLines].y2);
          particlesPointer->y[indx] =
              r * std::sin(theta - boundaryVector[nLines].y2);

          vx0 = vrMag * std::cos(theta - boundaryVector[nLines].y2) -
                vthetaMag * std::sin(theta - boundaryVector[nLines].y2);
          vy0 = vrMag * std::sin(theta - boundaryVector[nLines].y2) +
                vthetaMag * std::cos(theta - boundaryVector[nLines].y2);
          particlesPointer->vx[indx] = vx0;
          particlesPointer->vy[indx] = vy0;
        }
      }
#else
      if (boundaryVector[nLines].periodic) {
        if (particlesPointer->y[indx] < boundaryVector[nLines].y1) {
          particlesPointer->y[indx] =
              boundaryVector[nLines].y2 -
              (boundaryVector[nLines].y1 - particlesPointer->y[indx]);
          particlesPointer->yprevious[indx] =
              boundaryVector[nLines].y2 -
              (boundaryVector[nLines].y1 - particlesPointer->y[indx]);

        } else if (particlesPointer->y[indx] > boundaryVector[nLines].y2) {
          particlesPointer->y[indx] =
              boundaryVector[nLines].y1 +
              (particlesPointer->y[indx] - boundaryVector[nLines].y2);
          particlesPointer->yprevious[indx] =
              boundaryVector[nLines].y1 +
              (particlesPointer->y[indx] - boundaryVector[nLines].y2);
        }
      }
#endif
#if USE3DTETGEOM > 0

      float a = 0.0;
      float b = 0.0;
      float c = 0.0;
      float d = 0.0;
      float plane_norm = 0.0;
      float pointToPlaneDistance0 = 0.0;
      float pointToPlaneDistance1 = 0.0;
      float signPoint0 = 0.0;
      float signPoint1 = 0.0;
      float t = 0.0;
      float A[3] = {0.0, 0.0, 0.0};
      float B[3] = {0.0, 0.0, 0.0};
      float C[3] = {0.0, 0.0, 0.0};
      float AB[3] = {0.0, 0.0, 0.0};
      float AC[3] = {0.0, 0.0, 0.0};
      float BC[3] = {0.0, 0.0, 0.0};
      float CA[3] = {0.0, 0.0, 0.0};
      float p[3] = {0.0, 0.0, 0.0};
      float Ap[3] = {0.0, 0.0, 0.0};
      float Bp[3] = {0.0, 0.0, 0.0};
      float Cp[3] = {0.0, 0.0, 0.0};
      float normalVector[3] = {0.0, 0.0, 0.0};
      float crossABAp[3] = {0.0, 0.0, 0.0};
      float crossBCBp[3] = {0.0, 0.0, 0.0};
      float crossCACp[3] = {0.0, 0.0, 0.0};
      float signDot0 = 0.0;
      float signDot1 = 0.0;
      float signDot2 = 0.0;
      float totalSigns = 0.0;
      int nBoundariesCrossed = 0;
      int boundariesCrossed[6] = {0, 0, 0, 0, 0, 0};
      /*
      if(particlesPointer->xprevious[indx] < 0 ||
      particlesPointer->yprevious[indx] < 0 ||
              particlesPointer->zprevious[indx]< 0)
      {
      std::cout << "pos " << particlesPointer->xprevious[indx] << " "
          << particlesPointer->yprevious[indx]
          << " " << particlesPointer->zprevious[indx]  << std::endl;
      }
       */

      float p0[3] = {particlesPointer->xprevious[indx],
                     particlesPointer->yprevious[indx],
                     particlesPointer->zprevious[indx]};
      float p1[3] = {particlesPointer->x[indx], particlesPointer->y[indx],
                     particlesPointer->z[indx]};
#if GEOM_HASH > 0
      // find which hash
      int nHash = 0;
      int rHashInd = 0;
      int yHashInd = 0;
      int zHashInd = 0;
      int rHashInd1 = 0;
      int yHashInd1 = 0;
      int zHashInd1 = 0;
      float r_position = particlesPointer->xprevious[indx];
      for (int i = 0; i < nHashes; i++) {
        rHashInd1 = nR_closeGeom[i] - 1;
        yHashInd1 = nY_closeGeom[i] - 1;
        zHashInd1 = nZ_closeGeom[i] - 1;
        if (i > 0)
          rHashInd = nR_closeGeom[i - 1];
        if (i > 0)
          yHashInd = nY_closeGeom[i - 1];
        if (i > 0)
          zHashInd = nZ_closeGeom[i - 1];
        if (i > 0)
          rHashInd1 = nR_closeGeom[i - 1] + nR_closeGeom[i] - 1;
        if (i > 0)
          yHashInd1 = nY_closeGeom[i - 1] + nY_closeGeom[i] - 1;
        if (i > 0)
          zHashInd1 = nZ_closeGeom[i - 1] + nZ_closeGeom[i] - 1;
        // std::cout << "rpos " <<rHashInd<< " " << rHashInd1 << " " <<
        // closeGeomGridr[rHashInd] << " "
        //          << closeGeomGridr[rHashInd1] << std::endl;
        // std::cout << "ypos " << closeGeomGridy[yHashInd] << " "
        //          << closeGeomGridy[yHashInd1] << std::endl;
        // std::cout << "zpos " << closeGeomGridz[zHashInd] << " "
        //         << closeGeomGridz[zHashInd1] << std::endl;
        if (r_position < closeGeomGridr[rHashInd1] &&
            r_position > closeGeomGridr[rHashInd] &&
            particlesPointer->yprevious[indx] < closeGeomGridy[yHashInd1] &&
            particlesPointer->yprevious[indx] > closeGeomGridy[yHashInd] &&
            particlesPointer->zprevious[indx] < closeGeomGridz[zHashInd1] &&
            particlesPointer->zprevious[indx] > closeGeomGridz[zHashInd]) {
          nHash = i;
        }
      }
      // std::cout << "nHash " << nHash << std::endl;
      rHashInd = 0;
      yHashInd = 0;
      zHashInd = 0;
      if (nHash > 0)
        rHashInd = nR_closeGeom[nHash - 1];
      if (nHash > 0)
        yHashInd = nY_closeGeom[nHash - 1];
      if (nHash > 0)
        zHashInd = nZ_closeGeom[nHash - 1];
      float dr = closeGeomGridr[rHashInd + 1] - closeGeomGridr[rHashInd];
      float dz = closeGeomGridz[zHashInd + 1] - closeGeomGridz[zHashInd];
      float dy = closeGeomGridy[yHashInd + 1] - closeGeomGridy[yHashInd];
      int rInd = std::floor((r_position - closeGeomGridr[rHashInd]) / dr + 0.5f);
      int zInd = std::floor(
          (particlesPointer->zprevious[indx] - closeGeomGridz[zHashInd]) / dz +
          0.5f);
      int i = 0;
      int yInd = std::floor(
          (particlesPointer->yprevious[indx] - closeGeomGridy[yHashInd]) / dy +
          0.5f);
      // std::cout << "rHashInd " << rHashInd << " " << yHashInd << " " <<
      // zHashInd << std::endl; std::cout << "dr dy dz " << dr << " " << dy << "
      // " << dz << std::endl; std::cout << "rind y z " << rInd << " " << yInd <<
      // " " << zInd << std::endl;
      if (rInd < 0 || yInd < 0 || zInd < 0) {
        rInd = 0;
        yInd = 0;
        zInd = 0;
#if USE_CUDA
#else
        // std::cout << "WARNING: particle outside of geometry hash range (low)"
        // << std::endl;
#endif
      } else if (rInd > nR_closeGeom[nHash] - 1 ||
                 yInd > nY_closeGeom[nHash] - 1 ||
                 zInd > nZ_closeGeom[nHash] - 1) {
        rInd = 0;
        yInd = 0;
        zInd = 0;
      }
      int buffIndx = 0;
      if (nHash > 0)
        buffIndx = nR_closeGeom[nHash - 1] * nY_closeGeom[nHash - 1] *
                   nZ_closeGeom[nHash - 1] * n_closeGeomElements[nHash - 1];
      // std::cout << "buff Index " << buffIndx << std::endl;
      for (int j = 0; j < n_closeGeomElements[nHash]; j++) {
        i = closeGeom[buffIndx +
                      zInd * nY_closeGeom[nHash] * nR_closeGeom[nHash] *
                          n_closeGeomElements[nHash] +
                      yInd * nR_closeGeom[nHash] * n_closeGeomElements[nHash] +
                      rInd * n_closeGeomElements[nHash] + j];
        // std::cout << "i's " << i << std::endl;
#else
      for (int i = 0; i < nLines; i++) {
#endif
        a = boundaryVector[i].a;
        b = boundaryVector[i].b;
        c = boundaryVector[i].c;
        d = boundaryVector[i].d;
        plane_norm = boundaryVector[i].plane_norm;
        pointToPlaneDistance0 =
            (a * p0[0] + b * p0[1] + c * p0[2] + d) / plane_norm;
        pointToPlaneDistance1 =
            (a * p1[0] + b * p1[1] + c * p1[2] + d) / plane_norm;
        // std::cout << "plane coeffs "<< i << " " << a << " " << b << " " << c
        // << " " << d << " " << plane_norm << std::endl; std::cout << "point to
        // plane dists "<< i << " " << pointToPlaneDistance0 << " " <<
        // pointToPlaneDistance1 << std::endl;
        signPoint0 = std::copysign(1.0, pointToPlaneDistance0);
        signPoint1 = std::copysign(1.0, pointToPlaneDistance1);

        if (signPoint0 != signPoint1) {
          t = -(a * p0[0] + b * p0[1] + c * p0[2] + d) /
              (a * (p1[0] - p0[0]) + b * (p1[1] - p0[1]) + c * (p1[2] - p0[2]));
          vectorAssign(p0[0] + t * (p1[0] - p0[0]), p0[1] + t * (p1[1] - p0[1]),
                       p0[2] + t * (p1[2] - p0[2]), p);
          // std::cout << " p " << p[0] << " " << p[1] << " " << p[2] <<
          // std::endl;
          vectorAssign(boundaryVector[i].x1, boundaryVector[i].y1,
                       boundaryVector[i].z1, A);
          vectorAssign(boundaryVector[i].x2, boundaryVector[i].y2,
                       boundaryVector[i].z2, B);
          vectorAssign(boundaryVector[i].x3, boundaryVector[i].y3,
                       boundaryVector[i].z3, C);

          vectorSubtract(B, A, AB);
          vectorSubtract(C, A, AC);
          vectorSubtract(C, B, BC);
          vectorSubtract(A, C, CA);

          vectorSubtract(p, A, Ap);
          vectorSubtract(p, B, Bp);
          vectorSubtract(p, C, Cp);

          vectorCrossProduct(AB, AC, normalVector);
          vectorCrossProduct(AB, Ap, crossABAp);
          vectorCrossProduct(BC, Bp, crossBCBp);
          vectorCrossProduct(CA, Cp, crossCACp);
          // std::cout << "AB " << AB[0] << " " << AB[1] << " " << AB[2] <<
          // std::endl; std::cout << "Ap " << Ap[0] << " " << Ap[1] << " " <<
          // Ap[2] << std::endl; std::cout << "BC " << BC[0] << " " << BC[1] << "
          // " << BC[2] << std::endl; std::cout << "Bp " << Bp[0] << " " << Bp[1]
          // << " " << Bp[2] << std::endl; std::cout << "CA " << CA[0] << " " <<
          // CA[1] << " " << CA[2] << std::endl; std::cout << "Cp " << Cp[0] << "
          // " << Cp[1] << " " << Cp[2] << std::endl;
          signDot0 =
              std::copysign(1.0, vectorDotProduct(crossABAp, normalVector));
          signDot1 =
              std::copysign(1.0, vectorDotProduct(crossBCBp, normalVector));
          signDot2 =
              std::copysign(1.0, vectorDotProduct(crossCACp, normalVector));
          totalSigns = 1.0 * std::abs(signDot0 + signDot1 + signDot2);
          // std::cout << "signdots " << signDot0 << " " << signDot1 << " " <<
          // signDot2 << " " << totalSigns << " " << (totalSigns
          // == 3.0)<<std::endl;

          // std::cout << "before loop totalSigns hitSurface " << totalSigns <<
          // " " << hitSurface << std::endl;
          hitSurface = 0;
          if (totalSigns == 3.0) {
            // std::cout << "in loop totalSigns hitSurface " << totalSigns << "
            // " << hitSurface << std::endl;
            hitSurface = 1;
          }
          if (vectorNorm(crossABAp) == 0.0 || vectorNorm(crossBCBp) == 0.0 ||
              vectorNorm(crossCACp) == 0.0) {
            hitSurface = 1;
          }
          // std::cout << "totalSigns hitSurface " << totalSigns << " " <<
          // hitSurface << std::endl;
          if (hitSurface == 1) {
            boundariesCrossed[nBoundariesCrossed] = i;
            // std::cout << "boundary crossed " << i << std::endl;
            nBoundariesCrossed++;
            particlesPointer->hitWall[indx] = 1.0;
            particlesPointer->xprevious[indx] = p[0];
            particlesPointer->yprevious[indx] = p[1];
            particlesPointer->zprevious[indx] = p[2];
            particlesPointer->x[indx] = p[0];
            particlesPointer->y[indx] = p[1];
            particlesPointer->z[indx] = p[2];
            particlesPointer->wallHit[indx] = i;
            //#if USESURFACEMODEL == 0
            //  #if USE_CUDA > 0
            //    atomicAdd(&boundaryVector[i].impacts,
            //    particlesPointer->weight[indx]);
            //  #else
            //    boundaryVector[i].impacts = boundaryVector[i].impacts +
            //    particlesPointer->weight[indx];
            //  #endif
            //#endif
            float E0 =
                0.5 * particlesPointer->amu[indx] * 1.66e-27 *
                (particlesPointer->vx[indx] * particlesPointer->vx[indx] +
                 particlesPointer->vy[indx] * particlesPointer->vy[indx] +
                 particlesPointer->vz[indx] * particlesPointer->vz[indx]) /
                1.602e-19;
            // std::cout << "Energy of particle that hit surface " << E0 <<
            // std::endl;
#if USE_CUDA > 0
            // float surfNormal[3] = {0.0f};
            // float partNormal[3] = {0.0f};
            // float partDotNormal = 0.0f;
            // partNormal[0] = particlesPointer->vx[indx];
            // partNormal[1] = particlesPointer->vy[indx];
            // partNormal[2] = particlesPointer->vz[indx];
            // getBoundaryNormal(boundaryVector,i,surfNormal,particlesPointer->x[indx],particlesPointer->y[indx]);
            // vectorNormalize(partNormal,partNormal);
            // partDotNormal = vectorDotProduct(partNormal,surfNormal);
            // float thetaImpact = acos(partDotNormal)*180.0/3.1415;
            // if(E0 < surfaces->E && E0> surfaces->E0)
            //{
            //    int tally_index = floor((E0-surfaces->E0)/surfaces->dE);
            //    if(thetaImpact > surfaces->A0 && thetaImpact < surfaces->A)
            //    {
            //        int aTally =
            //        floor((thetaImpact-surfaces->A0)/surfaces->dA);
            //        atomicAdd(&surfaces->energyDistribution[i*surfaces->nE*surfaces->nA+
            //        aTally*surfaces->nE + tally_index],
            //        particlesPointer->weight[indx]);
            //    }
            //}
#else
#endif
          }
        }
        if (nBoundariesCrossed == 0) {
          particlesPointer->xprevious[indx] = particlesPointer->x[indx];
          particlesPointer->yprevious[indx] = particlesPointer->y[indx];
          particlesPointer->zprevious[indx] = particlesPointer->z[indx];
        }
      }
#else // 2D geometry
#if USECYLSYMM > 0
      float pdim1 = std::sqrt(particlesPointer->x[indx] * particlesPointer->x[indx] +
                         particlesPointer->y[indx] * particlesPointer->y[indx]);
      float pdim1previous = std::sqrt(particlesPointer->xprevious[indx] *
                                     particlesPointer->xprevious[indx] +
                                 particlesPointer->yprevious[indx] *
                                     particlesPointer->yprevious[indx]);
      float theta0 = std::atan2(particlesPointer->yprevious[indx],
                            particlesPointer->xprevious[indx]);
      float theta1 =
          std::atan2(particlesPointer->y[indx], particlesPointer->x[indx]);
      float thetaNew = 0;
      float rNew = 0;
      float xNew = 0;
      float yNew = 0;
#else
      float pdim1 = particlesPointer->x[indx];
      float pdim1previous = particlesPointer->xprevious[indx];
#endif
      float particle_slope =
          (particlesPointer->z[indx] - particlesPointer->zprevious[indx]) /
          (pdim1 - pdim1previous);
      float particle_intercept =
          -particle_slope * pdim1 + particlesPointer->z[indx];
      float intersectionx[2] = {};
      // intersectionx = new float[nPoints];
      float intersectiony[2] = {};
      // intersectiony = new float[nPoints];
      float distances[2] = {};
      // distances = new float[nPoints];
      int intersectionIndices[2] = {};
      float tol_small = 1e-12f;
      float tol = 1e12f;
      int nIntersections = 0;
      float signPoint;
      float signPoint0;
      float signLine1;
      float signLine2;
      float minDist = 1e12f;
      int minDistInd = 0;
 //std::cout << "particle slope " << particle_slope << " " << particle_intercept
 //<< std::endl; std::cout << "r " << boundaryVector[0].x1 << " " <<
 //boundaryVector[0].x1 << " " << boundaryVector[0].slope_dzdx << std::endl;
 //std::cout << "r0 " << particlesPointer->x[indx] << " " <<
 //particlesPointer->y[indx] << " " <<
 //particlesPointer->z[indx]<< std::endl;
#if GEOM_HASH > 0
#if USECYLSYMM > 0
      float r_position = std::sqrt(particlesPointer->xprevious[indx] *
                                   particlesPointer->xprevious[indx] +
                               particlesPointer->yprevious[indx] *
                                   particlesPointer->yprevious[indx]);
#else
      float r_position = particlesPointer->xprevious[indx];
#endif
      float dr = closeGeomGridr[1] - closeGeomGridr[0];
      float dz = closeGeomGridz[1] - closeGeomGridz[0];
      int rInd = std::floor((r_position - closeGeomGridr[0]) / dr + 0.5f);
      int zInd = std::floor(
          (particlesPointer->zprevious[indx] - closeGeomGridz[0]) / dz + 0.5f);
      if (rInd < 0 || rInd >= nR_closeGeom[0])
        rInd = 0;
      if (zInd < 0 || zInd >= nZ_closeGeom[0])
        zInd = 0;
      int i = 0;
      int closeIndx = 0;
      for (int j = 0; j < n_closeGeomElements[0]; j++) {
        closeIndx = zInd * nR_closeGeom[0] * n_closeGeomElements[0] +
                    rInd * n_closeGeomElements[0] + j;
        // if(zInd*nR_closeGeom[0]*n_closeGeomElements[0] +
        // rInd*n_closeGeomElements[0] + j < 0)
        //{
        //        zInd=0;
        //        rInd=0;
        //        j=0;
        //    //std::cout << "index " <<
        //    zInd*nR_closeGeom[0]*n_closeGeomElements[0] +
        //    rInd*n_closeGeomElements[0] + j << std::endl;
        //}
        //    if(zInd*nR_closeGeom[0]*n_closeGeomElements[0] +
        //    rInd*n_closeGeomElements[0] + j > 1309440)
        //    {
        //        zInd=0;
        //        rInd=0;
        //        j=0;
        //        //std::cout << "index " <<
        //        zInd*nR_closeGeom[0]*n_closeGeomElements[0] +
        //        rInd*n_closeGeomElements[0] + j << std::endl;
        //    }
        i = closeGeom[closeIndx];
      
#else
      for (int i = 0; i < nLines; i++) {
#endif
        // std::cout << "vert geom " << i << "  " <<
        // fabs(boundaryVector[i].slope_dzdx) << " " << tol << std::endl;
        if (std::abs(boundaryVector[i].slope_dzdx) >= tol * 0.75f) 
        {
          signPoint = std::copysign(1.0, pdim1 - boundaryVector[i].x1);
          signPoint0 = std::copysign(1.0, pdim1previous - boundaryVector[i].x1);
          // std::cout << "signpoint1 " << signPoint << " " << signPoint0 <<
          // std::endl;
        } 
        else 
        {
          signPoint =
              std::copysign(1.0, particlesPointer->z[indx] -
                                     pdim1 * boundaryVector[i].slope_dzdx -
                                     boundaryVector[i].intercept_z);
          signPoint0 = std::copysign(1.0, particlesPointer->zprevious[indx] -
                                              pdim1previous *
                                                  boundaryVector[i].slope_dzdx -
                                              boundaryVector[i].intercept_z);
           //std::cout << "signpoint2 " << signPoint << " " << signPoint0 <<
           //std::endl;
        }

        if (signPoint != signPoint0) 
        {
          if (std::abs(particle_slope) >= tol * 0.75f) 
          {
            // std::cout << " isinf catch " << std::endl;
            particle_slope = tol;
          }
          if (std::abs(particle_slope) >= tol * 0.75f) 
          {
            signLine1 = std::copysign(1.0, boundaryVector[i].x1 - pdim1);
            signLine2 = std::copysign(1.0, boundaryVector[i].x2 - pdim1);
            // std::cout << "signlines3 " << signLine1 << " " << signLine2 <<
            // std::endl;
          }
          else 
          {
            signLine1 =
                std::copysign(1.0, boundaryVector[i].z1 -
                                       boundaryVector[i].x1 * particle_slope -
                                       particle_intercept);
            signLine2 =
                std::copysign(1.0, boundaryVector[i].z2 -
                                       boundaryVector[i].x2 * particle_slope -
                                       particle_intercept);
            //std::cout << "signline 1 and 2 " << signLine1 << " " << signLine2 << std::endl;
          }

          ////if (signPoint != signPoint0) 
          ////{
          ////  if (std::abs(particle_slope) >= tol * 0.75f) 
          ////  {
          ////    // std::cout << " isinf catch " << std::endl;
          ////    particle_slope = tol;
          ////  }
          ////  if (std::abs(particle_slope) >= tol * 0.75f) 
          ////  {
          ////    signLine1 = std::copysign(1.0, boundaryVector[i].x1 - pdim1);
          ////    signLine2 = std::copysign(1.0, boundaryVector[i].x2 - pdim1);
          ////    // std::cout << "signlines3 " << signLine1 << " " << signLine2 <<
          ////    // std::endl;
          ////  } 
          ////  else 
          ////  {
          ////    signLine1 =
          ////        std::copysign(1.0, boundaryVector[i].z1 -
          ////                               boundaryVector[i].x1 * particle_slope -
          ////                               particle_intercept);
          ////    signLine2 =
          ////        std::copysign(1.0, boundaryVector[i].z2 -
          ////                               boundaryVector[i].x2 * particle_slope -
          ////                               particle_intercept);
          ////  }
            // std::cout << "signLines " << signLine1 << " " << signLine2 <<
            // std::endl; std::cout << "bound vec points " <<
            // boundaryVector[i].z1 << " " << boundaryVector[i].x1 <<
            // " " << boundaryVector[i].z2 << " " << boundaryVector[i].x2 <<
            // std::endl;
            if (signLine1 != signLine2) 
            {
              intersectionIndices[nIntersections] = i;
              nIntersections++;

              // std::cout << "nintersections " << nIntersections << std::endl;
              // std::cout << fabs(particlesPointer->x[indx] -
              // particlesPointer->xprevious[indx]) << tol_small << std::endl;
              if (std::abs(pdim1 - pdim1previous) < tol_small) 
              {
                //  std::cout << "vertical line" << std::cout;
                intersectionx[nIntersections - 1] = pdim1previous;
                intersectiony[nIntersections - 1] =
                    intersectionx[nIntersections - 1] *
                        boundaryVector[i].slope_dzdx +
                    boundaryVector[i].intercept_z;
              } 
              else 
              {
                // std::cout << "not vertical line" << std::endl;
                // std::cout << 0.0*7.0 << " " << i << " " << nParam << " " <<
                // lines[i*nParam+4] << "  " <<tol << std::endl; std::cout <<
                // "boundaryVector slope " << boundaryVector[i].slope_dzdx << " "
                // << tol*0.75 <<std::endl;
                if (std::abs(boundaryVector[i].slope_dzdx) >= tol * 0.75f) 
                {
                  intersectionx[nIntersections - 1] = boundaryVector[i].x1;
                } 
                else 
                {
                  intersectionx[nIntersections - 1] =
                      (boundaryVector[i].intercept_z - particle_intercept) /
                      (particle_slope - boundaryVector[i].slope_dzdx);
                  //  std::cout << "in this else "<<
                  //  intersectionx[nIntersections -1] << std::endl;
                }
                intersectiony[nIntersections - 1] =
                    intersectionx[nIntersections - 1] * particle_slope +
                    particle_intercept;
                    //std::cout << "intersectionx and y"<<
                    //intersectionx[nIntersections -1] << " " << intersectiony[0] << std::endl;
              }
            }
          ////}
        }
      }
      //std::cout << " nIntersections " << nIntersections << std::endl;
        // if(particlesPointer->hitWall[indx] == 0.0)
        // {
        if (nIntersections == 0) 
        {
          particlesPointer->distTraveled[indx] =
              particlesPointer->distTraveled[indx] + dpath;
          particlesPointer->xprevious[indx] = particlesPointer->x[indx];
          particlesPointer->yprevious[indx] = particlesPointer->y[indx];
          particlesPointer->zprevious[indx] = particlesPointer->z[indx];
          // particlesPointer->test0[indx] = -50.0;

          // std::cout << "r " << particlesPointer->x[indx] << " " <<
          // particlesPointer->y[indx] << " " << particlesPointer->z[indx] <<
          // std::endl; std::cout << "r0 " << particlesPointer->xprevious[indx]
          // << " " << particlesPointer->yprevious[indx] << " " <<
          // particlesPointer->zprevious[indx] << std::endl;
        } 
        else if (nIntersections == 1) 
        {
          particlesPointer->hitWall[indx] = 1.0f;
          particlesPointer->wallIndex[indx] = intersectionIndices[0];
          particlesPointer->wallHit[indx] = intersectionIndices[0];
          // particlesPointer->test0[indx] = -100.0;
          if (particle_slope >= tol * 0.75f) 
          {
#if USECYLSYMM > 0
            float x0 = particlesPointer->xprevious[indx];
            float x1 = particlesPointer->x[indx];
            float y0 = particlesPointer->yprevious[indx];
            float y1 = particlesPointer->y[indx];
            float tt = findT(x0, x1, y0, y1, intersectionx[0]);
            xNew = x0 + (x1 - x0) * tt;
            yNew = y0 + (y1 - y0) * tt;
            rNew = std::sqrt(xNew * xNew + yNew * yNew);
            thetaNew = theta0 +
                       (intersectiony[0] - particlesPointer->zprevious[indx]) /
                           (particlesPointer->z[indx] -
                            particlesPointer->zprevious[indx]) *
                           (theta1 - theta0);
            particlesPointer->y[indx] = yNew;
            particlesPointer->yprevious[indx] = yNew;
#else
            // std::cout << "Particle index " << indx << " hit wall and is
            // calculating y point " << particlesPointer->y[indx] << std::endl;
            particlesPointer->y[indx] =
                particlesPointer->yprevious[indx] +
                (intersectiony[0] - particlesPointer->zprevious[indx]) /
                    (particlesPointer->z[indx] -
                     particlesPointer->zprevious[indx]) *
                    (particlesPointer->y[indx] -
                     particlesPointer->yprevious[indx]);
            // std::cout << "yprev,intersectiony,zprevious,z,y " <<
            // particlesPointer->yprevious[indx] << " " << intersectiony[0] << "
            // "<< particlesPointer->zprevious[indx] << " " <<
            // particlesPointer->z[indx] << " " << particlesPointer->y[indx] <<
            // std::endl;

#endif
          } 
          else 
          {
#if USECYLSYMM > 0
            float x0 = particlesPointer->xprevious[indx];
            float x1 = particlesPointer->x[indx];
            float y0 = particlesPointer->yprevious[indx];
            float y1 = particlesPointer->y[indx];
            float tt = findT(x0, x1, y0, y1, intersectionx[0]);
            xNew = x0 + (x1 - x0) * tt;
            yNew = y0 + (y1 - y0) * tt;
            rNew = std::sqrt(xNew * xNew + yNew * yNew);
            // particlesPointer->test0[indx] = -200.0;
            thetaNew = theta0 + (intersectionx[0] - pdim1previous) /
                                    (pdim1 - pdim1previous) * (theta1 - theta0);
            //std::cout << " tt xnew ynew " << tt << " " << xNew << " " << yNew << std::endl;
            particlesPointer->yprevious[indx] = yNew;
            particlesPointer->y[indx] = yNew;
            // float rrr  =
            // std::sqrt(particlesPointer->x[indx]*particlesPointer->x[indx] +
            // particlesPointer->y[indx]*particlesPointer->y[indx]);
            // if(particlesPointer->z[indx]< -4.1 & rrr > 5.5543)
            //{
            //  std::cout <<" positions of intersection 2" <<
            //  particlesPointer->x[indx] << " " << particlesPointer->y[indx]<<
            //  std::endl; std::cout <<" r " << rrr << " " <<
            //  boundaryVector[particlesPointer->wallHit[indx]].x1 << " " <<
            //  boundaryVector[particlesPointer->wallHit[indx]].x2<< std::endl;
            // std::cout << "x0 x1 y0 y1 rNew "  << " "<< x0 << " " << x1 << " "
            // << y0 << " " << y1 << " " << rNew << std::endl; std::cout << "xNew
            // yNew " << xNew << " " << yNew << std::endl; std::cout <<
            // "intersectionx " << intersectionx[0] << std::endl;
            //}
#else
            // std::cout << "Particle index " << indx << " hit wall and is
            // calculating y point " << particlesPointer->y[indx] << std::endl;
            particlesPointer->y[indx] =
                particlesPointer->yprevious[indx] +
                (intersectionx[0] - particlesPointer->xprevious[indx]) /
                    (particlesPointer->x[indx] -
                     particlesPointer->xprevious[indx]) *
                    (particlesPointer->y[indx] -
                     particlesPointer->yprevious[indx]);
            // std::cout << "yprev,intersectiony,zprevious,z,y " <<
            // particlesPointer->yprevious[indx] << " " << intersectiony[0] << "
            // "<< particlesPointer->zprevious[indx] << " " <<
            // particlesPointer->z[indx] << " " << particlesPointer->y[indx] <<
            // std::endl;
#endif
          }
#if USECYLSYMM > 0
          particlesPointer->xprevious[indx] = xNew;
          particlesPointer->x[indx] = particlesPointer->xprevious[indx];
#else
          particlesPointer->x[indx] = intersectionx[0];
          particlesPointer->xprevious[indx] = intersectionx[0];
#endif
          particlesPointer->zprevious[indx] = intersectiony[0];
          particlesPointer->z[indx] = intersectiony[0];
          // std::cout << "nInt = 1 position " << intersectionx[0] << " " <<
          // intersectiony[0]  << std::endl;
        } 
        else 
        {
          // std::cout << "nInts greater than 1 " << nIntersections <<
          // std::endl;
          for (int i = 0; i < nIntersections; i++) {
            distances[i] =
                (pdim1previous - intersectionx[i]) *
                    (pdim1previous - intersectionx[i]) +
                (particlesPointer->zprevious[indx] - intersectiony[i]) *
                    (particlesPointer->zprevious[indx] - intersectiony[i]);
            if (distances[i] < minDist) {
              minDist = distances[i];
              minDistInd = i;
            }
          }

          particlesPointer->wallIndex[indx] = intersectionIndices[minDistInd];
          particlesPointer->wallHit[indx] = intersectionIndices[minDistInd];
          particlesPointer->hitWall[indx] = 1.0f;
#if USECYLSYMM > 0
          thetaNew = theta0 + (intersectionx[minDistInd] - pdim1previous) /
                                  (pdim1 - pdim1previous) * (theta1 - theta0);
          particlesPointer->yprevious[indx] =
              intersectionx[minDistInd] * std::sin(thetaNew);
          particlesPointer->y[indx] = particlesPointer->yprevious[indx];
          // particlesPointer->y[indx] =
          // intersectionx[minDistInd]*cosf(thetaNew);
          particlesPointer->x[indx] =
              intersectionx[minDistInd] * std::cos(thetaNew);
#else
          // std::cout << "Particle index " << indx << " hit wall and is
          // calculating y point " << particlesPointer->yprevious[indx] << " " <<
          // particlesPointer->y[indx] << std::endl;
          // particlesPointer->y[indx] = particlesPointer->yprevious[indx] +
          // (intersectionx[minDistInd] - pdim1previous)/(pdim1 -
          // pdim1previous)*(particlesPointer->y[indx] -
          // particlesPointer->yprevious[indx]); std::cout <<
          // "intersectionx,pdp,pd,y " << intersectionx[0] << " "<< pdim1previous
          // << " " << pdim1 << " " << particlesPointer->y[indx] << std::endl;
          particlesPointer->y[indx] =
              particlesPointer->yprevious[indx] +
              (intersectiony[0] - particlesPointer->zprevious[indx]) /
                  (particlesPointer->z[indx] -
                   particlesPointer->zprevious[indx]) *
                  (particlesPointer->y[indx] -
                   particlesPointer->yprevious[indx]);
          // std::cout << "yprev,intersectiony,zprevious,z,y " <<
          // particlesPointer->yprevious[indx] << " " << intersectiony[0] << "
          // "<< particlesPointer->zprevious[indx] << " " <<
          // particlesPointer->z[indx] << " " << particlesPointer->y[indx] <<
          // std::endl;
          particlesPointer->x[indx] = intersectionx[minDistInd];
#endif
          particlesPointer->z[indx] = intersectiony[minDistInd];
        }
      ////}
        // else
        //{
        //    if (particlesPointer->y[indx] < boundaryVector[nLines].y1)
        //    {
        //        particlesPointer->hitWall[indx] = 1.0f;
        //    }
        //    else if (particlesPointer->y[indx] > boundaryVector[nLines].y2)
        //    {
        //        particlesPointer->hitWall[indx] = 1.0f;
        //    }
        //}
#endif
      if (particlesPointer->hitWall[indx] == 1.0) {

#if (FLUX_EA > 0 && USESURFACEMODEL == 0)
        float E0 = 0.0;
        float thetaImpact = 0.0;
        float particleTrackVector[3] = {0.0f};
        float surfaceNormalVector[3] = {0.0f};
        float norm_part = 0.0;
        float partDotNormal = 0.0;
        particleTrackVector[0] = particlesPointer->vx[indx];
        particleTrackVector[1] = particlesPointer->vy[indx];
        particleTrackVector[2] = particlesPointer->vz[indx];
        norm_part = std::sqrt(particleTrackVector[0] * particleTrackVector[0] +
                         particleTrackVector[1] * particleTrackVector[1] +
                         particleTrackVector[2] * particleTrackVector[2]);
        E0 = 0.5 * particlesPointer->amu[indx] * 1.6737236e-27 *
             (norm_part * norm_part) / 1.60217662e-19;
        int wallHitP = particlesPointer->wallHit[indx];
        boundaryVector[particlesPointer->wallHit[indx]].getSurfaceNormal(
            surfaceNormalVector, particlesPointer->y[indx],
            particlesPointer->x[indx]);
        particleTrackVector[0] = particleTrackVector[0] / norm_part;
        particleTrackVector[1] = particleTrackVector[1] / norm_part;
        particleTrackVector[2] = particleTrackVector[2] / norm_part;

        partDotNormal =
            vectorDotProduct(particleTrackVector, surfaceNormalVector);
        thetaImpact = std::acos(partDotNormal);
        if (thetaImpact > 3.14159265359 * 0.5) {
          thetaImpact = std::abs(thetaImpact - (3.14159265359));
        }
        thetaImpact = thetaImpact * 180.0 / 3.14159265359;
        EdistInd = std::floor((E0 - E0dist) / dEdist);
        AdistInd = std::floor((thetaImpact - A0dist) / dAdist);
        int surfaceHit =
            boundaryVector[particlesPointer->wallHit[indx]].surfaceNumber;
        int surface = boundaryVector[particlesPointer->wallHit[indx]].surface;
        float weight = particlesPointer->weight[indx];
        // particlesPointer->test[indx] = norm_part;
        // particlesPointer->test0[indx] = partDotNormal;
        // particlesPointer->test1[indx] = particleTrackVector[0];
        // particlesPointer->test2[indx] = particleTrackVector[1];
        // particlesPointer->test3[indx] = particleTrackVector[2];
        // particlesPointer->test4[indx] = particles;
        // std::cout << "impact energy and angle " << E0 << " " << thetaImpact
        // << std::endl; std::cout << "surface EAinds " <<surface<< " " <<
        // EdistInd << " " << AdistInd << std::endl;
        if (surface) {
          if ((EdistInd >= 0) && (EdistInd < nEdist) && (AdistInd >= 0) &&
              (AdistInd < nAdist)) {
#if USE_CUDA > 0
            atomicAdd(
                &surfaces->energyDistribution[surfaceHit * nEdist * nAdist +
                                              EdistInd * nAdist + AdistInd],
                particlesPointer->weight[indx]);
            atomicAdd(&surfaces->grossDeposition[surfaceHit],
                      particlesPointer->weight[indx]);
            atomicAdd(&surfaces->sumWeightStrike[surfaceHit],
                      particlesPointer->weight[indx]);
            atomicAdd(&surfaces->sumParticlesStrike[surfaceHit], 1);
#else

            surfaces->energyDistribution[surfaceHit * nEdist * nAdist +
                                         EdistInd * nAdist + AdistInd] =
                surfaces->energyDistribution[surfaceHit * nEdist * nAdist +
                                             EdistInd * nAdist + AdistInd] +
                weight;
            surfaces->sumWeightStrike[surfaceHit] =
                surfaces->sumWeightStrike[surfaceHit] + weight;
            surfaces->sumParticlesStrike[surfaceHit] =
                surfaces->sumParticlesStrike[surfaceHit] + 1;
            surfaces->grossDeposition[surfaceHit] =
                surfaces->grossDeposition[surfaceHit] + weight;
#endif
          }
        }
#elif (FLUX_EA == 0 && USESURFACEMODEL == 0)
          particlesPointer->weight[indx] = 0.0;
#endif
        // particlesPointer->transitTime[indx] = tt*dt;
      }
    }

    // std::cout << "2geometry check particle x" << particlesPointer->x[indx] <<
    // particlesPointer->x[indx]previous <<std::endl; std::cout << "2geometry
    // check particle y" << particlesPointer->y[indx] <<
    // particlesPointer->y[indx]previous <<std::endl; std::cout << "2geometry
    // check particle z" << particlesPointer->z[indx] <<
    // particlesPointer->z[indx]previous <<std::endl;
  }
};

#endif
