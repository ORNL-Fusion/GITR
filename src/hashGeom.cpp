#include "hashGeom.h"

hashGeom::hashGeom( int _nLines,int _nHashes,
                Boundary* _boundary,
                gitr_precision* _x,
                gitr_precision* _y, 
                gitr_precision* _z, 
                int* _n_closeGeomElements,//gitr_precision *_minDist,
                int *_closeGeom,
                int* _nR, int* _nY, int* _nZ, int use_3d_geom_ )
               :  nLines(_nLines),nHashes(_nHashes),boundary(_boundary), x(_x), y(_y), z(_z), 
               n_closeGeomElements(_n_closeGeomElements), 
               closeGeom(_closeGeom), nR(_nR), nY(_nY), nZ(_nZ), use_3d_geom( use_3d_geom_ )
{ }

CUDA_CALLABLE_MEMBER_DEVICE 
void hashGeom::operator()(std::size_t indx) {
  //printf(" indx %i \n", indx);
  int nHash=0;
  //std::cout << "nHashes "<<nHashes << std::endl;
  int hashSum=0;
  int nRhashSum=0;
  int nYhashSum=0;
  int nZhashSum=0;
  int nHashPoints=0;
  /* Captain! Iterate over the number of hash files */
  for(int i=0;i<nHashes;i++)
  {  
    nRhashSum = nRhashSum + nR[i];
    nYhashSum = nYhashSum + nY[i];
    nZhashSum = nZhashSum + nZ[i];
    nHashPoints = nHashPoints+nR[i]*nY[i]*nZ[i];
    hashSum = hashSum + nR[i]*nY[i]*nZ[i]*n_closeGeomElements[i];
    if(indx >= nHashPoints)
    {nHash = nHash +1;}
  }
  //std::cout << "nHash " << nHash << std::endl;
  hashSum=0;
  nRhashSum=0;
  nYhashSum=0;
  nZhashSum=0;
  nHashPoints=0;
  for(int i=0;i<nHash;i++)
  {  
    nRhashSum = nRhashSum + nR[i];
    nYhashSum = nYhashSum + nY[i];
    nZhashSum = nZhashSum + nZ[i];
    nHashPoints = nHashPoints+nR[i]*nY[i]*nZ[i];
    hashSum = hashSum + nR[i]*nY[i]*nZ[i]*n_closeGeomElements[i];
  }
  //std::cout << "index " << indx << std::endl;
  //std::cout << "hashSum " << hashSum << std::endl;
  //std::cout << "nR[nHash] " << nR[nHash] << std::endl;
  //std::cout << "nY[nHash] " << nY[nHash] << std::endl;
  //std::cout << "nRhashSum " << nRhashSum << std::endl;
  //std::cout << "nYhashSum " << nYhashSum << std::endl;
  //std::cout << "nZhashSum " << nZhashSum << std::endl;
  gitr_precision x0;

  gitr_precision y0;

  gitr_precision z0;

  int xyzIndx;

  int buffIndx;

    if( use_3d_geom > 0 )
    {
  gitr_precision kk = (indx-nHashPoints)/(nR[nHash]*nY[nHash]);
  //std::cout << "kk " << kk << std::endl;

  int k = std::floor(kk);
  //std::cout << "k " << k << std::endl;
  int jjj = (indx-nHashPoints) - k*nR[nHash]*nY[nHash];
  //std::cout << "jjj " << jjj << std::endl;
  gitr_precision jj = 1.0*jjj/nR[nHash];
  //std::cout << "jj " << jj << std::endl;
  int j = std::floor(jj);
  //std::cout << "j " << j << std::endl;
  int i = (indx-nHashPoints)- j*nR[nHash] - k*(nR[nHash]*nY[nHash]);
  //std::cout << "i " << i << std::endl;

  //gitr_precision jj = indx/nR;
  //int j = floor(jj);
  //int i = indx - j*nR;
  //int xyzIndx = k*nR*nY + indx;
  //if( i > nR || i < 0){ std::cout << "i out of range " << i << std::endl; exit(0);}
  //if( j > nY || j < 0){ std::cout << "j out of range " << j << std::endl; exit(0);}
  //if( k > nZ || k < 0){ std::cout << "k out of range " << k  << "indx " << indx<< std::endl; exit(0);}
  //std::cout << "ijk " << i << " " << j << " "<< k << std::endl;
  xyzIndx = indx;
  buffIndx = hashSum+(k*(nR[nHash]*nY[nHash])+j*nR[nHash]+i)*n_closeGeomElements[nHash] ;
  x0 = x[nRhashSum+i];
  y0 = y[nYhashSum+j];
  z0 = z[nZhashSum+k];
  //std::cout << "point "  << nHash << " " <<   x0 << " " <<  y0 << " "
  //     <<  z0 << std::endl;
    }
    else
    {
  nHash=0;
  hashSum=0;
  nRhashSum=0;
  nYhashSum=0;
  nZhashSum=0;
  nHashPoints=0;
  gitr_precision kk = indx/(nR[0]);
  int k = std::floor(kk);
  int i = indx - k*(nR[0]);
  x0 = x[i];
  y0 = 0.0;
  z0 = z[k];
  xyzIndx = indx;
  buffIndx=(k*(nR[0])+ i)*n_closeGeomElements[0];
  //std::cout << "point "  <<nHash<< " " <<   x0 << " " <<  z0 << " "
  //     <<  buffIndx << std::endl;

    }
  //gitr_precision minDist[n_closeGeomElements] = {0.0};
  //for(int i1=0;i1<n_closeGeomElements; i1++)
  //{
  //  minDist[i1] = 1.0e6;
  //  //closeGeom[indx*n_closeGeomElements + i1] = indx;
  //}
  gitr_precision A[3] = {0.0,0.0,0.0};
  gitr_precision B[3] = {0.0,0.0,0.0};
  gitr_precision C[3] = {0.0,0.0,0.0};
  gitr_precision AB[3] = {0.0,0.0,0.0};
  gitr_precision AC[3] = {0.0,0.0,0.0};
  gitr_precision BC[3] = {0.0,0.0,0.0};
  gitr_precision CA[3] = {0.0,0.0,0.0};
  gitr_precision p[3] = {0.0,0.0,0.0};
  gitr_precision Ap[3] = {0.0,0.0,0.0};
  gitr_precision Bp[3] = {0.0,0.0,0.0};
  gitr_precision Cp[3] = {0.0,0.0,0.0};
  gitr_precision normalVector[3] = {0.0,0.0,0.0};
  gitr_precision crossABAp[3] = {0.0,0.0,0.0};
  gitr_precision crossBCBp[3] = {0.0,0.0,0.0};
  gitr_precision crossCACp[3] = {0.0,0.0,0.0};
  gitr_precision signDot0 = 0.0;
  gitr_precision signDot1 = 0.0;
  gitr_precision signDot2 = 0.0;
  gitr_precision totalSigns = 0.0;
#if USE_CUDA
  gitr_precision *minDist  = new gitr_precision[n_closeGeomElements[nHash]];
  for(int i1=0;i1<n_closeGeomElements[nHash];i1++){ minDist[i1] = 1.0e6;}
  //gitr_precision minDist[n_closeGeomElements[nHash]] = {1.0e6};
#else
  sim::Array<gitr_precision> minDist(n_closeGeomElements[nHash],1e6);      
#endif
  for(int l=0; l<nLines; l++)
  {
    //  if(indx ==1)
    //  {std::cout << "l minDist1" << l << " " << minDist[0]<< std::endl;}
    //std::cout << " line l xyz " << l << " " <<  boundary[l].x1 << " " <<  boundary[l].y1 << " "
    //     <<  boundary[l].z1 << std::endl;
    //std::cout << " xyz 2 " <<  boundary[l].x2 << " " <<  boundary[l].y2 << " "
    //     <<  boundary[l].z2 << std::endl;
    //std::cout << "xyz 3 "  <<  boundary[l].x3 << " " <<  boundary[l].y3 << " "
    //     <<  boundary[l].z3 << std::endl;
    gitr_precision a = boundary[l].a;
    gitr_precision b = boundary[l].b;
    gitr_precision c = boundary[l].c;
    gitr_precision d = boundary[l].d;
    //std::cout << "abcd "  << a << " " <<b << " "
    //    <<  c << " " <<d << std::endl;
    gitr_precision perpDist;
    if( use_3d_geom > 0 )
    {
    gitr_precision plane_norm = boundary[l].plane_norm;
    gitr_precision t = -(a*x0 + b*y0 + c*z0 + d)/(a*a + b*b + c*c);
    p[0] = a*t + x0;
    p[1] = b*t + y0;
    p[2] = c*t + z0;
    perpDist = std::sqrt((x0-p[0])*(x0-p[0]) + (y0-p[1])*(y0-p[1]) + (z0-p[2])*(z0-p[2]));
    }

    vectorAssign(boundary[l].x1, boundary[l].y1, 
        boundary[l].z1, A);    
    vectorAssign(boundary[l].x2, boundary[l].y2, 
        boundary[l].z2, B);    
    if( use_3d_geom > 0 )
    {
    vectorAssign(boundary[l].x3, boundary[l].y3, 
        boundary[l].z3, C); 
    }
    vectorSubtract(B,A,AB);
    if( use_3d_geom > 0 )
    {
    vectorSubtract(C,A,AC);
    vectorSubtract(C,B,BC);
    vectorSubtract(A,C,CA);

    vectorSubtract(p,A,Ap);
    vectorSubtract(p,B,Bp);
    vectorSubtract(p,C,Cp);

    vectorCrossProduct(AB,AC,normalVector);
    vectorCrossProduct(AB,Ap,crossABAp);
    vectorCrossProduct(BC,Bp,crossBCBp);
    vectorCrossProduct(CA,Cp,crossCACp);

    signDot0 = std::copysign(1.0,vectorDotProduct(crossABAp, normalVector));
    signDot1 = std::copysign(1.0,vectorDotProduct(crossBCBp, normalVector));
    signDot2 = std::copysign(1.0,vectorDotProduct(crossCACp, normalVector));
    totalSigns = std::abs(signDot0 + signDot1 + signDot2);

    if (totalSigns == 3.0)
    {
    }
    else perpDist = 1.0e6;
    }
    //std::cout << "perpDist " << perpDist << std::endl;
    //Edge checking
    p[0] = x0;
    p[1] = y0;
    p[2] = z0;
    gitr_precision pA[3] = {0.0};
    gitr_precision cEdge1[3] = {0.0};
    gitr_precision dEdge1[3] = {0.0};
    vectorSubtract(A,p,pA);
    gitr_precision cEdge1mag = vectorDotProduct(pA,AB)/vectorDotProduct(AB,AB);
    gitr_precision distE1 = 1.0e6;
    if(cEdge1mag < 0.0 && cEdge1mag > -1.0)
    {
      vectorScalarMult(cEdge1mag,AB,cEdge1);
      vectorSubtract(pA,cEdge1,dEdge1);
      distE1 = std::sqrt(vectorDotProduct(dEdge1,dEdge1));
    }
    //std::cout << "edge1 comp " << pA[0] << " " << pA[1] << " " << pA[2] <<
    //   " " << cEdge1mag << " " << cEdge1[0] << " " << cEdge1[1] << " " << cEdge1[2] << " "
    //   << dEdge1[0] << " " <<dEdge1[1] << " " << dEdge1[2] << std::endl;
    gitr_precision minEdge;
    if( use_3d_geom > 0 )
    {
    gitr_precision pB[3] = {0.0};
    gitr_precision cEdge2[3] = {0.0};
    gitr_precision dEdge2[3] = {0.0};
    vectorSubtract(B,p,pB);
    gitr_precision cEdge2mag = vectorDotProduct(pB,BC)/vectorDotProduct(BC,BC);
    gitr_precision distE2 = 1.0e6;
    if(cEdge2mag < 0.0 && cEdge2mag > -1.0)
    {
      vectorScalarMult(cEdge2mag,BC,cEdge2);
      vectorSubtract(pB,cEdge2,dEdge2);
      distE2 = std::sqrt(vectorDotProduct(dEdge2,dEdge2));
    }
    gitr_precision pC[3] = {0.0};
    gitr_precision cEdge3[3] = {0.0};
    gitr_precision dEdge3[3] = {0.0};
    vectorSubtract(C,p,pC);
    gitr_precision cEdge3mag = vectorDotProduct(pC,CA)/vectorDotProduct(CA,CA);
    gitr_precision distE3 = 1.0e6;
    if(cEdge3mag < 0.0 && cEdge3mag > -1.0)
    {
      vectorScalarMult(cEdge3mag,CA,cEdge3);
      vectorSubtract(pC,cEdge3,dEdge3);
      distE3 = std::sqrt(vectorDotProduct(dEdge3,dEdge3));
    }
    minEdge = std::min(distE1,distE2);
    minEdge = std::min(distE3,minEdge);
    }
    else
    {
    //
    minEdge = distE1;
    }
    //std::cout << "edgeDistances " << distE1 << " " << distE2 << " " << distE3 << std::endl;
    gitr_precision d1 =std::sqrt((x0 - boundary[l].x1)*(x0 - boundary[l].x1)
        +  (y0 - boundary[l].y1)*(y0 - boundary[l].y1)
        +  (z0 - boundary[l].z1)*(z0 - boundary[l].z1));
    gitr_precision d2 =std::sqrt((x0 - boundary[l].x2)*(x0 - boundary[l].x2)
        +  (y0 - boundary[l].y2)*(y0 - boundary[l].y2)
        +  (z0 - boundary[l].z2)*(z0 - boundary[l].z2));
    gitr_precision d3;
    if( use_3d_geom > 0 )
    {
      d3 = std::sqrt((x0 - boundary[l].x3)*(x0 - boundary[l].x3)
        +  (y0 - boundary[l].y3)*(y0 - boundary[l].y3)
        +  (z0 - boundary[l].z3)*(z0 - boundary[l].z3));
    }
    //std::cout << " point Distances " << d3 << " " << d2 << " " << d1 << std::endl;
    gitr_precision minOf3 = std::min(d1,d2);
    minOf3 = std::min(minOf3,minEdge);
    //std::cout << "min of two " << minOf3 << std::endl;
    if( use_3d_geom > 0 )
    {
    minOf3 = std::min(minOf3,perpDist);
    minOf3 = std::min(minOf3,d3);
    }
    //std::cout << "mindist "  << minOf3 << " " <<  std::endl;
    //  if(indx ==1)
    //  {std::cout << "minof3" << perpDist <<  " " << minEdge << " " << minOf3<< std::endl;}
    int minIndClose = n_closeGeomElements[nHash];
    for(int m=0; m< n_closeGeomElements[nHash]; m++)
    {
      //  if(indx ==1)
      //  {std::cout << "minDist" << minDist[m] << std::endl;}
      //if(minDist[xyzIndx*n_closeGeomElements + m] > minOf3)
      if(minDist[m] > minOf3)
      {
        minIndClose = minIndClose-1;
      }
    }

    if((minIndClose < n_closeGeomElements[nHash]) && (minIndClose > -1))
    {
      //std::cout << "min INd close " << l << std::endl;
      //%shift numbers down
      for(int n=n_closeGeomElements[nHash]-1; n>minIndClose; n--)
      {
        //minDist[xyzIndx*n_closeGeomElements + n] = 
        //minDist[xyzIndx*n_closeGeomElements + n-1];  
        minDist[n] = 
          minDist[n-1];  
        closeGeom[buffIndx+ n] =    
          closeGeom[buffIndx + n-1];
      }
      //minDist[xyzIndx*n_closeGeomElements + minIndClose] = minOf3;
      minDist[minIndClose] = minOf3;
      closeGeom[buffIndx + minIndClose] = l;
      //if(indx ==1)
      //{std::cout << "l minof3" << l << " " << minOf3<< std::endl;}
      //    if((indx*n_closeGeomElements + minIndClose) ==10)
      //    {
      //        if(indx > 1) std::cout << "this is the mess up " << indx << " " << n_closeGeomElements << " " << minIndClose << std::endl;
      //  }
    }
    /*     
           if(l == nLines - 1)
           {
           for(int o=0;o<n_closeGeomElements;o++)
           {
           std::cout << closeGeom[xyzIndx*n_closeGeomElements + o] << " ";
           }
           std::cout << std::endl;
           }
     */
  }
#if USE_CUDA
  delete[] minDist;
#endif
}
