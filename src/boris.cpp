#include "boris.h"
#include "constants.h"

/* Are these preprocessor defines necessary? */
#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

CUDA_CALLABLE_MEMBER
void vectorAdd(gitr_precision A[], gitr_precision B[],gitr_precision C[])
{
    C[0] = A[0] + B[0];
    C[1] = A[1] + B[1];
    C[2] = A[2] + B[2];
}

CUDA_CALLABLE_MEMBER
void vectorSubtract(gitr_precision A[], gitr_precision B[],gitr_precision C[])
{
    C[0] = A[0] - B[0];
    C[1] = A[1] - B[1];
    C[2] = A[2] - B[2];
}

CUDA_CALLABLE_MEMBER
void vectorScalarMult(gitr_precision a, gitr_precision B[],gitr_precision C[])
{
    C[0] = a*B[0];
    C[1] = a*B[1];
    C[2] = a*B[2];
}

CUDA_CALLABLE_MEMBER
void vectorAssign(gitr_precision a, gitr_precision b,gitr_precision c, gitr_precision D[])
{
    D[0] = a;
    D[1] = b;
    D[2] = c;
}

CUDA_CALLABLE_MEMBER
gitr_precision vectorNorm(gitr_precision A[])
{
    gitr_precision norm = 0.0;
    norm = std::sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);

    return norm;
}

CUDA_CALLABLE_MEMBER
void vectorNormalize(gitr_precision A[],gitr_precision B[])
{
    gitr_precision norm = 0.0;
    norm = std::sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
    B[0] = A[0]/norm;
    B[1] = A[1]/norm;
    B[2] = A[2]/norm;

}

CUDA_CALLABLE_MEMBER
gitr_precision vectorDotProduct(gitr_precision A[], gitr_precision B[])
{
    gitr_precision c = A[0]*B[0] +  A[1]*B[1] + A[2]*B[2];
    return c;
}

CUDA_CALLABLE_MEMBER
void vectorCrossProduct(gitr_precision A[], gitr_precision B[], gitr_precision C[])
{
    gitr_precision tmp[3] = {0.0,0.0,0.0};
    tmp[0] = A[1]*B[2] - A[2]*B[1];
    tmp[1] = A[2]*B[0] - A[0]*B[2];
    tmp[2] = A[0]*B[1] - A[1]*B[0];

    C[0] = tmp[0];
    C[1] = tmp[1];
    C[2] = tmp[2];
}

CUDA_CALLABLE_MEMBER
gitr_precision getE ( gitr_precision x0, gitr_precision y, gitr_precision z, gitr_precision E[], Boundary *boundaryVector, int nLines,
       int nR_closeGeom, int nY_closeGeom,int nZ_closeGeom, int n_closeGeomElements, 
       gitr_precision *closeGeomGridr,gitr_precision *closeGeomGridy, gitr_precision *closeGeomGridz, int *closeGeom, 
         int&  closestBoundaryIndex, int biased_surface, int use_3d_geom, 
         int geom_hash_sheath, int cylsymm  ) 
    {

    gitr_precision pot = 0.0;
    int minIndex = 0;
    gitr_precision minDistance = 1.0e12;
    gitr_precision Emag = 0.0;
    gitr_precision angle = 0.0;
    gitr_precision fd = 0.0;
    gitr_precision directionUnitVector[ 3 ] = { 0.0, 0.0, 0.0 };
    gitr_precision Er = 0.0;
    gitr_precision Et = 0.0;

    if( use_3d_geom > 0 )
    {
    gitr_precision p0[3] = {x0,y,z};
      gitr_precision a = 0.0;
      gitr_precision b = 0.0;
      gitr_precision c = 0.0;
      gitr_precision d = 0.0;
      gitr_precision plane_norm = 0.0;
      gitr_precision pointToPlaneDistance0 = 0.0;
      gitr_precision pointToPlaneDistance1 = 0.0;
      gitr_precision signPoint0 = 0.0;
      gitr_precision signPoint1 = 0.0;
      gitr_precision t = 0.0;
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
      gitr_precision p0A[3] = {0.0,0.0,0.0};
      gitr_precision p0B[3] = {0.0,0.0,0.0};
      gitr_precision p0C[3] = {0.0,0.0,0.0};
      gitr_precision p0AB[3] = {0.0,0.0,0.0};
      gitr_precision p0BC[3] = {0.0,0.0,0.0};
      gitr_precision p0CA[3] = {0.0,0.0,0.0};
      gitr_precision p0Anorm = 0.0;
      gitr_precision p0Bnorm = 0.0;
      gitr_precision p0Cnorm = 0.0;
      gitr_precision normalVector[3] = {0.0,0.0,0.0};
      gitr_precision crossABAp[3] = {0.0,0.0,0.0};
      gitr_precision crossBCBp[3] = {0.0,0.0,0.0};
      gitr_precision crossCACp[3] = {0.0,0.0,0.0};
      gitr_precision dot0 = 0.0;
      gitr_precision dot1 = 0.0;
      gitr_precision dot2 = 0.0;

      gitr_precision normAB = 0.0;
      gitr_precision normBC = 0.0;
      gitr_precision normCA = 0.0;
      gitr_precision ABhat[3] = {0.0,0.0,0.0};
      gitr_precision BChat[3] = {0.0,0.0,0.0};
      gitr_precision CAhat[3] = {0.0,0.0,0.0};
      gitr_precision tAB = 0.0;
      gitr_precision tBC = 0.0;
      gitr_precision tCA = 0.0;
      gitr_precision projP0AB[3] = {0.0,0.0,0.0};
      gitr_precision projP0BC[3] = {0.0,0.0,0.0};
      gitr_precision projP0CA[3] = {0.0,0.0,0.0};
      gitr_precision p0ABdist = 0.0;
      gitr_precision p0BCdist = 0.0;
      gitr_precision p0CAdist = 0.0;
      gitr_precision perpDist = 0.0;
      gitr_precision signDot0 = 0.0;
      gitr_precision signDot1 = 0.0;
      gitr_precision signDot2 = 0.0;
      gitr_precision totalSigns = 0.0;
      minDistance = 1.0e12;
      int nBoundariesCrossed = 0;
      int boundariesCrossed[6] = {0,0,0,0,0,0};
      gitr_precision distances[7] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
      gitr_precision normals[21] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                           0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                           0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  int top_limit = -1;

  gitr_precision dr;
  gitr_precision dy;
  gitr_precision dz;

  int rInd;
  int yInd;
  int zInd;

  if( geom_hash_sheath > 0 )
  {
    dr = closeGeomGridr[1] - closeGeomGridr[0];
    dy = closeGeomGridy[1] - closeGeomGridy[0];
    dz = closeGeomGridz[1] - closeGeomGridz[0];
    rInd = std::floor((x0 - closeGeomGridr[0])/dr + 0.5);
    yInd = std::floor((y - closeGeomGridy[0])/dy + 0.5);
    zInd = std::floor((z - closeGeomGridz[0])/dz + 0.5);

    if(rInd < 0 || rInd >= nR_closeGeom) rInd =0;
    
    if(yInd < 0 || yInd >= nY_closeGeom) yInd =0;

    if(zInd < 0 || zInd >= nZ_closeGeom) zInd =0;

    top_limit = n_closeGeomElements;
  }

  else top_limit = nLines;

  for (int k=0; k < top_limit; k++) //n_closeGeomElements
  {
    int i = -1;

    if( geom_hash_sheath > 0 )
    {
       i = closeGeom[zInd*nY_closeGeom*nR_closeGeom*n_closeGeomElements 
                   + yInd*nR_closeGeom*n_closeGeomElements
                   + rInd*n_closeGeomElements + k];
    }

    else
    {
      i = k;
    }

    a = boundaryVector[i].a;
    b = boundaryVector[i].b;
    c = boundaryVector[i].c;
    d = boundaryVector[i].d;
    plane_norm = boundaryVector[i].plane_norm;
    pointToPlaneDistance0 = (a * p0[0] + b * p0[1] + c * p0[2] + d) / plane_norm;
    vectorAssign(a / plane_norm, b / plane_norm, c / plane_norm, normalVector);
    vectorAssign(p0[0] - pointToPlaneDistance0 * normalVector[0],
                 p0[1] - pointToPlaneDistance0 * normalVector[1],
                 p0[2] - pointToPlaneDistance0 * normalVector[2], p);

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
    signDot0 = std::copysign(1.0,vectorDotProduct(crossABAp, normalVector));
    signDot1 = std::copysign(1.0,vectorDotProduct(crossBCBp, normalVector));
    signDot2 = std::copysign(1.0,vectorDotProduct(crossCACp, normalVector));
         totalSigns = std::abs(signDot0 + signDot1 + signDot2);
         vectorSubtract(A,p0,p0A);
         vectorSubtract(B,p0,p0B);
         vectorSubtract(C,p0,p0C);
         
         p0Anorm = vectorNorm(p0A);   
         p0Bnorm = vectorNorm(p0B);   
         p0Cnorm = vectorNorm(p0C);
         distances[1] = p0Anorm;   
         distances[2] = p0Bnorm;   
         distances[3] = p0Cnorm;   
             normals[3] = p0A[0]/p0Anorm;
             normals[4] = p0A[1]/p0Anorm;
             normals[5] = p0A[2]/p0Anorm;
             normals[6] = p0B[0]/p0Bnorm;
             normals[7] = p0B[1]/p0Bnorm;
             normals[8] = p0B[2]/p0Bnorm;
             normals[9] = p0C[0]/p0Cnorm;
             normals[10] = p0C[1]/p0Cnorm;
             normals[11] = p0C[2]/p0Cnorm;
         normAB = vectorNorm(AB);
         normBC = vectorNorm(BC);
         normCA = vectorNorm(CA);
         vectorAssign(AB[0]/normAB,AB[1]/normAB,AB[2]/normAB,ABhat);
         vectorAssign(BC[0]/normBC,BC[1]/normBC,BC[2]/normBC,BChat);
         vectorAssign(CA[0]/normCA,CA[1]/normCA,CA[2]/normCA,CAhat);
         
         tAB = vectorDotProduct(p0A,ABhat);
         tBC = vectorDotProduct(p0B,BChat);
         tCA = vectorDotProduct(p0C,CAhat);
         tAB = -1.0*tAB;
         tBC = -1.0*tBC;
         tCA = -1.0*tCA;
         if((tAB > 0.0) && (tAB < normAB))
         {
             vectorScalarMult(tAB,ABhat,projP0AB);
             vectorAdd(A,projP0AB,projP0AB);
             vectorSubtract(projP0AB,p0,p0AB);
             p0ABdist = vectorNorm(p0AB);
             distances[4] = p0ABdist;   
             normals[12] = p0AB[0]/p0ABdist;
             normals[13] = p0AB[1]/p0ABdist;
             normals[14] = p0AB[2]/p0ABdist;

         }
         else
         {
             p0ABdist = 1.0e12;
             distances[4] = p0ABdist;   
         } 
         
         
         if((tBC > 0.0) && (tBC < normBC))
         {
             vectorScalarMult(tBC,ABhat,projP0BC);
             vectorAdd(B,projP0BC,projP0BC);
             vectorSubtract(projP0BC,p0,p0BC);
             p0BCdist = vectorNorm(p0BC);
             distances[5] = p0BCdist;   
             normals[15] = p0BC[0]/p0BCdist;
             normals[16] = p0BC[1]/p0BCdist;
             normals[17] = p0BC[2]/p0BCdist;

         }
         else
         {
             p0BCdist = 1.0e12;
             distances[5] = p0BCdist;   

         } 
         
         if((tCA > 0.0) && (tCA < normCA))
         {
             vectorScalarMult(tCA,CAhat,projP0CA);
             vectorAdd(C,projP0CA,projP0CA);
             //std::cout << "projP0CA " << projP0CA[0] << " " << projP0CA[1] << " " << projP0CA[2] << std::endl; 
             vectorSubtract(projP0CA,p0,p0CA);
             p0CAdist = vectorNorm(p0CA);
             distances[6] = p0CAdist;   
             normals[18] = p0CA[0]/p0CAdist;
             normals[19] = p0CA[1]/p0CAdist;
             normals[20] = p0CA[2]/p0CAdist;
             //std::cout << "p0CA " << p0CA[0] << " " << p0CA[1] << " " << p0CA[2] << std::endl; 
         }
         else
         {
             p0CAdist = 1.0e12;
             distances[6] = p0CAdist;   
         } 

         if (totalSigns == 3.0)
         {
             //if (fabs(pointToPlaneDistance0) < minDistance)
             //{
                perpDist = std::abs(pointToPlaneDistance0); 
                vectorSubtract(p,p0 ,normalVector);
                vectorNormalize(normalVector,normalVector);
             distances[0] = perpDist;   
             normals[0] = boundaryVector[i].unit_vec0; //normalVector[0];
             normals[1] = boundaryVector[i].unit_vec1; //normalVector[1];
             normals[2] = boundaryVector[i].unit_vec2; //normalVector[2];
             //}
         }
         else
         {
             perpDist = 1.0e12;
             distances[0] = perpDist;   
         }
         int index = 0;
         for(int j = 0; j < 7; j++)
         {
            if(distances[j] < distances[index])
            index = j;              
         }

         if (distances[index] < minDistance)
         {
                 minDistance = distances[index];
                 vectorAssign(normals[index*3], normals[index*3+1],normals[index*3+2], directionUnitVector);
                 //std::cout << "min dist " << minDistance << std::endl;
                 //std::cout << "min normal " << normals[index*3] << " " 
                 //   <<normals[index*3+1] << " " << normals[index*3+2] << std::endl;
               //closestBoundaryIndex = i;
          closestBoundaryIndex = i;
          minIndex = i;
         }
  }

    if(isnan(directionUnitVector[0]) || isnan(directionUnitVector[1]) || isnan(directionUnitVector[2])){
	    //printf("minDist %f \n", minDistance);
	    //printf("directionV %f %f %f \n", directionUnitVector[0],directionUnitVector[1],directionUnitVector[2]);
	    directionUnitVector[0] = 0.0;
	    directionUnitVector[1] = 0.0;
	    directionUnitVector[2] = 0.0;
    }
      //vectorScalarMult(-1.0,directionUnitVector,directionUnitVector);
      //std::cout << "min dist " << minDistance << std::endl;
    }
    else
    {
                
    int direction_type;
    gitr_precision tol = 1e12;
    gitr_precision point1_dist;
    gitr_precision point2_dist;
    gitr_precision perp_dist;
    gitr_precision vectorMagnitude;
    gitr_precision max = 0.0;
    gitr_precision min = 0.0;
    gitr_precision Bfabsfperp = 0.0;
    gitr_precision distanceToParticle = 0.0;
    int pointLine=0;
    gitr_precision x;
     if( cylsymm > 0 )
     {
    x = std::sqrt(x0*x0 + y*y);
    }
    else
    {
    x = x0;
    }

    int top_limit = -1;
    gitr_precision dr;
    gitr_precision dz;

    int rInd;
    int zInd;

  if( geom_hash_sheath > 0 )
  {
  dr = closeGeomGridr[1] - closeGeomGridr[0];

  dz = closeGeomGridz[1] - closeGeomGridz[0];

  rInd = std::floor((x - closeGeomGridr[0])/dr + 0.5);

  zInd = std::floor((z - closeGeomGridz[0])/dz + 0.5);

  if(rInd >= nR_closeGeom) rInd = nR_closeGeom -1;

  if(zInd >= nZ_closeGeom) zInd = nZ_closeGeom -1;

  if(rInd < 0) rInd = 0;

  if(zInd < 0) zInd = 0;

  top_limit = n_closeGeomElements;
  }

  else top_limit = nLines;
  
  for( int k = 0; k < top_limit; k++) //n_closeGeomElements
    {
      int j = -1;

      if( geom_hash_sheath > 0 )
       j = closeGeom[zInd*nR_closeGeom*n_closeGeomElements + rInd*n_closeGeomElements + k];

      else j = k;

       gitr_precision boundZhere = boundaryVector[j].Z;
       
        if (boundZhere != 0.0)
        {
            point1_dist = std::sqrt((x - boundaryVector[j].x1)*(x - boundaryVector[j].x1) + 
                    (z - boundaryVector[j].z1)*(z - boundaryVector[j].z1));
            point2_dist = std::sqrt((x - boundaryVector[j].x2)*(x - boundaryVector[j].x2) + 
                                        (z - boundaryVector[j].z2)*(z - boundaryVector[j].z2));
            perp_dist = (boundaryVector[j].slope_dzdx*x - z + boundaryVector[j].intercept_z)/
                std::sqrt(boundaryVector[j].slope_dzdx*boundaryVector[j].slope_dzdx + 1.0);   
	
	
          if (std::abs(boundaryVector[j].slope_dzdx) >= tol*0.75)
	  {
	   perp_dist = x0 - boundaryVector[j].x1;
	  }
            if (point1_dist > point2_dist)
            {
                max = point1_dist;
                min = point2_dist;
            }
            else
            {
                max = point2_dist;
                min = point1_dist;
            }
            if (boundaryVector[j].length*boundaryVector[j].length + perp_dist*perp_dist >=
                    max*max)
            {
                distanceToParticle = std::abs(perp_dist);
                pointLine = 1;
            }
            else
            {
                distanceToParticle = min;
                if (boundaryVector[j].distanceToParticle == point1_dist)
                {
                    pointLine = 2;
                }
                else
                {
                    pointLine = 3;
                }
            }

            if (distanceToParticle < minDistance)
            {
                minDistance = distanceToParticle;
                minIndex = j;
                closestBoundaryIndex = j;
                direction_type = pointLine;
            }
        }
        else
        {
            distanceToParticle = tol;
        }
    }
    if (direction_type == 1)
    {
        if (boundaryVector[minIndex].slope_dzdx == 0)
        {
            directionUnitVector[0] = 0.0;
            directionUnitVector[1] = 0.0;
            directionUnitVector[2] = 1.0 * std::copysign(1.0,boundaryVector[minIndex].z1 - z);
        }
        else if (std::abs(boundaryVector[minIndex].slope_dzdx)>= 0.75*tol)
        {
            
            directionUnitVector[0] = boundaryVector[minIndex].x1 - x;
            directionUnitVector[1] = 0.0;
            directionUnitVector[2] = 0.0;
        }
        else
        {
            directionUnitVector[0] = 1.0 * std::copysign(1.0,(z - boundaryVector[minIndex].intercept_z)/(boundaryVector[minIndex].slope_dzdx) - x0);
            directionUnitVector[1] = 0.0;
            directionUnitVector[2] = 1.0 * std::copysign(1.0,perp_dist)/(boundaryVector[minIndex].slope_dzdx);
        }
    }
    else if (direction_type == 2)
    {
        directionUnitVector[0] = (boundaryVector[minIndex].x1 - x);
        directionUnitVector[1] = 0.0;
        directionUnitVector[2] = (boundaryVector[minIndex].z1 - z);
    }
    else
    {
        directionUnitVector[0] = (boundaryVector[minIndex].x2 - x);
        directionUnitVector[1] = 0.0;
        directionUnitVector[2] = (boundaryVector[minIndex].z2 - z);
    }

    vectorMagnitude = std::sqrt(directionUnitVector[0]*directionUnitVector[0] + directionUnitVector[1]*directionUnitVector[1]
                                + directionUnitVector[2]*directionUnitVector[2]);
    directionUnitVector[0] = directionUnitVector[0]/vectorMagnitude;
    directionUnitVector[1] = directionUnitVector[1]/vectorMagnitude;
    directionUnitVector[2] = directionUnitVector[2]/vectorMagnitude;
    }
    if( biased_surface > 0 )
    {
    pot = boundaryVector[minIndex].potential;
    Emag = pot/(2.0*boundaryVector[minIndex].ChildLangmuirDist)*exp(-minDistance/(2.0*boundaryVector[minIndex].ChildLangmuirDist));
    }
    else
    {
    angle = boundaryVector[minIndex].angle;    
    fd  = boundaryVector[minIndex].fd;
    pot = boundaryVector[minIndex].potential;
        gitr_precision debyeLength = boundaryVector[minIndex].debyeLength;
        gitr_precision larmorRadius = boundaryVector[minIndex].larmorRadius;
        Emag = pot*(fd/(2.0 * boundaryVector[minIndex].debyeLength)*exp(-minDistance/(2.0 * boundaryVector[minIndex].debyeLength))+ (1.0 - fd)/(boundaryVector[minIndex].larmorRadius)*exp(-minDistance/boundaryVector[minIndex].larmorRadius) );
        gitr_precision part1 = pot*(fd/(2.0 * boundaryVector[minIndex].debyeLength)*exp(-minDistance/(2.0 * boundaryVector[minIndex].debyeLength)));
        gitr_precision part2 = pot*(1.0 - fd)/(boundaryVector[minIndex].larmorRadius)*exp(-minDistance/boundaryVector[minIndex].larmorRadius);
    }
        /* Captain! This appears to be skipped? */
    if(minDistance == 0.0 || boundaryVector[minIndex].larmorRadius == 0.0)
    {
        Emag = 0.0;
        directionUnitVector[0] = 0.0;
        directionUnitVector[1] = 0.0;
        directionUnitVector[2] = 0.0;

    }
        
	Er = Emag*directionUnitVector[0];
        Et = Emag*directionUnitVector[1];
        E[2] = Emag*directionUnitVector[2];
        //std::cout << "Emag " << Emag << std::endl;
        //std::cout << "Min dist " << minDistance << std::endl;
        //std::cout << "r " << x << "z " << z << std::endl;
        //std::cout << "E components " << Er << " " << Et << " " << E[2] << std::endl;
        //std::cout << "direction unit vector " << directionUnitVector[0] << " " << directionUnitVector[1] << " " << directionUnitVector[2] << std::endl;
    
    //std::cout << "pos " << x << " " << y << " "<< z << " min Dist" << minDistance << "Efield " << Emag << std::endl;
    if( use_3d_geom > 0 )
    {
            E[0] = Er;
            E[1] = Et;
    }
    else
    {
     if( cylsymm > 0 )
     {
            //if cylindrical geometry
            gitr_precision theta = std::atan2(y,x0);
  
            E[0] = std::cos(theta)*Er - std::sin(theta)*Et;
            E[1] = std::sin(theta)*Er + std::cos(theta)*Et;
    }
    else
    {
            E[0] = Er;
            E[1] = Et;
    }
    }
            //std::cout << "Ex and Ey and Ez " << E[0] << " " << E[1] << " " << E[2] << std::endl;
   
      return minDistance;
}

move_boris::move_boris(

  Particles *_particlesPointer,
  gitr_precision _span,
  Boundary *_boundaryVector,
  int _nLines,
  int _nR_Bfield,
  int _nZ_Bfield,
  gitr_precision * _BfieldGridRDevicePointer,
  gitr_precision * _BfieldGridZDevicePointer,
  gitr_precision * _BfieldRDevicePointer,
  gitr_precision * _BfieldZDevicePointer,
  gitr_precision * _BfieldTDevicePointer,
  int _nR_Efield,
  int _nY_Efield,
  int _nZ_Efield,
  gitr_precision * _EfieldGridRDevicePointer,
  gitr_precision * _EfieldGridYDevicePointer,
  gitr_precision * _EfieldGridZDevicePointer,
  gitr_precision * _EfieldRDevicePointer,
  gitr_precision * _EfieldZDevicePointer,
  gitr_precision * _EfieldTDevicePointer,
  int _nR_closeGeom,
  int _nY_closeGeom,
  int _nZ_closeGeom,
  int _n_closeGeomElements,
  gitr_precision *_closeGeomGridr,
  gitr_precision *_closeGeomGridy,
  gitr_precision *_closeGeomGridz,
  int *_closeGeom,
  Flags* _gitr_flags,
  int sheath_efield_,
  int presheath_efield_,
  int biased_surface_,
  int geom_hash_sheath_,
  int use_3d_geom_,
  int cylsymm_,
  gitr_precision _max_dt)

  : 
  particlesPointer(_particlesPointer),
        boundaryVector(_boundaryVector),
        nR_Bfield(_nR_Bfield),
        nZ_Bfield(_nZ_Bfield),
        BfieldGridRDevicePointer(_BfieldGridRDevicePointer),
        BfieldGridZDevicePointer(_BfieldGridZDevicePointer),
        BfieldRDevicePointer(_BfieldRDevicePointer),
        BfieldZDevicePointer(_BfieldZDevicePointer),
        BfieldTDevicePointer(_BfieldTDevicePointer),
        nR_Efield(_nR_Efield),
        nY_Efield(_nY_Efield),
        nZ_Efield(_nZ_Efield),
        EfieldGridRDevicePointer(_EfieldGridRDevicePointer),
        EfieldGridYDevicePointer(_EfieldGridYDevicePointer),
        EfieldGridZDevicePointer(_EfieldGridZDevicePointer),
        EfieldRDevicePointer(_EfieldRDevicePointer),
        EfieldZDevicePointer(_EfieldZDevicePointer),
        EfieldTDevicePointer(_EfieldTDevicePointer),
        nR_closeGeom_sheath(_nR_closeGeom),
        nY_closeGeom_sheath(_nY_closeGeom),
        nZ_closeGeom_sheath(_nZ_closeGeom),
        n_closeGeomElements_sheath(_n_closeGeomElements),
        closeGeomGridr_sheath(_closeGeomGridr),
        closeGeomGridy_sheath(_closeGeomGridy),
        closeGeomGridz_sheath(_closeGeomGridz),
        closeGeom_sheath(_closeGeom),
	gitr_flags(_gitr_flags),
        max_dt(_max_dt),
        span(_span),
        nLines(_nLines),
        magneticForce{0.0, 0.0, 0.0},
        electricForce{0.0, 0.0, 0.0},
        sheath_efield( sheath_efield_ ),
        presheath_efield( presheath_efield_ ),
        biased_surface( biased_surface_ ),
        geom_hash_sheath( geom_hash_sheath_ ),
        use_3d_geom( use_3d_geom_ ),
        cylsymm( cylsymm_ )
        {}

CUDA_CALLABLE_MEMBER    
void move_boris::operator()(std::size_t indx)
{
  gitr_precision v_minus[3]= {0.0, 0.0, 0.0};
  gitr_precision v_prime[3]= {0.0, 0.0, 0.0};
  gitr_precision v_h[3]= {0.0, 0.0, 0.0};
  gitr_precision position0[3]= {0.0, 0.0, 0.0};
  gitr_precision position[3]= {0.0, 0.0, 0.0};
  gitr_precision v0[3]= {0.0, 0.0, 0.0};
  gitr_precision v_dt[3]= {0.0, 0.0, 0.0};
  gitr_precision v_half_dt[3]= {0.0, 0.0, 0.0};
  gitr_precision v[3]= {0.0, 0.0, 0.0};
  gitr_precision E[3] = {0.0, 0.0, 0.0};
  gitr_precision PSE[3] = {0.0, 0.0, 0.0};
  gitr_precision B[3] = {0.0,0.0,0.0};
  gitr_precision Bmag = 0.0;
  gitr_precision gyrofrequency = 0.0;
  gitr_precision q_prime = 0.0;
  gitr_precision coeff = 0.0;
  gitr_precision minDist = 0.0;
  int closestBoundaryIndex;

#if ODEINT ==	0 
  gitr_precision qpE[3] = {0.0,0.0,0.0};
  gitr_precision vmxB[3] = {0.0,0.0,0.0};
  gitr_precision vpxB[3] = {0.0,0.0,0.0};
  gitr_precision qp_vmxB[3] = {0.0,0.0,0.0};
  gitr_precision c_vpxB[3] = {0.0,0.0,0.0};
  vectorAssign(particlesPointer->xprevious[indx], particlesPointer->yprevious[indx], particlesPointer->zprevious[indx],position);
  vectorAssign(particlesPointer->xprevious[indx], particlesPointer->yprevious[indx], particlesPointer->zprevious[indx],position0);
  vectorAssign(particlesPointer->vx[indx], particlesPointer->vy[indx], particlesPointer->vz[indx],v0);
            
  gitr_precision dt = particlesPointer->dt[indx];
  gitr_precision vMag_dt = 0.0;
  gitr_precision vMag_half_dt = 0.0;
  gitr_precision half_dt = 0.5*dt;
  gitr_precision new_dt = 0.0;
  gitr_precision new_advance = false;
  if( sheath_efield > 0 )
  {
  minDist = getE(position[0], position[1], position[2],
		  E,boundaryVector,nLines,nR_closeGeom_sheath,
                  nY_closeGeom_sheath,nZ_closeGeom_sheath,
                  n_closeGeomElements_sheath,closeGeomGridr_sheath,
                  closeGeomGridy_sheath,
                  closeGeomGridz_sheath,closeGeom_sheath, closestBoundaryIndex,
                  biased_surface, use_3d_geom, geom_hash_sheath, cylsymm  );
  }

  if( presheath_efield > 0 )
  {
/*
#if LC_INTERP==3
              
	        //gitr_precision PSE2[3] = {0.0, 0.0, 0.0};
                 interp3dVector(PSE,position[0], position[1], position[2],nR_Efield,nY_Efield,nZ_Efield,
                     EfieldGridRDevicePointer,EfieldGridYDevicePointer,EfieldGridZDevicePointer,EfieldRDevicePointer,
                     EfieldZDevicePointer,EfieldTDevicePointer);
//E[0]= E[0] + PSE[0];
//E[1]= E[1] + PSE[1];
//E[2]= E[2] + PSE[2];
                 vectorAdd(E,PSE,E);
              //gitr_precision a = interp3d(position[0], position[1], position[2],nR_Efield,nY_Efield,nZ_Efield,
                //                EfieldGridRDevicePointer,EfieldGridYDevicePointer,EfieldGridZDevicePointer,EfieldRDevicePointer);
              //PSE[0] = 1.23;

#else
*/
  interp2dVector(&PSE[0],position[0], position[1], position[2],nR_Efield,nZ_Efield,
                     EfieldGridRDevicePointer,EfieldGridZDevicePointer,EfieldRDevicePointer,
                     EfieldZDevicePointer,EfieldTDevicePointer, cylsymm );
                 
  vectorAdd(E,PSE,E);
              //std::cout << "Efield in boris " <<E[0] << " " << E[1] << " " <<  E[2] << std::endl;
  }
  interp2dVector(&B[0],position[0], position[1], position[2],nR_Bfield,nZ_Bfield,
                    BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,
                    BfieldZDevicePointer,BfieldTDevicePointer, cylsymm );        
  Bmag = vectorNorm(B);
  gyrofrequency = particlesPointer->charge[indx]*1.60217662e-19*Bmag/(particlesPointer->amu[indx]*1.6737236e-27);

  //q_prime = 9.572528104401468e7*particlesPointer->charge[indx] / particlesPointer->amu[indx] * dt * 0.5;
  /* Captain! original code above, new code below. q_prime = q * dt / ( 2 * m ) */
  q_prime = particlesPointer->charge[ indx ] * gitr_constants::electron_volt * dt * 0.5 /
            ( particlesPointer->amu[ indx ] * gitr_constants::dalton );

    coeff = 2.0*q_prime/(1.0+(q_prime*Bmag)*(q_prime*Bmag));

    vectorAssign(particlesPointer->vx[indx], particlesPointer->vy[indx], particlesPointer->vz[indx],v);

        //vectorScalarMult(q_prime,E,qpE);
    qpE[0] = q_prime*E[0];
    qpE[1] = q_prime*E[1];
    qpE[2] = q_prime*E[2];
    //v_minus = v + q_prime*E
    //vectorAdd(v,qpE,v_minus);
     v_minus[0] = v[0] + qpE[0];
     v_minus[1] = v[1] + qpE[1];
     v_minus[2] = v[2] + qpE[2];
    //this->electricForce[0] = 2.0*qpE[0];
    //this->electricForce[1] = 2.0*qpE[1];
    //this->electricForce[2] = 2.0*qpE[2];
    
    //v_prime = v_minus + q_prime*(v_minus x B)
    vectorCrossProduct(v_minus,B,vmxB);
    vectorScalarMult(q_prime,vmxB,qp_vmxB);
    vectorAdd(v_minus,qp_vmxB,v_prime);       
    //this->magneticForce[0] = qp_vmxB[0];
    //this->magneticForce[1] = qp_vmxB[1];
    //this->magneticForce[2] = qp_vmxB[2];
     
     //v = v_minus + coeff*(v_prime x B)
     vectorCrossProduct(v_prime, B, vpxB);
     vectorScalarMult(coeff,vpxB,c_vpxB);
     vectorAdd(v_minus, c_vpxB, v_h);
     
     //v = v + q_prime*E
     //vectorAdd(v,qpE,v);
     v[0] = v_h[0] + qpE[0];
     v[1] = v_h[1] + qpE[1];
     v[2] = v_h[2] + qpE[2];
  //if(indx == 0){
  //printf("qprime qpexyz %.16e %.16e %.16e %.16e ",q_prime, qpE[0], qpE[1],qpE[2]);
  //////printf("vhxyz %.16e %.16e %.16e ",v_h[0], v_h[1], v_h[2]);
  //////printf("vmxyz %.16e %.16e %.16e ",v_minus[0], v_minus[1], v_minus[2]);
  ////printf("xyz %.16e %.16e %.16e ",position0[0], position0[1], position0[2]);
  ////printf("vxyz %.16e %.16e %.16e ",v[0], v[1], v[2]);
  ////printf("vmxyz %.16e %.16e %.16e ",v_minus[0], v_minus[1], v_minus[2]);
  ////printf("vpxyz %.16e %.16e %.16e ",v_prime[0], v_prime[1], v_prime[2]);
  ////printf("vpxBxyz %.16e %.16e %.16e ",vpxB[0], vpxB[1], vpxB[2]);
  ////printf("c_vpxBxyz %.16e %.16e %.16e \n",c_vpxB[0], c_vpxB[1], c_vpxB[2]);
  //}
	       
  if(gitr_flags->USE_ADAPTIVE_DT)
  {
    vectorAssign(v[0],v[1],v[2],v_dt);
    vMag_dt = vectorNorm(v_dt);
    
    vectorAssign(v0[0],v0[1],v0[2],v);
   
    // First step of half_dt
    q_prime = particlesPointer->charge[indx]*1.60217662e-19/(particlesPointer->amu[indx]*1.6737236e-27)*half_dt*0.5;
    coeff = 2.0*q_prime/(1.0+(q_prime*Bmag)*(q_prime*Bmag));
    
    //v = v + q_prime*E
    vectorScalarMult(q_prime,E,qpE);
    vectorAdd(v,qpE,v_minus);
    
    //v_prime = v_minus + q_prime*(v_minus x B)
    vectorCrossProduct(v_minus,B,vmxB);
    vectorScalarMult(q_prime,vmxB,qp_vmxB);
    vectorAdd(v_minus,qp_vmxB,v_prime);       
     
     //v = v_minus + coeff*(v_prime x B)
     vectorCrossProduct(v_prime, B, vpxB);
     vectorScalarMult(coeff,vpxB,c_vpxB);
     vectorAdd(v_minus, c_vpxB, v);
     
     //v = v + q_prime*E
     vectorAdd(v,qpE,v);
     
     position[0] = position[0] + v[0] * half_dt;
     position[1] = position[1] + v[1] * half_dt;
     position[2] = position[2] + v[2] * half_dt;
    
     // second step of half_dt
  if( sheath_efield > 0 )
  {
  minDist = getE(position[0], position[1], position[2],
		  E,boundaryVector,nLines,nR_closeGeom_sheath,
                  nY_closeGeom_sheath,nZ_closeGeom_sheath,
                  n_closeGeomElements_sheath,closeGeomGridr_sheath,
                  closeGeomGridy_sheath,
                  closeGeomGridz_sheath,closeGeom_sheath, closestBoundaryIndex,
                  biased_surface, use_3d_geom, geom_hash_sheath, cylsymm  );
  }

  if( presheath_efield > 0 )
  {
  interp2dVector(&PSE[0],position[0], position[1], position[2],nR_Efield,nZ_Efield,
                     EfieldGridRDevicePointer,EfieldGridZDevicePointer,EfieldRDevicePointer,
                     EfieldZDevicePointer,EfieldTDevicePointer, cylsymm );
                 
  vectorAdd(E,PSE,E);
  }
  interp2dVector(&B[0],position[0], position[1], position[2],nR_Bfield,nZ_Bfield,
                    BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,
                    BfieldZDevicePointer,BfieldTDevicePointer, cylsymm );        
  Bmag = vectorNorm(B);
    q_prime = particlesPointer->charge[indx]*1.60217662e-19/(particlesPointer->amu[indx]*1.6737236e-27)*half_dt*0.5;
    coeff = 2.0*q_prime/(1.0+(q_prime*Bmag)*(q_prime*Bmag));
    
    //v = v + q_prime*E
    vectorScalarMult(q_prime,E,qpE);
    vectorAdd(v,qpE,v_minus);
    
    //v_prime = v_minus + q_prime*(v_minus x B)
    vectorCrossProduct(v_minus,B,vmxB);
    vectorScalarMult(q_prime,vmxB,qp_vmxB);
    vectorAdd(v_minus,qp_vmxB,v_prime);       
     
     //v = v_minus + coeff*(v_prime x B)
     vectorCrossProduct(v_prime, B, vpxB);
     vectorScalarMult(coeff,vpxB,c_vpxB);
     vectorAdd(v_minus, c_vpxB, v);
     
     //v = v + q_prime*E
     vectorAdd(v,qpE,v);
     
     position[0] = position[0] + v[0] * half_dt;
     position[1] = position[1] + v[1] * half_dt;
     position[2] = position[2] + v[2] * half_dt;
    
     vMag_half_dt = vectorNorm(v);
     gitr_precision error_tolerance = 1.0e-7;
     gitr_precision velocity_error = std::abs(vMag_dt - vMag_half_dt)/vMag_half_dt;
     if(velocity_error <= error_tolerance)
     {
       new_advance = true;
     }
     else
     {
       new_advance = false;
     }

     gitr_precision max_term =
     std::max(static_cast<double>(std::sqrt(error_tolerance/(2*velocity_error))), 0.3);

     gitr_precision min_term = std::min(static_cast<double>(max_term), 2.0);
     new_dt = 0.9*dt*min_term;
     if (new_dt< span)
     {
      new_dt = span;
      new_advance = true;
     }

     if ( new_dt > 0.34/gyrofrequency)
     {
      new_dt = 0.34/gyrofrequency;
      new_advance = true;
     }

     if ( new_dt > max_dt)
     {
      new_dt = max_dt;
     }
     
     if (particlesPointer->charge[indx]==0)
     {
      new_dt = span;
      new_advance = true;
     }

     if((new_advance == true) && (particlesPointer->hitWall[indx] == 0.0))
     {
          particlesPointer->x[indx] = position[0];
          particlesPointer->y[indx] = position[1];
          particlesPointer->z[indx] = position[2];
          particlesPointer->vx[indx] = v[0];
          particlesPointer->vy[indx] = v[1];
          particlesPointer->vz[indx] = v[2];    
          particlesPointer->time[indx] = particlesPointer->time[indx]+dt;    
          particlesPointer->dt[indx] = new_dt;    
          particlesPointer->advance[indx] = new_advance;    
     }
     else if((new_advance == false) && (particlesPointer->hitWall[indx] == 0.0))
     {
          particlesPointer->dt[indx] = new_dt;    
          particlesPointer->advance[indx] = new_advance;    
     }
  }
  else
  {
    
     if(particlesPointer->hitWall[indx] == 0.0)
      {
          particlesPointer->x[indx] = position[0] + v[0] * dt;
          particlesPointer->y[indx] = position[1] + v[1] * dt;
          particlesPointer->z[indx] = position[2] + v[2] * dt;
          particlesPointer->vx[indx] = v[0];
          particlesPointer->vy[indx] = v[1];
          particlesPointer->vz[indx] = v[2];    
      }
  }

#endif

#if ODEINT == 1
        gitr_precision m = particlesPointer->amu[indx]*1.6737236e-27;
        gitr_precision q_m = particlesPointer->charge[indx]*1.60217662e-19/m;
        gitr_precision r[3]= {0.0, 0.0, 0.0};
        gitr_precision r2[3]= {0.0, 0.0, 0.0};
        gitr_precision r3[3]= {0.0, 0.0, 0.0};
        gitr_precision r4[3]= {0.0, 0.0, 0.0};
        gitr_precision v2[3]= {0.0, 0.0, 0.0};
        gitr_precision v3[3]= {0.0, 0.0, 0.0};
        gitr_precision v4[3]= {0.0, 0.0, 0.0};
        gitr_precision k1r[3]= {0.0, 0.0, 0.0};
        gitr_precision k2r[3]= {0.0, 0.0, 0.0};
        gitr_precision k3r[3]= {0.0, 0.0, 0.0};
        gitr_precision k4r[3]= {0.0, 0.0, 0.0};
        gitr_precision k1v[3]= {0.0, 0.0, 0.0};
        gitr_precision k2v[3]= {0.0, 0.0, 0.0};
        gitr_precision k3v[3]= {0.0, 0.0, 0.0};
        gitr_precision k4v[3]= {0.0, 0.0, 0.0};
        gitr_precision dtqm = dt*q_m;
        gitr_precision vxB[3] = {0.0,0.0,0.0};
        gitr_precision EplusvxB[3] = {0.0,0.0,0.0};
        gitr_precision halfKr[3] = {0.0,0.0,0.0};
        gitr_precision halfKv[3] = {0.0,0.0,0.0};
        gitr_precision half = 0.5;
                v[0] = particlesPointer->vx[indx];
                v[1] = particlesPointer->vy[indx];
	              v[2] = particlesPointer->vz[indx];

                r[0] = particlesPointer->xprevious[indx];
                r[1] = particlesPointer->yprevious[indx];
	              r[2] = particlesPointer->zprevious[indx];
#ifdef __CUDACC__
#else
#endif
for ( int s=0; s<nSteps; s++ ) 
    {
#ifdef __CUDACC__
#else
#endif
    if( sheath_efield > 0 )
    {
    minDist = getE(r[0],r[1],r[2],E,boundaryVector,nLines);
    }

    if( presheath_efield > 0 )
    {
    interparticlesPointer->dVector(&particlesPointer->E[0],particlesPointer->xparticlesPointer->evious,particlesPointer->yparticlesPointer->evious,particlesPointer->zparticlesPointer->evious,nR_Efield,nZ_Efield,
          EfieldGridRDeviceparticlesPointer->inter,EfieldGridZDeviceparticlesPointer->inter,EfieldRDeviceparticlesPointer->inter,
          EfieldZDeviceparticlesPointer->inter,EfieldTDeviceparticlesPointer->inter);
                 
    vectorAdd(E,particlesPointer->E,E);
    }
#ifdef __CUDACC__
#else
#endif
    interp2dVector(&B[0],r[0],r[1],r[2],nR_Bfield,nZ_Bfield,
               BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,
               BfieldZDevicePointer,BfieldTDevicePointer);        
#ifdef __CUDACC__
#else
#endif
    //k1r = dt*v
    vectorScalarMult(dt,v,k1r);
    /*
    k1r[0] = v[0]*dt;
    k1r[1] = v[1]*dt;
    k1r[2] = v[2]*dt;
    */
    //k1v = dt*q_m * (E + (v x B))
    vectorCrossProduct(v,B,vxB);
    vectorAdd(E,vxB,EplusvxB);
    vectorScalarMult(dtqm,EplusvxB,k1v);
    /*
    k1v[0] = dt*q_m*(E[0] + (v[1]*B[2] - v[2]*B[1]));
    k1v[1] = dt*q_m*(E[1] + (v[2]*B[0] - v[0]*B[2]));
    k1v[2] = dt*q_m*(E[2] + (v[0]*B[1] - v[1]*B[0]));
    */
    //r2 = r + 0.5*k1r
    vectorScalarMult(half,k1r,halfKr);
    vectorAdd(r,k1r,r2);
    /*
    r2[0] = r[0] + k1r[0]*0.5;
    r2[1] = r[1] + k1r[1]*0.5;
    r2[2] = r[2] + k1r[2]*0.5;
    */

    //v2 = v + 0.5*k1v
    vectorScalarMult(half,k1v,halfKv);
    vectorAdd(v, halfKv,v2);
        /*
    v2[0] = v[0] + k1v[0]*0.5;
    v2[1] = v[1] + k1v[1]*0.5;
    v2[2] = v[2] + k1v[2]*0.5;
    */
#ifdef __CUDACC__
#else
#endif

if( sheath_efield > 0 )
{
    minDist = getE(r2[0],r2[1],r2[2],E,boundaryVector,nLines);
}
if( presheath_efield > 0 )
{
    interparticlesPointer->dVector(&particlesPointer->E[0],particlesPointer->xparticlesPointer->evious,particlesPointer->yparticlesPointer->evious,particlesPointer->zparticlesPointer->evious,nR_Efield,nZ_Efield,
               EfieldGridRDeviceparticlesPointer->inter,EfieldGridZDeviceparticlesPointer->inter,EfieldRDeviceparticlesPointer->inter,
               EfieldZDeviceparticlesPointer->inter,EfieldTDeviceparticlesPointer->inter);
    vectorAdd(E,particlesPointer->E,E);
}
#ifdef __CUDACC__
#else
#endif


    interp2dVector(&B[0],r2[0],r2[1],r2[2],nR_Bfield,nZ_Bfield,
             BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,
             BfieldZDevicePointer,BfieldTDevicePointer);        
#ifdef __CUDACC__
#else
#endif
    //k2r = dt*v2
    vectorScalarMult(dt,v2,k2r);
    /*
    k2r[0] = v2[0]*dt;
    k2r[1] = v2[1]*dt;
    k2r[2] = v2[2]*dt;
    */
    //k2v = dt*q_m*(E + (v x B))
    vectorCrossProduct(v2,B,vxB);
    vectorAdd(E,vxB,EplusvxB);
    vectorScalarMult(dtqm,EplusvxB,k2v);
    /*
    k2v[0] = dt*q_m*(E[0] + (v2[1]*B[2] - v2[2]*B[1]));
    k2v[1] = dt*q_m*(E[1] + (v2[2]*B[0] - v2[0]*B[2]));
    k2v[2] = dt*q_m*(E[2] + (v2[0]*B[1] - v2[1]*B[0]));
    */
    //r3 = r + 0.5*k2r
    vectorScalarMult(half,k2r,halfKr);
    vectorAdd(r,k2r,r3);
    /*
    r3[0] = r[0] + k2r[0]*0.5;
    r3[1] = r[1] + k2r[1]*0.5;
    r3[2] = r[2] + k2r[2]*0.5;
    */
    //v3 = v + 0.5*k2v
    vectorScalarMult(half,k2v,halfKv);
    vectorAdd(v, halfKv,v3);
    /*
    v3[0] = v[0] + k2v[0]*0.5;
    v3[1] = v[1] + k2v[1]*0.5;
    v3[2] = v[2] + k2v[2]*0.5;
    */
#ifdef __CUDACC__
#else
#endif

if( sheath_efield > 0 )
{
    minDist = getE(r3[0],r3[1],r3[2],E,boundaryVector,nLines);
}

if( presheath_efield > 0 )
{
    interparticlesPointer->dVector(&particlesPointer->E[0],particlesPointer->xparticlesPointer->evious,particlesPointer->yparticlesPointer->evious,particlesPointer->zparticlesPointer->evious,nR_Efield,nZ_Efield,
               EfieldGridRDeviceparticlesPointer->inter,EfieldGridZDeviceparticlesPointer->inter,EfieldRDeviceparticlesPointer->inter,
               EfieldZDeviceparticlesPointer->inter,EfieldTDeviceparticlesPointer->inter);
    vectorAdd(E,particlesPointer->E,E);
}

#ifdef __CUDACC__
#else
#endif
    interp2dVector(&B[0],r3[0],r3[1],r3[2],nR_Bfield,nZ_Bfield,
                 BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,
                 BfieldZDevicePointer,BfieldTDevicePointer);        
                
#ifdef __CUDACC__
#else
#endif
    //k3r = dt*v3
    vectorScalarMult(dt,v3,k3r);
    /*
    k3r[0] = v3[0]*dt;
    k3r[1] = v3[1]*dt;
    k3r[2] = v3[2]*dt;
    */
    //k3v = dt*qm*(E + (v x B))
    vectorCrossProduct(v3,B,vxB);
    vectorAdd(E,vxB,EplusvxB);
    vectorScalarMult(dtqm,EplusvxB,k3v);
    /*
    k3v[0] = dt*q_m*(E[0] + (v3[1]*B[2] - v3[2]*B[1]));
    k3v[1] = dt*q_m*(E[1] + (v3[2]*B[0] - v3[0]*B[2]));
    k3v[2] = dt*q_m*(E[2] + (v3[0]*B[1] - v3[1]*B[0]));
    */
    //r4 = r + k3r
    vectorAdd(r, k3r,r4);
    /*
    r4[0] = r[0] + k3r[0];
    r4[1] = r[1] + k3r[1];
    r4[2] = r[2] + k3r[2];
    */
    //v4 = v + k3v
    vectorAdd(v, k3v, v4);
        /*
    v4[0] = v[0] + k3v[0];
    v4[1] = v[1] + k3v[1];
    v4[2] = v[2] + k3v[2];
    */
#ifdef __CUDACC__
#else
#endif

  if( sheath_efield > 0 )
  {
	minDist = getE(r4[0],r4[1],r4[2],E,boundaryVector,nLines);
  }
  if( presheath_efield > 0 )
  {
   interp2dVector(&particlesPointer->E[0],particlesPointer->xparticlesPointer->evious,particlesPointer->yparticlesPointer->evious,particlesPointer->zparticlesPointer->evious,nR_Efield,nZ_Efield,
               EfieldGridRDeviceparticlesPointer->inter,EfieldGridZDeviceparticlesPointer->inter,EfieldRDeviceparticlesPointer->inter,
               EfieldZDeviceparticlesPointer->inter,EfieldTDeviceparticlesPointer->inter);
    vectorAdd(E,particlesPointer->E,E);
  }
#ifdef __CUDACC__
#else
#endif

    interp2dVector(&B[0],r4[0],r4[1],r4[2],nR_Bfield,nZ_Bfield,
                        BfieldGridRDevicePointer,BfieldGridZDevicePointer,
                        BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer);        
#ifdef __CUDACC__
#else
#endif

    //k4r = dt*v4
    vectorScalarMult(dt,v4,k4r);
    /*
   k4r[0] = v4[0]*dt;
   k4r[1] = v4[1]*dt;
   k4r[2] = v4[2]*dt;
   */
    //k4v = dt*q_m*(E + (v x B))
    vectorCrossProduct(v4,B,vxB);
    vectorAdd(E,vxB,EplusvxB);
    vectorScalarMult(dtqm,EplusvxB,k4v);
    /*
   k4v[0] = dt*q_m*(E[0] + (v4[1]*B[2] - v4[2]*B[1]));
   k4v[1] = dt*q_m*(E[1] + (v4[2]*B[0] - v4[0]*B[2]));
   k4v[2] = dt*q_m*(E[2] + (v4[0]*B[1] - v4[1]*B[0]));
   */
   particlesPointer->x[indx] = r[0] + (k1r[0] + 2*k2r[0] + 2*k3r[0] + k4r[0])/6;
   particlesPointer->y[indx] = r[1] + (k1r[1] + 2*k2r[1] + 2*k3r[1] + k4r[1])/6;
   particlesPointer->z[indx] = r[2] + (k1r[2] + 2*k2r[2] + 2*k3r[2] + k4r[2])/6;
   particlesPointer->vx[indx] = v[0] + (k1v[0] + 2*k2v[0] + 2*k3v[0] + k4v[0])/6;
   particlesPointer->vy[indx] = v[1] + (k1v[1] + 2*k2v[1] + 2*k3v[1] + k4v[1])/6;
   particlesPointer->vz[indx] = v[2] + (k1v[2] + 2*k2v[2] + 2*k3v[2] + k4v[2])/6;
#ifdef __CUDACC__
#else
#endif
//std::cout << "OparticlesPointer->rations Time: " << oparticlesPointer->rationsTime <<std::endl;
//std::cout << "Efield InterparticlesPointer->lation Time: " << interparticlesPointer->Time <<std::endl;
//std::cout << "Bfield InterparticlesPointer->lation Time: " << interparticlesPointer->Time <<std::endl;
//std::cout << "Init Time: " << initTime <<std::endl;
            }
#endif
} 
