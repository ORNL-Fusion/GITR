#include "h1.cuh"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <netcdf>
#include "Boundary.h"
#include "Particle.h"

#ifdef __CUDACC__
#include <thrust/host_vector.h>
#endif

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

using namespace std;
using namespace netCDF;
using namespace exceptions;

//IO

int read_ar2Input( string fileName) {

    cout << "Reading "<<fileName<<endl;

    // Check input file exists

    ifstream file(fileName.c_str());
    if(!file.good()) {
        cout<<"ERROR: Cannot file input file ... "<<fileName<<endl;
        exit(1);
    }

    NcFile nc(fileName.c_str(), NcFile::read);

    NcDim nc_nR(nc.getDim("nR"));
    NcDim nc_nZ(nc.getDim("nZ"));
    
    int nR = nc_nR.getSize(); 
    int nZ = nc_nZ.getSize(); 

    NcVar nc_r(nc.getVar("r"));
    NcVar nc_z(nc.getVar("z"));

    vector<float> r;
    r.resize(nR);
    nc_r.getVar(&r[0]);

    vector<float> z;
    z.resize(nZ);
    nc_z.getVar(&z[0]);

    for(int i=0;i<nR;i++){
        cout<<r[i]<<endl;
    }

    for(int j=0;j<nZ;j++){
        cout<<z[j]<<endl;
    }

    return(0);

}


//OUTPUT
void OUTPUT(char outname[],int nX, int nY, double **array2d)
{
       ofstream outfile;
				//Output


			outfile.open (outname );
			
				 for(int i=1 ; i<=nX ; i++)
				{
				outfile << "Dep( " << i<< ",:) = [ " ;
					for(int j=0 ; j<nY ; j++)
					{
					outfile << array2d[i-1][j] << "  " ;
					//std::cout << r[i] << std::endl;
					}
					outfile << "  ];" << std::endl;
				}
			outfile.close();	
		
		
}
