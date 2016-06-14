#include "h1.cuh"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <netcdf>
#include "Boundary.h"
#include "Particle.h"
#include "boost/multi_array.hpp"
#include "io.hpp"
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


    // Allocate contiguous 2D array for netcdf to work
    float **br = new float*[nR];
    br[0] = new float[nR*nZ];
    for(int i=0; i<nR; i++){
        br[i] = &br[0][i*nZ];
    }


    NcVar nc_br(nc.getVar("br"));

    nc_br.getVar(br[0]);

    return(0);

}


int read_profileNs( string fileName,int &n_x,int &n_z ) {

    // Check input file exists

    ifstream file(fileName.c_str());
    if(!file.good()) {
        cout<<"ERROR: Cannot file input file ... "<<fileName<<endl;
        exit(1);
    }

    NcFile nc(fileName.c_str(), NcFile::read);

    NcDim nc_nx(nc.getDim("n_x"));
    NcDim nc_nz(nc.getDim("n_z"));
    
    n_x = nc_nx.getSize(); 
    n_z = nc_nz.getSize(); 


    return(0);

}


int read_profiles( string fileName, int &n_x, int &n_z,std::vector<double>& gridx, 
                    std::vector<double>& gridz, std::vector<double>& data) {

    // Check input file exists

    ifstream file(fileName.c_str());
    if(!file.good()) {
        cout<<"ERROR: Cannot file input file ... "<<fileName<<endl;
        exit(1);
    }

    NcFile nc(fileName.c_str(), NcFile::read);

    NcVar nc_gridx(nc.getVar("gridx"));
    NcVar nc_gridz(nc.getVar("gridz"));

    nc_gridx.getVar(&gridx[0]);
    nc_gridz.getVar(&gridz[0]);

    NcVar nc_ne(nc.getVar("ne"));
    std::cout << " about to read data from nc object" << std::endl;
    nc_ne.getVar(&data[0]);

    return(0);

}
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
