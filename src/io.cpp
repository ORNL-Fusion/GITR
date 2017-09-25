#include "h1.cuh"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <netcdf>
#include "Boundary.h"
#include "Particle.h"
#include "libconfig.h++"
#if USE_BOOST
	#include "boost/multi_array.hpp"
	#include "boost/filesystem.hpp"
#endif
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

//void handle_error(int status) {
//    if (status != NC_NOERR) {
//           fprintf(stderr, "%s\n", nc_strerror(status));
//              exit(-1);
//                 }
//}
using namespace std;
using namespace netCDF;
using namespace exceptions;
using namespace netCDF::exceptions;
using namespace libconfig;
//IO
int importGeometry(libconfig::Config &cfg_geom, sim::Array<Boundary> &boundaries)
{
    Setting& geom = cfg_geom.lookup("geom");
    std::cout << "Boundary import routine " << int(boundaries.size()) << std::endl;
    int nLines = boundaries.size() - 1;
  std::string geom_outname = "geom.m";
  std::string geom_folder = "output/geometry";
  ofstream outfile;

  #if USE_BOOST
    boost::filesystem::path dirOut("output");
    
    if(!(boost::filesystem::exists(dirOut)))
    {
       std::cout<<"Doesn't Exists"<<std::endl;

       if (boost::filesystem::create_directory(dirOut))
       {
          std::cout << " Successfully Created " << std::endl;
       }
    }
    //Output
    boost::filesystem::path dir(geom_folder);
    
    if(!(boost::filesystem::exists(dir)))
    {
       std::cout<<"Doesn't Exists"<<std::endl;

       if (boost::filesystem::create_directory(dir))
       {
          std::cout << " Successfully Created " << std::endl;
       }
    }
  #endif

  std::string full_path = geom_folder + "/" + geom_outname;
  outfile.open (full_path );
  #if USE3DTETGEOM > 0
    for(int i=0 ; i<nLines ; i++)
    {
       boundaries[i].x1 = geom["x1"][i];
       boundaries[i].y1 = geom["y1"][i];
       boundaries[i].z1 = geom["z1"][i];
       boundaries[i].x2 = geom["x2"][i];
       boundaries[i].y2 = geom["y2"][i];
       boundaries[i].z2 = geom["z2"][i];
       boundaries[i].x3 = geom["x3"][i];
       boundaries[i].y3 = geom["y3"][i];
       boundaries[i].z3 = geom["z3"][i];
       boundaries[i].Z = geom["Z"][i];
       boundaries[i].a = geom["a"][i];
       boundaries[i].b = geom["b"][i];
       boundaries[i].c = geom["c"][i];
       boundaries[i].d = geom["d"][i];
       boundaries[i].plane_norm = geom["plane_norm"][i];
       boundaries[i].area = geom["area"][i];
     }   

    outfile.close();
  #else

    //int nMaterials = geom["nMaterials"];
    //std::cout << "nmat " << nMaterials << std::endl;
    for(int i=0 ; i<nLines ; i++)
    {
       boundaries[i].x1 = geom["x1"][i];
       boundaries[i].z1 = geom["z1"][i];
       boundaries[i].x2 = geom["x2"][i];
       boundaries[i].z2 = geom["z2"][i];
       boundaries[i].Z = geom["Z"][i];
       boundaries[i].slope_dzdx = geom["slope"][i];
       boundaries[i].intercept_z = geom["intercept"][i];
       boundaries[i].length = geom["length"][i];

       outfile << "geom(" << i+1 << ",:) = ["<<boundaries[i].x1 << ", " <<
          boundaries[i].z1 << ", " <<
          boundaries[i].x2 << ", " << boundaries[i].z2 << ", " <<
          boundaries[i].slope_dzdx << ", " << boundaries[i].intercept_z << ", " <<
          boundaries[i].length << ", " << boundaries[i].Z << "];" << std::endl;
    }   

    outfile.close();
    boundaries[nLines].Z = geom["Z"][nLines];
    boundaries[nLines].y1 = geom["y1"];
    boundaries[nLines].y2 = geom["y2"];
    boundaries[nLines].periodic = geom["periodic"];
  #endif
    return 0;
}
int read_ar2Input( string fileName, float *Bfield[]) {

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

    for(int i=0; i<nR; i++){
        for(int j=0; j<nZ; j++){
           Bfield[i][j] = br[j][i]; 
        }
    }

    return(0);

}


int read_profileNs( string fileName, string nxName, string nzName,int &n_x,int &n_z ) {

    // Check input file exists

    ifstream file(fileName.c_str());
    if(!file.good()) {
        cout<<"ERROR: Cannot file input file ... "<<fileName<<endl;
        exit(1);
    }

    NcFile nc(fileName.c_str(), NcFile::read);

    NcDim nc_nx(nc.getDim(nxName));
    NcDim nc_nz(nc.getDim(nzName));
    
    n_x = nc_nx.getSize(); 
    n_z = nc_nz.getSize(); 


    return(0);

}

int read_profileNsChar(const char *fileName,const char *nxName,const char *nzName,int &n_x,int &n_z ) {

    // Check input file exists

    //ifstream file(fileName.c_str());
    ifstream file(fileName);
    if(!file.good()) {
        cout<<"ERROR: Cannot file input file ... "<<fileName<<endl;
        exit(1);
    }

    //NcFile nc(fileName.c_str(), NcFile::read);
    NcFile nc(fileName, NcFile::read);

    NcDim nc_nx(nc.getDim(nxName));
    NcDim nc_nz(nc.getDim(nzName));
    
    n_x = nc_nx.getSize(); 
    n_z = nc_nz.getSize(); 


    return(0);

}

int read_profiles( string fileName, int &n_x, int &n_z,string gridxName, sim::Array<float>& gridx,string gridzName,
          sim::Array<float>& gridz,string dataName, sim::Array<float>& data) {

    // Check input file exists

    ifstream file(fileName.c_str());
    if(!file.good()) {
        cout<<"ERROR: Cannot file input file ... "<<fileName<<endl;
        exit(1);
    }

    NcFile nc(fileName.c_str(), NcFile::read);

    NcVar nc_gridx(nc.getVar(gridxName));
    NcVar nc_gridz(nc.getVar(gridzName));

    nc_gridx.getVar(&gridx[0]);
    nc_gridz.getVar(&gridz[0]);

    NcVar nc_ne(nc.getVar(dataName));
    nc_ne.getVar(&data[0]);

    return(0);

}

int read_profile2d( string fileName,string dataName, sim::Array<float>& data) {
    std::cout << "reading 2d profile" << std::endl;
    //NcError err(2);
    //NcError::Behavior bb= (NcError::Behavior) 0;
    //NcError err(NcError::silent_nonfatal);
    // Check input file exists

    ifstream file(fileName.c_str());
    if(!file.good()) {
        cout<<"ERROR: Cannot file input file ... "<<fileName<<endl;
        exit(1);
    }

    NcFile nc(fileName.c_str(), NcFile::read);


    //NcVar nc_ne(nc.getVar(dataName));
    static const int NC_ERR = 2;
    NcVar nc_ne;
    try{
        std::cout << "inside try " << std::endl;
        nc_ne = nc.getVar(dataName);
    if(nc_ne.isNull()){std::cout << "erororororo" << std::endl;}
        std::cout << "finishing try " << std::endl;
    }
    catch(int an)//NcException& e
    {
        
        std::cout << " problem in file: " << fileName << " variable: " << dataName << std::endl;
        //e.what();
        std::cout << " problem in file: " << fileName << " variable: " << dataName << std::endl;
        //return NC_ERR;
    }
    nc_ne.getVar(&data[0]);

    return(0);

}
int read_profile3d( string fileName,string dataName, sim::Array<int>& data) {

    // Check input file exists

    ifstream file(fileName.c_str());
    if(!file.good()) {
        cout<<"ERROR: Cannot file input file ... "<<fileName<<endl;
        exit(1);
    }

    NcFile nc(fileName.c_str(), NcFile::read);


    NcVar nc_ne(nc.getVar(dataName));
    nc_ne.getVar(&data[0]);

    return(0);

}
int read_profile1d( string fileName,string gridxName, sim::Array<float>& gridx) {

    // Check input file exists

    ifstream file(fileName.c_str());
    if(!file.good()) {
        cout<<"ERROR: Cannot file input file ... "<<fileName<<endl;
        exit(1);
    }

    NcFile nc(fileName.c_str(), NcFile::read);

    NcVar nc_gridx(nc.getVar(gridxName));

    nc_gridx.getVar(&gridx[0]);


    return(0);

}
void OUTPUT(char outname[],int nX, int nY, float **array2d)
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

void OUTPUT2d(std::string folder,std::string outname,int nX, int nY, float *array2d)
{
       ofstream outfile;
#if USE_BOOST	
			//Output
        boost::filesystem::path dir(folder);

            if(!(boost::filesystem::exists(dir)))
            {
              std::cout<<"Doesn't Exists"<<std::endl;
              if (boost::filesystem::create_directory(dir))
              {
              std::cout << " Successfully Created " << std::endl;
              }
            }
#endif
            std::string full_path = folder + "/" + outname;
			outfile.open (full_path );
			
				 for(int i=1 ; i<=nY ; i++)
				{
				outfile << "val2d( :," << i<< ") = [ " ;
					for(int j=0 ; j<nX ; j++)
					{
					outfile << array2d[(i-1)*nX + j] << "  " ;
					//std::cout << r[i] << std::endl;
					}
					outfile << "  ];" << std::endl;
				}
			outfile.close();	
		
		
}

void OUTPUT1d(std::string folder,std::string outname,int nX, float *array2d)
{
       ofstream outfile;
#if USE_BOOST
				//Output
        boost::filesystem::path dir(folder);

            if(!(boost::filesystem::exists(dir)))
            {
             // std::cout<<"Doesn't Exists"<<std::endl;
              if (boost::filesystem::create_directory(dir))
              {
              //std::cout << " Successfully Created " << std::endl;
              }
            }
#endif
            std::string full_path = folder + "/" + outname;
			outfile.open (full_path );
			
				outfile << "val1d " << "  = [ " ;
				 for(int i=0 ; i<nX ; i++)
				{
					outfile << array2d[i] << "  " ;
				}
					outfile << "  ];" << std::endl;
			outfile.close();	
		
		
}

void OUTPUT3d(std::string folder,std::string outname,int nX, int nY, int nZ, float *array3d)
{
       ofstream outfile;
#if USE_BOOST	
			//Output
        boost::filesystem::path dir(folder);

            if(!(boost::filesystem::exists(dir)))
            {
              std::cout<<"Doesn't Exists"<<std::endl;
              if (boost::filesystem::create_directory(dir))
              {
              std::cout << " Successfully Created " << std::endl;
              }
            }
#endif
            std::string full_path = folder + "/" + outname;
			outfile.open (full_path );
			for(int k=1; k<=nZ; k++)
            {
				 for(int i=1 ; i<=nY ; i++)
				{
				outfile << "val3d( :," << i<< "," << k << ") = [ " ;
					for(int j=0 ; j<nX ; j++)
					{
					outfile << array3d[(k-1)*nX*nY + (i-1)*nX + j] << "  " ;
					//std::cout << r[i] << std::endl;
					}
					outfile << "  ];" << std::endl;
				}
            }
			outfile.close();	
		
		
}
void OUTPUT2d(std::string folder,std::string outname,int nX, int nY, int *array2d)
{
       ofstream outfile;
#if USE_BOOST	
			//Output
        boost::filesystem::path dir(folder);

            if(!(boost::filesystem::exists(dir)))
            {
              std::cout<<"Doesn't Exists"<<std::endl;
              if (boost::filesystem::create_directory(dir))
              {
              std::cout << " Successfully Created " << std::endl;
              }
            }
#endif
            std::string full_path = folder + "/" + outname;
			outfile.open (full_path );
			
				 for(int i=1 ; i<=nY ; i++)
				{
				outfile << "val2d( :," << i<< ") = [ " ;
					for(int j=0 ; j<nX ; j++)
					{
					outfile << array2d[(i-1)*nX + j] << "  " ;
					//std::cout << r[i] << std::endl;
					}
					outfile << "  ];" << std::endl;
				}
			outfile.close();	
		
		
}

void OUTPUT1d(std::string folder,std::string outname,int nX, int *array2d)
{
       ofstream outfile;
#if USE_BOOST
				//Output
        boost::filesystem::path dir(folder);

            if(!(boost::filesystem::exists(dir)))
            {
             // std::cout<<"Doesn't Exists"<<std::endl;
              if (boost::filesystem::create_directory(dir))
              {
              //std::cout << " Successfully Created " << std::endl;
              }
            }
#endif
            std::string full_path = folder + "/" + outname;
			outfile.open (full_path );
			
				outfile << "val1d " << "  = [ " ;
				 for(int i=0 ; i<nX ; i++)
				{
					outfile << array2d[i] << "  " ;
				}
					outfile << "  ];" << std::endl;
			outfile.close();	
		
		
}

void OUTPUT3d(std::string folder,std::string outname,int nX, int nY, int nZ, int *array3d)
{
       ofstream outfile;
#if USE_BOOST	
			//Output
        boost::filesystem::path dir(folder);

            if(!(boost::filesystem::exists(dir)))
            {
              std::cout<<"Doesn't Exists"<<std::endl;
              if (boost::filesystem::create_directory(dir))
              {
              std::cout << " Successfully Created " << std::endl;
              }
            }
#endif
            std::string full_path = folder + "/" + outname;
			outfile.open (full_path );
			for(int k=1; k<=nZ; k++)
            {
				 for(int i=1 ; i<=nY ; i++)
				{
				outfile << "val3d( :," << i<< "," << k << ") = [ " ;
					for(int j=0 ; j<nX ; j++)
					{
					outfile << array3d[(k-1)*nX*nY + (i-1)*nX + j] << "  " ;
					//std::cout << r[i] << std::endl;
					}
					outfile << "  ];" << std::endl;
				}
            }
			outfile.close();	
		
		
}
int readFileDim(const std::string& fileName,const std::string& varName)
{
           NcFile nc(fileName, NcFile::read);

                  if(nc.isNull()){
                             std::cout << "ERROR: Failed to open " << fileName << std::endl;
                                    }
                         NcDim nc_nx(nc.getDim(varName));

                                int n_x = nc_nx.getSize();
                                       return n_x;

}
int ncdfIO(int rwCode,const std::string& fileName,vector< std::string> dimNames,vector<int> dims,
        vector< std::string> gridNames,vector<int> gridMapToDims,
         vector<float*> pointers,vector< std::string> intVarNames,vector<vector<int>> intVarDimMap, vector<int*> intVarPointers)
{
    std::cout << "readWritecode " << rwCode << std::endl;//0 is read, 1 is write
    if(rwCode == 1)
    {
       NcFile thisFile;
       try{
             thisFile.open(fileName+".nc",NcFile::newFile);
             if(thisFile.isNull())
             {
                std::cout << "ERROR: Failed to open " << fileName << std::endl;
             }
             else
             {   
                std::cout << "successfully opened " << fileName << std::endl;
             }
          }
       catch(NcException& e)
       {
          std::cout << " caught exists exception " << thisFile.isNull() << std::endl;
          int fileInt = 0;
          while (thisFile.isNull())
          {
            std::cout << "filename " << fileName << " is taken " <<std::endl;
            try{
                  thisFile.open(fileName+std::to_string(fileInt)+".nc",NcFile::newFile);
               }
            catch(NcException& e){
            std::cout << "Filename " << fileName+std::to_string(fileInt)+".nc" <<
                                 " is taken " << std::endl;
            }
                         fileInt++;
          }
       }

       vector<NcDim> theseDims(dims.size());
       for(int i=0;i<theseDims.size();i++)
       {
           theseDims[i] = thisFile.addDim(dimNames[i],dims[i]);
       }

       vector<NcVar> theseVars(dims.size());
       for(int i=0;i<pointers.size();i++)
       {
           theseVars[i] = thisFile.addVar(gridNames[i],ncDouble,theseDims[gridMapToDims[i]]);
           theseVars[i].putVar(pointers[i]);
       }
       vector<NcVar> intVars(intVarNames.size());
       vector<vector<NcDim>> intVarDims(intVarNames.size());
       for(int i=0;i<intVarNames.size();i++)
       {
           for(int j=0;j<intVarDimMap[i].size();j++)
           {
               intVarDims[i].push_back(theseDims[intVarDimMap[i][j]]);
           }
           intVars[i] = thisFile.addVar(intVarNames[i],NC_INT,intVarDims[i]);
           intVars[i].putVar(intVarPointers[i]);
       }
      thisFile.close();
    }
    std::cout << "filename " << fileName << std::endl;

    for(int i=0; i<dimNames.size();i++)
    {
       std::cout << "dimname " <<i << " " <<  dimNames[i] << std::endl;
    }
    for(int i=0; i<dims.size();i++)
    {
       std::cout << "dimension " <<i << " " <<  dims[i] << std::endl;
    }
}
//template <class T>
//int readFileVar1d(const char *fileName,const char *varName,T &x ) {
//
//    // Check input file exists
//
//    ifstream file(fileName);
//    if(!file.good()) {
//        cout<<"ERROR: Cannot file input file ... "<<fileName<<endl;
//        exit(1);
//    }
//
//    NcFile nc(fileName, NcFile::read);
//
//    NcDim nc_xdim(nc.getDim(varName));
//    NcVar nc_x;
//
//    int n_x=0;
//   
//    try{
//         n_x = nc_xdim.getSize(); 
//         if(n_x==0){
//             printf("ERROR reading filename \n");
//             }
//         else
//         {
//            try{
//                  nc_x = nc.getVar(varName);
//                      if(nc_x.isNull()){std::cout << "ERROR:   stuff" << std::endl;}
//               }
//            catch(NcException& e)
//                 {}
//         }
//    }
//    catch(NcException& e){}
//
//    nc_x.getVar(&x[0]);
//    return n_x;
//
//}

