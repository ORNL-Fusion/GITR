#ifndef _IO_
#define _IO_

using namespace std;

int read_ar2Input( std::string fileName, float *Bfield[]);

int read_profileNs( std::string fileName,std::string nzName,std::string nxName,int &n_x,int &n_z );
int read_profile2d( string fileName,string dataName, sim::Array<float>& data);
int read_profile1d( string fileName,string gridxName, sim::Array<float>& gridx);
int read_profile3d( string fileName,string dataName, sim::Array<int>& data);

int read_profiles( std::string fileName, int &n_x, int &n_z,std::string gridxName, sim::Array<float>& gridx,std::string gridzName,
                            sim::Array<float>& gridz, std::string dataName, sim::Array<float>& data);
void OUTPUT(char outname[],int nX, int nY, float **array2d);
void OUTPUT2d(std::string folder,std::string outname,int nX, int nY, float *array2d);
void OUTPUT1d(std::string folder,std::string outname,int nX, float *array2d);
#endif


