#ifndef _IO_
#define _IO_

using namespace std;

int read_ar2Input( std::string fileName, double *Bfield[]);

int read_profileNs( std::string fileName,std::string nzName,std::string nxName,int &n_x,int &n_z );
int read_profile2d( string fileName,string dataName, std::vector<double>& data);
int read_profile1d( string fileName,string gridxName, std::vector<double>& gridx);

int read_profiles( std::string fileName, int &n_x, int &n_z,std::string gridxName, std::vector<double>& gridx,std::string gridzName,
                            std::vector<double>& gridz, std::string dataName, std::vector<double>& data);
void OUTPUT(char outname[],int nX, int nY, double **array2d);
void OUTPUT2d(std::string folder,std::string outname,int nX, int nY, double *array2d);
void OUTPUT1d(std::string folder,std::string outname,int nX, double *array2d);
#endif


