#ifndef _IO_
#define _IO_

using namespace std;

int read_ar2Input( std::string fileName);

int read_profileNs( std::string fileName,int &n_x,int &n_z );

int read_profiles( std::string fileName, int &n_x, int &n_z,std::vector<double>& gridx,
                            std::vector<double>& gridz, std::vector<double>& data);
#endif


