#include <stdlib.h>
#include <stdio.h>
//#include <netcdf.h>
#include <netcdf>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
using namespace std;
using namespace netCDF;

    int main()
{
    std::cout << "hello" << std::endl;
    std::string fileName = "positions1e7.nc";

    ifstream file(fileName.c_str());
        if(!file.good()) {
                    cout<<"ERROR: Cannot find input file ... "<<fileName<<endl;
                            exit(1);
              }

            NcFile nc(fileName.c_str(), NcFile::read);

                NcDim nc_nR(nc.getDim("nP"));
                    //NcDim nc_nZ(nc.getDim("nZ"));

                        int nP = nc_nR.getSize();
                            //int nZ = nc_nZ.getSize();

                                NcVar nc_x(nc.getVar("x"));
                                    vector<float> x;
                                        x.resize(nP);
                                            nc_x.getVar(&x[0]);
                                   // NcVar nc_z(nc.getVar("z"));
                                NcVar nc_y(nc.getVar("y"));
                                    vector<float> y;
                                        y.resize(nP);
                                            nc_y.getVar(&y[0]);
                                NcVar nc_z(nc.getVar("z"));
                                    vector<float> z;
                                        z.resize(nP);
                                            nc_z.getVar(&z[0]);
                                NcVar nc_h(nc.getVar("hitWall"));
                                    vector<float> h;
                                        h.resize(nP);
                                            nc_h.getVar(&h[0]);
                     std::cout << nP << std::endl;
                     std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;
  vector<float> r(nP);
  for(int i=0;i<nP;i++)
  {
    r[i] = sqrt(x[i]*x[i] + y[i]*y[i]);
  }
  vector<int> beads(13,0);
  int ii = 0;
  for(int i=0;i<nP;i++)
  {
      if(r[i] < 0.07 && z[i] > 0.01275 && h[i] == 1.0 && z[i] < 0.19)
      {
           ii = floor((z[i]-0.01275)/0.01);
           beads[ii] = beads[ii]+1;
      }
  }
  for(int i=0; i<13;i++)
      std::cout << "Bead " << i << " " << beads[i] << std::endl;
}
