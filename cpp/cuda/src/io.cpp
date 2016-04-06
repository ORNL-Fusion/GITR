#include "h1.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <cstring>
using namespace std;



//IO


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
