

#include <iostream>
#include <stdlib.h>
#include "h1.h"
#include <mpi.h>
using namespace std;

int main( int argc, char *argv[] )
{
int nProc = 16;
int myID = 1;
int mster;
int nWRs;
double tt0;
double tt1;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myID);
    mster=0;
    nWRs = nProc-1;
    
    //std::cout << "MAIN running on " << nWRs << " WRs" << std::endl;
    
    if (myID == mster)
    	{
    	tt0 = MPI_Wtime();
    	MASTER(nWRs,mster);
    	tt1 = MPI_Wtime();
    	std::cout << "MAIN MR timing = " << tt1-tt0 << "sec on " <<nWRs << "WRs" << std::endl;	
    	}
    else
    	{
    	WORKER( nWRs, myID );
    
    	}
	MPI_Finalize();
  return 0;
}