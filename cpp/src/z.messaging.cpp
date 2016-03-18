#include "h1.h"
#include <cmath>
#include <mpi.h>




void SEND_2doutput_MPI(int myID,int nX, int nY,double **arr)
{
	   int buffer[3];
   MPI_Status status;

   /* Send first point, number of points and results to master */
   buffer[0] = myID;
   buffer[1] = nX;
   buffer[2] = nY;
   
   
   MPI_Send(&buffer, 3, MPI_INT, 0, 30, MPI_COMM_WORLD);

   MPI_Send(arr[0], nX*nY, MPI_DOUBLE, 0, 40, MPI_COMM_WORLD);
   

}
void RECV_output_MPI(int nWRs, int M, double global[])
{
int start;
int npts;
int buffer[2];
double *arr1;
for (int i = 1; i <= nWRs; i++) {
		MPI_Recv(&buffer, 2, MPI_INT, i, 30, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		start = buffer[0];
		npts = buffer[1];

		arr1 = new double[npts];
		MPI_Recv(&arr1[0], M, MPI_DOUBLE, i, 40, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		
		 	for (int j = 0; j < M; j++)
 				{	
   			     global[j+start] = arr1[j];
   			     
    		    }
		}
		

}
void RECV_2doutput_MPI(int nWRs, int nX, int nY,double **local, double **global)
{
int ID;
int npts;
int nX_rec;
int nY_rec;
int buffer[3];



for (int i = 1; i <= nWRs; i++) {
		MPI_Recv(&buffer, 3, MPI_INT, i, 30, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		ID = buffer[0];
		nX_rec = buffer[1];
		nY_rec = buffer[2];
		//std::cout << start << Mz_rec << Mr_rec << std::endl;


		MPI_Recv(local[0], nX*nY, MPI_DOUBLE, i, 40, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		
		//std::cout << "Start" << start << std::endl;
			 for(int j=0 ; j<nX ; j++)
				{	
					for(int k=0 ; k<nY ; k++)
					{
					global[j][k] = global[j][k] + local[j][k];
					}
				}
		}
		
		
}

