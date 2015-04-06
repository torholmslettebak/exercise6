#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common.h"
#include "poissoncommon.h"


double source(double x, double y)
{
	return x*y;
  // return -30.0*pow(y,4)*x*(pow(x,5.0)-1)-30.0*pow(x,4)*y*(pow(y,5)-1);
}

double exact(double x, double y)
{
  return sin(M_PI*x)*sin(2*M_PI*y);
}


// Generally faster for big matrices >256*256
void transposeOmp(Matrix A_t, Matrix A)
{
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < A->rows; i++)
	{
		for (int j = 0; j < A->rows; j++)
		{
			A_t ->data[j][i] = A -> data[i][j];
		}
	}
}

void transposeMPI(Matrix A_t, Matrix A, int rank, int size)
{
	for(int i = rank; i < (A->rows); i + size)
	{
		// printf("here\n");
		for (int j = 0; j < (A->rows); j++)
		{
			A_t->data[j][i] = A->data[i][j];
		}
	}
}

void fillMatrix(Matrix A, int rank)
{
	double counter = 1;
	for(int i = 0; i < A->rows; i++)
	{
		for(int j = 0; j < A->rows; j++)
		{
			A->data[i][j] = counter; 
			// printf("I am rank %d, and my A[i][j]: %lf\n", rank, A->data[i][j]);
			counter++;
		}
	}
}

void printArrayOfDouble(double *array, int m, int rank)
{
	printf("I AM RANK: %d\n", rank);
	for (int i = 0; i < m; i++)
	{
		printf(" %lf ", array[i]);
	}
	printf("\n");
}

void printArrayOfInts(int *array, int m)
{
	for (int i = 0; i < m; i++)
	{
		printf(" %d ", array[i]);
	}
	printf("\n");
}
// void insertSubMatrix(Matrix subMat, Matrix toMatrix, int level, int startCol)
// {
// 	for (int i = 0; i < )
// 	{

// 	}
// }
// appendToBuffVec(Matrix submat, Vector vec, rank, blockSize, N)
// {
// 	Vector subAsVec = submat->as_vec;

// }

void assemblePartMatrix(Matrix A_t, Vector recvbuf, int rank, int blockSize, int size)
{
	int row = 0;
	int block = 0;
	int numberOfBlocks = size;
	// for (int i = rank; i < rank+blockSize; i++) // loops over the columns where the vector is to be placed.
	// {
	// 	for (int j = 0;)
	// }
	
	
	for (int i = 0; i < numberOfBlocks; i ++)
	{
		int startIndex = i * blockSize*blockSize; 
		for(int j = startIndex; j < startIndex + blockSize*blockSize; j++)
		{
			//  int j er nå index i recvbuf
			// blocksize*blocksize antall elementer skal nå inn i A_t
			// recvbuf->data[j];
			A_t -> data[rank][i] = recvbuf->data[j];
			A_t -> data[rank][i+1] = recvbuf->data[]
		}
	}
	
}





int main(int argc, char **argv)
{	
	Matrix A, A_t;
	int rank, size, tag;
	double t1, t2;
	int N = atoi(argv[1]);
	MPI_Status status;
	int mpi_top_coords[2];
	int mpi_top_sizes[2];
	tag = 1;
	Matrix recvMat;
	
	// Vector recvbuf, recvcounts, rdispls;
	int *recvcounts, *rdispls;
	// Vector recvbuf = createVector(N);
	// recvbuf = calloc(N,sizeof(double));

	init_app(argc, argv, &rank, &size);

	int *len, *displ;

	A = createMatrix(N, N);
	A_t = createMatrix(N,N);
	fillMatrix(A, rank);
	// printf("The displ vector: \n");
	// printArrayOfInts(displ, size);
	// Now to split matrix into blocks
	// Assuming even number of processes
	int blockSize = N/size;
	// each processor is responsible for a column of blocks
	// this gives balanced workload for the processors
	// equal number of sends and receives, as well as number of transposes
	// number if blocks for each process = number of processes
	// if (rank == 0)
	// {
	Vector bufferVec = createVector(N*blockSize);
	Vector recvbuf = createVector(N*blockSize);

	// if (rank==0)
	// {
	// 	printf("bufferVec size: %d\n", bufferVec->);
	// }
	for (int i = 0; i < size; i++) // sizeequals number of blocks per proc
	{

		Matrix submat = subMatrix(A, i*blockSize, blockSize, rank*blockSize, blockSize);
		Matrix submat_t = createMatrix(blockSize, blockSize);
			// printMatrix(submat);
		transposeMatrix(submat_t, submat);
		// if (rank == 2)
		// {
		// 	printMatrix(submat_t);
		// }
			// insertSubMatrix(submat_t, A_t, i, rank);
			// buffermatrix->data[i]
		// KK, store each processes block column in a vector, and use gather or alltoall to merge together
		appendVector(bufferVec, submat_t->as_vec, i*blockSize*blockSize, blockSize);
	}
	if(rank==0)
	{
		printVector(A->as_vec);
		// printVector(bufferVec);
	} 

	// if(rank==1)
	// {
	// 	printVector(bufferVec);
	// } 
	// if(rank==2)
	// {
	// 	printVector(bufferVec);
	// } 
	// if(rank==3)
	// {
	// 	printVector(bufferVec);
	// } 
	// Now each process has its own buffer vector. perform alltoall to send the parts that are to be swapped to each other
	splitVector(N*blockSize, size, &len, &displ);
	// printArrayOfInts(len, size);
	// printArrayOfInts(displ, size);
	MPI_Alltoallv(bufferVec->data, len, displ, MPI_DOUBLE, recvbuf->data,len, displ, MPI_DOUBLE, WorldComm);

	if(rank==0)
	{
		printVector(recvbuf);
	} 

	if(rank==1)
	{
		printVector(recvbuf);
	} 
	if(rank==2)
	{
		printVector(recvbuf);
	} 
	if(rank==3)
	{
		printVector(recvbuf);
	} 

	// Now every column is transposed, and ready to be sent to assembly
	// MPI_Gather(bufferVec->data, bufferVec->len, MPI_DOUBLE, recvMat->data[rank], bufferVec->len, MPI_DOUBLE, 0, WorldComm);
	// printf("recvbuf->data: \n");
	// printVector(recvbuf);

	// Make each process assemble its own A_t with its own part of recvbuff
	// Then send this with MPI_Gather(v) too rank 0 
	//  if this works transpose is done and return A_t
	assemblePartMatrix(A_t, recvbuf, rank, blockSize);
	MPI_Gather(recvbuf->data, N*blockSize, MPI_DOUBLE, A_t->as_vec->data, N*blockSize, MPI_DOUBLE, 0, WorldComm);	
	if(rank==0)
	{
		printMatrix(A_t);
	}
	// recvMat = createMatrix(size, N*blockSize);
	// printMatrix(recvMat);
	// MPI_Barrier(WorldComm);
	// if (rank == 0)
	// {
	// 	printVector(recvMat->as_vec);
	// 	// printArrayOfDouble(recvMat->data[1], blockSize*N, rank);
	// 	printf("THe number of rows %d, the number of cols %d\n", recvMat->rows, recvMat->cols);
	// }

	// }
}

