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

void fillMatrix(Matrix A)
{
	double counter = 1;
	for(int i = 0; i < A->rows; i++)
	{
		for(int j = 0; j < A->rows; j++)
		{
			A->data[i][j] = counter; 
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

int main(int argc, char **argv)
{	
	Matrix A, A_t;
	int rank, size, tag;
	double t1, t2;
	int N = atoi(argv[1]);
	MPI_Status status;
	tag = 1;
	
	// Vector recvbuf, recvcounts, rdispls;
	int *recvcounts, *rdispls;
	double *recvbuf;
	recvbuf = calloc(N,sizeof(double));

	init_app(argc, argv, &rank, &size);

	int *len, *displ;
	splitVector(N, size, &len, &displ);

	if (size == 1)
	{
	// printMatrix(A);
		A = createMatrix(N,N);
		A_t = createMatrix(N,N);
		fillMatrix(A);

		t1 = WallTime();
		transposeMatrix(A_t, A);
		t2 = WallTime();
		printf("The common.c transpose performed in: %lf\n", t2-t1);

		t1 = WallTime();
		transposeOmp(A_t, A);
		t2 = WallTime();
		printf("The transposeOmp performs with %d threads, at %lf \n", getMaxThreads(), t2-t1);
	// printMatrix(A_t);
	}

	else
	{
		Vector matrixAsVec;
		A = createMatrix(N,N);
		A_t = createMatrix(N,N);
		fillMatrix(A);
		matrixAsVec = A->as_vec;
		// printArrayOfDouble(A->data[rank], N, rank);
			// transposeMPI(A_t, A, rank, size, l,en, displ);

		// if (rank == 0)
		// {
		// 	printf("I am rank %d, and len = %d, displ = %d \n", rank, len[rank], displ[rank]);
		// }

		// if (rank !=0)
		// {
		// 	printf("I am rank %d, and len = %d, displ = %d \n", rank, len[rank], displ[rank]);
		// }
		// printf("The send buffer: \n");
		// printVector(matrixAsVec);
		// printf("Before\n");
		// printf("Rank before alltoallv: %d\n", rank);
		MPI_Alltoallv((A->data[rank]), len, displ, MPI_DOUBLE, recvbuf, len, displ, MPI_DOUBLE, WorldComm);
		// printf("After\n");
		if (rank!=0)
		{
			MPI_Send(recvbuf, N, MPI_DOUBLE, 0, tag, WorldComm);
		}
		// MPI_Barrier(WorldComm);
		if (rank==0)
		{
			A_t->data[rank] = recvbuf;
			Vector rec;
			for (int i = 1; i<size;i++)
			{
				rec = createVector(N);
				MPI_Recv(rec->data, N, MPI_DOUBLE, i, tag, WorldComm, MPI_STATUS_IGNORE);
				A_t->data[i]=rec->data;
				printVector(rec);
			}
			printMatrix(A_t);
		}

	}

}