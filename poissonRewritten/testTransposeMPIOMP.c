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

void fillMatrix(Matrix A, int N, int rank, int rows)
{
	for(int i = 0; i < N; i++) { 
        /* Initialize data on the rows of the matrix owned by this rank */ 
		if (i >= rank*rows && i < (rank+1)*rows) { 
			for(int j = 0; j < N; j++) { 
				A->data[i][j] = i*N + j; 
			} 
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
	int size, rank, cols; 
    int i, j, k, N; 
    int row, col;
    int *displ, *len;
    double temp;
    N = atoi(argv[1]);
    MPI_Datatype mpi_all_unaligned_t, mpi_all_t; 
    init_app(argc, argv, &rank, &size);
    Matrix transposed = createMatrix(N,N);
    Matrix A = createMatrixMPI(N, N, N, N, &SelfComm);
    Matrix A_t = createMatrixMPI(N, N, N, N, &SelfComm);
    cols = N/size;
    MPI_Type_vector(N, cols, N, MPI_DOUBLE, &mpi_all_unaligned_t);
    MPI_Type_create_resized(mpi_all_unaligned_t, 0, cols*sizeof(double), &mpi_all_t);
    MPI_Type_free(&mpi_all_unaligned_t); 
    MPI_Type_commit(&mpi_all_t); 

    fillMatrix(A, N, rank, cols);
    // MPI_Barrier(WorldComm);
    // printMatrix(A);
    transposeMatrix(A_t, A);
    // printMatrix(A_t);
    splitVector(N*N, size, &len, &displ);

    displ = displ+(2*rank);
    if (rank==3)
    {
    	// printArrayOfInts(len, size);
    	// printArrayOfInts(displ, size);
    	// printVector(A_t->as_vec);
    	// printMatrix(A_t);
    	// printArrayOfDouble(A_t->data[0], N, rank);
    }

    MPI_Allgather(A_t->as_vec->data, N*N, MPI_DOUBLE, transposed->as_vec->data, N*N, MPI_DOUBLE, SelfComm);
    // MPI_Alltoall(A_t->data, N*, displ, MPI_DOUBLE, transposed->data, len, displ, MPI_DOUBLE, SelfComm);
    // MPI_Allreduce(A_t->as_vec->data, transposed->as_vec->data, N*N, MPI_DOUBLE, MPI_SUM, SelfComm);
    // MPI_Alltoall(A_t->as_vec->data, N*N, displ, MPI_DOUBLE, transposed->data, len, displ, MPI_DOUBLE, SelfComm);
    
    if(rank == 0)
    {
    	printMatrix(transposed);
    }
}

