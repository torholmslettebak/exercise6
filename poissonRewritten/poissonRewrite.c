/*
	C-program to solve the two-dimensional Poisson equation on 
	a unit square using one-dimensional eigenvalue decompositions
	and fast sine transforms

	einar m. ronquist
	ntnu, october 2000
	revised, october 2001
*/

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include "common.h"
#include "poissoncommon.h"

typedef double Real;

Real source(Real x, Real y)
{
	return (5*M_PI*M_PI*sin(M_PI*x)*sin(2*M_PI*y));
}

Real exact(Real x, Real y)
{
	return (sin(x*M_PI)*sin(2*y*M_PI));
}

// void printVector2(Real *vec, int m)
// {
// 	int i = 0;
// 	for (i = 0; i < m; i++)
// 	{
// 		printf("	%lf	", vec[i]);
// 	}
// 	printf("\n");
// }
// void printVector3(int *vec, int m)
// {
// 	int i = 0;
// 	for (i = 0; i < m; i++)
// 	{
// 		printf("	%d	", vec[i]);
// 	}
// 	printf("\n");
// }

// void printMatrix2(Real **matrix, int m, int n)
// {
//   int i,j;
// 	for(i = 0; i < m; i++)
// 	{
// 		for(j = 0; j < n; j++)
// 		{
// 			printf("	%lf	", matrix[i][j]);
// 		}
// 		printf("\n");
// 	}
// }

void assembleMatrix()
{

}

void matrixToVec(Real *toVec, Real **matrix, int rows, int cols, int rank, int *displ)
{
	int startpoint = 0;
	for (int i = 0; i < rows; i++)
	{
		appendArray(toVec, matrix[i], startpoint, cols);
		startpoint = startpoint+cols;
		// printf("inside matrixToVec\n");
		// printVector2(toVec, rows*cols);
	}
}


Real* myEquidistantMesh(Real x0, Real x1, int N)
{
  double h = (x1-x0)/N;
  Real *result = calloc(N+1, sizeof(double));

  int i;

  for (i=0;i<N+1;++i)
    result[i] = x0+i*h;

  return result;
}

void DiagonalizationPoisson2Dfst(int n, int rank, int size)
{
	double t1 = WallTime();
	int m, nn;
	Real **b, **b_t, **z, *diag, h, umax;
	int *len, *displ;
	
	m = n-1;
	nn = 4*n;
	splitVector(m, size, &len, &displ);
	printVector3(len, size);
	printVector3(displ, size);

	diag = createRealArray (m);
	b    = createReal2DArray (len[rank], m);
	z    = createReal2DArray (omp_get_max_threads(), nn);

	h    = 1./(Real)n;
	// printf("stepsize h = %lf\n", h);

	Real* grid = myEquidistantMesh(0.0, 1.0, n);
	myEvalMeshInternal(b, grid, source, n, len[rank], rank, displ);
	// printf("finished with my myEvalMeshInternal\n");
	myScaleMatrix(b, len[rank], h*h, m);
	// printf("Finished with myScaleMatrix \n");
	// printVector2(b[0], m);
	// printf("after printmatrix\n");
	  // Generates the eigenvalues and stores it in diagonal
	for (int i=0; i < m; i++) 
	{
		diag[i] = 2.*(1.-cos((i+1)*M_PI/(Real)n));
	}

	// printf("After eigenvalues\n");
#pragma omp parallel for schedule(static)
	for (int j=0; j < len[rank]; j++) 
	{
		fst_(b[j], &n, z[omp_get_thread_num()], &nn);
	}
	// printf("After fst\n");
	// printf("Before first transpose\n");

	// printMatrix2(b, len[rank], m);
	myTranspose(b, rank, size, m, len, displ);
	// printMatrix2(b, len[rank], m);
	// printf("After first transpose\n");
#pragma omp parallel for schedule(static)
	for (int i=0; i < len[rank]; i++) 
	{
		fstinv_(b[i], &n, z[omp_get_thread_num()], &nn);
	}
	// printf("after fstinv, i am rank: %d\n", rank);

#pragma omp parallel for schedule(static)
	for (int j=0; j < len[rank]; j++) 
	{
		for (int i=0; i < m; i++) 
		{
			b[j][i] = b[j][i]/(diag[i]+diag[j+displ[rank]]);
		}
	}
	// printf("after insertion of eigenvalues\n");
#pragma omp parallel for schedule(static)	  
	for (int i=0; i < len[rank]; i++) 
	{
		fst_(b[i], &n, z[omp_get_thread_num()], &nn);
	}
	// printf("after fst number 2\n");
	myTranspose(b, rank, size, m, len, displ);

	// printf("after another transpose, and before fstinv number 2\n");
#pragma omp parallel for schedule(static)
	for (int j=0; j < len[rank]; j++) 
	{
		fstinv_(b[j], &n, z[omp_get_thread_num()], &nn);
	}
	// printf("after fstinv number 2\n");
	// printMatrix2(b, len[rank], m);
	Real *b_as_vec = createRealArray(m*len[rank]);
	matrixToVec(b_as_vec, b, len[rank], m, rank, displ);
	printVector2(b_as_vec, m*len[rank]);
	// every rank now has its own vector representation its matrix
	// Now to gather it all at root
	Real *recvbuf = createRealArray(m*m);
	int *recvcounts, *recvdispl;
	splitVector(m*m, size, &recvcounts, &recvdispl);
	// printVector3(recvcounts, size);
	// MPI_Gatherv(b_as_vec, m*len[rank], MPI_DOUBLE, recvbuf, recvcounts, recvdispl, MPI_DOUBLE, 0, WorldComm);
	if (rank = 0)
	{
		// printf("I am rank %d\n", rank);
		// printVector2(recvbuf, m*m);
	}
	// if (rank == 0)
	// {
	// 	umax = 0.0;
	// 	for (int j=0; j < m; j++) 
	// 	{
	// 		for (int i=0; i < m; i++) 
	// 		{
	// 			if (b[j][i] > umax) umax = b[j][i];
	// 		}
	// 	}
	// 	printf (" umax = %e \n",umax);

	// }	
}


int main(int argc, char *argv[]) 
{
	int rank, size;

	Real *diag, **b, **bt, *z;
	Real pi, h, umax;
	int i, j, n, m, nn;
	// start MPI
	MPI_Status status;
	init_app(argc, argv, &rank, &size);

	  /* the total number of grid points in each spatial direction is (n+1) */
	  /* the total number of degrees-of-freedom in each spatial direction is (n-1) */
	  /* this version requires n to be a power of 2 */

	if( argc < 2 )
	{
		printf("need a problem size\n");
		close_app();
		return 1;
	}
	else
	{
		n  = atoi(argv[1]);
		DiagonalizationPoisson2Dfst(n, rank, size);
		close_app();
	}
} 