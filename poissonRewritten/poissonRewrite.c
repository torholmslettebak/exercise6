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


Real **assembleMatrix(Real *b_as_vec, Real *finalMatrix_as_vec, int *len, int *displ, int tag, int m, int size)
{
	// b_as_vec is roots own partMatrix
	if(len[0]>0)
	{
		appendArray(finalMatrix_as_vec, b_as_vec, 0, len[0]*m);
	}
	for (int i = 1; i < size; i++)
	{
		int recvCount = len[i]*m;
		int firstIndex = displ[i]*m;
		MPI_Recv(&finalMatrix_as_vec[firstIndex], recvCount, MPI_DOUBLE, i, tag, WorldComm, MPI_STATUS_IGNORE);
	}

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

	double generateError(Real *localVec, Real *exactVec, int m, int rank, int *len, int *displ)
	{
		// now to find local max local error
		Real localError, globalError;
		double alpha = -1.0;
		int stride = 1;
		myaxpy(localVec, exactVec, alpha, m*len[rank], stride);
		localError = myMaxNorm(localVec, m*len[rank], 1, 1, SelfComm);
		// Now to find global error through MPI_Reduce
		MPI_Reduce(&localError, &globalError, 1, MPI_DOUBLE, MPI_MAX, 0, WorldComm);
		return globalError;
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
	double time1 = WallTime();
	int m, nn;
	Real **b, **e, **z, *diag, h;
	int *len, *displ;
	
	m = n-1;
	nn = 4*n;
	if(rank == 0)
		printf("time  = %lf\n", WallTime()-time1);
	splitVector(m, size, &len, &displ);

	diag = createRealArray (m);
	b    = createReal2DArray (len[rank], m);
	z    = createReal2DArray (omp_get_max_threads(), nn);
	h    = 1./(Real)n;
	e 	 = createReal2DArray (len[rank], m);
	Real* grid = myEquidistantMesh(0.0, 1.0, n);
	myEvalMeshInternal(b, grid, source, n, len[rank], rank, displ);
	myScaleMatrix(b, len[rank], h*h, m);

	// create solution matrix for each processor
	myEvalMeshInternal(e, grid, exact, n, len[rank], rank, displ);

	// Generates the eigenvalues and stores it in diagonal
	for (int i=0; i < m; i++) 
	{
		diag[i] = 2.*(1.-cos((i+1)*M_PI/(Real)n));
	}

#pragma omp parallel for schedule(static)
	for (int j=0; j < len[rank]; j++) 
	{
		fst_(b[j], &n, z[omp_get_thread_num()], &nn);
	}

	myTranspose(b, rank, size, m, len, displ);

#pragma omp parallel for schedule(static)
	for (int i=0; i < len[rank]; i++) 
	{
		fstinv_(b[i], &n, z[omp_get_thread_num()], &nn);
	}

#pragma omp parallel for schedule(static)
	for (int j=0; j < len[rank]; j++) 
	{
		for (int i=0; i < m; i++) 
		{
			b[j][i] = b[j][i]/(diag[i]+diag[j+displ[rank]]);
		}
	}

#pragma omp parallel for schedule(static)	  
	for (int i=0; i < len[rank]; i++) 
	{
		fst_(b[i], &n, z[omp_get_thread_num()], &nn);
	}

	myTranspose(b, rank, size, m, len, displ);

#pragma omp parallel for schedule(static)
	for (int j=0; j < len[rank]; j++) 
	{
		fstinv_(b[j], &n, z[omp_get_thread_num()], &nn);
	}


	Real *b_as_vec = createRealArray(m*len[rank]);
	Real *e_as_vec = createRealArray(m*len[rank]);
	matrixToVec(b_as_vec, b, len[rank], m, rank, displ);
	matrixToVec(e_as_vec, e, len[rank], m, rank, displ);

	// every rank now has its own vector representation its matrix
	// Now to gather it all at root
	int tag = 1;
	// One variant is to send each solution vector to gathering at root
	// but since this already is a part solution, it should be possible to calculate max error locally
	// and then calculate max error in total over MPI_Allreduce with MPI_MAX
	// if(rank!=0)
	// {
	// 	MPI_Send(b_as_vec, m*len[rank], MPI_DOUBLE, 0, tag, WorldComm);		
	// }

	// Lets try the other variant to see if there can be a speedup  this way
	// Each process needs it's own part of solution matrix
	
	double error = generateError(b_as_vec, e_as_vec, m, rank, len, displ);


	if (rank == 0)
	{
		// double time2 = WallTime();
		// Real *finalMatrix_as_vec = createRealArray(m*m);
		// assembleMatrix(b_as_vec, finalMatrix_as_vec,len, displ, tag, m, size);
		// // Generate exact solution matrix
		// Real **exactSol = createReal2DArray(m, m);
		// evalMeshInternal2Arrays(exactSol, grid, exact, n);
		// Real *exactSol_as_vec = createRealArray(m*m);
		// matrixToVec(exactSol_as_vec, exactSol, m, m, rank, displ);
		// // Compute error
		// double alpha = -1.0;
		// int stride = 1;
		// int vecLen = m*m;
		// myaxpy(finalMatrix_as_vec, exactSol_as_vec, alpha, vecLen, stride);
		// double error = myMaxNorm(finalMatrix_as_vec, m*m, 1, 1, SelfComm);

		// printf("Time = %lf\n", time2-time1);
		// printf("Error = %e\n", error);
		printf("runtime = %lf\n", WallTime()-time1);
		printf("Error = %e\n", error);
		
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
	int rank, size, n;
	double time1;
	// start MPI
	MPI_Status status;
	init_app(argc, argv, &rank, &size);
	if( argc < 2 )
	{
		printf("need a problem size\n");
		close_app();
		return 1;
	}
	else
	{
		n  = atoi(argv[1]);
		Real *finalMatrix_as_vec = createRealArray((n-1)*(n-1));
		
		DiagonalizationPoisson2Dfst(n, rank, size);
				
		close_app();
	}
} 