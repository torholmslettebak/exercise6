/*
	C-program to solve the two-dimensional Poisson equation on a unit square,
	using one dimensional eigenvalue decompositions and fast sine transforms

	Based on einar m. ronquist's version, and "modded" with Arne Morten Kvarvigs code
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common.h"
#include "poissoncommon.h"

double alpha = 0.0;

void printArray(int *arr, int l);

double source(double x, double y)
{
	return 5*pow(M_PI,2)*sin(M_PI*x)*sin(2*M_PI*y);
  // return -30.0*pow(y,4)*x*(pow(x,5.0)-1)-30.0*pow(x,4)*y*(pow(y,5)-1);
}

double exact(double x, double y)
{
  return sin(M_PI*x)*sin(2*M_PI*y);
}

void DiagonalizationPoisson2Dfst(Matrix b, const Vector lambda)
{
	Matrix ut = createMatrix(b->rows, b->cols);
	Vector buf;
	int N=b->rows+1;
	int NN=4*N;
	printf("here 1\n");
	#pragma omp parallel for schedule(static) private(buf)
	for (int i=0;i<b->cols;++i)
	{
		buf = createVector(4*(b->rows+1));
		fst(b->data[i], &N, buf->data, &NN); 
		freeVector(buf);
	}
	transposeMatrix(ut, b);
	printf("here 2\n");

	#pragma omp parallel for schedule(static) private(buf)
	for (int i=0;i<ut->cols;++i)
	{
		buf = createVector(4*(b->rows+1));
		fstinv(ut->data[i], &N, buf->data, &NN);
		freeVector(buf);
	}
	// #pragma omp parallel for schedule(static) private(buf)
	printf("here 3\n");

  	for (int j=0;j<b->cols;++j)
    {
		for (int i=0;i<b->rows;++i)
		{
			ut->data[j][i] /= lambda->data[i]+lambda->data[j];
		}
	}
	printf("here 4\n");

	#pragma omp parallel for schedule(static) private(buf)
	for (int i=0;i<b->cols;++i)
  	{
  		buf = createVector(4*(b->rows+1));
    	fst(ut->data[i], &N, buf->data, &NN);
    	freeVector(buf);
  	}
  	transposeMatrix(b, ut);
	printf("here 5\n");

  	#pragma omp parallel for schedule(static) private(buf)
  	for (int i=0;i<ut->cols;++i)
  	{
  		buf = createVector(4*(b->rows+1));
  		fstinv(b->data[i], &N, buf->data, &NN);
  		freeVector(buf);
  	}
  	printf("here 6\n");

  	freeMatrix(ut);
  	printf("here 7\n");
  	// freeVector(buf);
}

void DiagonalizationPoisson2DfstParallell(Matrix b, const Vector lambda)
{
	Matrix ut = createMatrix(b->rows, b->cols);
	Vector buf = createVector(4*(b->rows+1));
	int N=b->rows+1;
	int NN=4*N;
	// #pragma omp parallell for schedule(static)
	for (int i=0;i<b->cols;++i)
	{
		fst(b->data[i], &N, buf->data, &NN); 
	}
	transposeMatrix(ut, b);
	for (int i=0;i<ut->cols;++i)
	{
		fstinv(ut->data[i], &N, buf->data, &NN);
	}
  	for (int j=0;j<b->cols;++j)
    {
		for (int i=0;i<b->rows;++i)
		{
			ut->data[j][i] /= lambda->data[i]+lambda->data[j];
		}
	}
	for (int i=0;i<b->cols;++i)
  	{
    	fst(ut->data[i], &N, buf->data, &NN);
  	}
  	transposeMatrix(b, ut);
  	for (int i=0;i<ut->cols;++i)
  	{
  		fstinv(b->data[i], &N, buf->data, &NN);
  	}
  freeMatrix(ut);
  freeVector(buf);
}
// 
void printArray(int *array, int m)
{
	for (int i = 0; i < m; i++)
	{
		printf(" %d ", array[i]);
	}
	printf("\n");
}

void printResults(int omp, int mpi, int N, double runtime, double error)
{
	printf("omp %d\nmpi %d\nn %d\nruntime %f\nerror %e\n", omp, mpi, N, runtime, error);
}

int main(int argc, char **argv)
{
	Vector grid, lambda=NULL;
	Matrix b, b_t, e;
	int local, N, n, m;
	double h, t1, t2;

	if(argc < 2 )
	{
		printf("need a problem size (in each direction)\n");
		return 1;
	}

	local = 0;
	N = atoi(argv[1]);
	m  = N-1;
	// MPI Initalization and declaration
	int rank, size, tag;
	MPI_Status status;
	init_app(argc, argv, &rank, &size);

	printf("The size = %d, The rank = %d\n", size, rank);
	int *len, *displ;
	splitVector(m, size-1, &len, &displ);
	// printf("Test of len and displ array: \n");
	// printArray(len, m);
	// printArray(displ, m);
	// printf("\n");
	if (size==1)
	{
		// This means program has been run without additional processors
		// it will execute serial within this block
		b = createMatrix(m, m);
		e = createMatrix(m,m);
		grid = equidistantMesh(0.0, 1.0, N); // Create grid

		evalMeshInternal2(b, grid, source, local);
		h = grid->data[1]-grid->data[0];	// find step-value h from grid
		scaleVector(b->as_vec, pow(h,2));
		evalMeshInternal2(e, grid, exact, local);
		axpy(b->as_vec,e->as_vec, alpha);
		lambda = generateEigenValuesP1D(m);	// generate eigenvalues

		t1 = WallTime();
		DiagonalizationPoisson2Dfst(b, lambda);	// Run solution routine
		t2 = WallTime();
		axpy(b->as_vec,e->as_vec,-1.0);  // Calculating errors
		double error = maxNorm(b->as_vec);

		printResults(getMaxThreads(), size, N, t2-t1, error);
	}

	else
	{
		// The program is to be executed with several processors
		b = createMatrixMPI(m, -1, m, m, &WorldComm);
		// printMatrix(b);
		// b_t = createMatrixMPI(m, -1, m, m, &WorldComm);
		// e = createMatrixMPI(m, -1, m, m, &WorldComm);
		if (rank == 0)
		{
			printf("I am %d\n", rank);
			printMatrix(b);
			// printf("Here I am! rank: %d\n", rank);
			grid = equidistantMesh(0.0, 1.0, N); // Create grid
			printf(" Rock me like a hurricane! \n");
			h = grid->data[1]-grid->data[0];	// find step-value h from grid
			lambda = generateEigenValuesP1D(m);	// generate eigenvalues
			// printf("something\n");

			evalMeshInternal2(b, grid, source, local);
			printf("something dark\n");
			scaleVector(b->as_vec, pow(h,2));
			printMatrix(b);
		}
		// time = WallTime();
		// DiagonalizationPoisson2Dfst(b, lambda);	// Run solution routine
		// printf("The b matrix: \n");
		// printMatrix(b);
		// printf("elapsed: %f\n", WallTime()-time);
	}



	close_app();
}