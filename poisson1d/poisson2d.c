#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "common.h"
#include "poissoncommon.h"


double source(double x, double y)
{
	return -30.0*pow(y,4)*x*(pow(x,5.0)-1)-30.0*pow(x,4)*y*(pow(y,5)-1);
}

double exact(double x, double y)
{
	return x*(pow(x,5)-1.0)*y*(pow(y,5)-1.0);
}

void DiagonalizationPoisson2Dfst(Matrix b, const Vector lambda)
{
	Vector buf = createVector(4*(b->rows+1));
	int i,j;
	Matrix ut = createMatrix(b->rows, b->cols);
	int N=b->rows+1;
	int NN=4*N;

	for (i=0;i<b->cols;++i)
	{
		fst(b->data[i], &N, buf->data, &NN);
	}
	
	transposeMatrix(ut, b);
	
	for (i=0;i<ut->cols;++i)
	{
		fstinv(ut->data[i], &N, buf->data, &NN);
	}

	for (j=0;j<b->cols;++j)
	{
		for (i=0;i<b->rows;++i)
		{
			ut->data[j][i] /= lambda->data[i]+lambda->data[j];
		}
	}
	for (i=0;i<b->cols;++i)
	{
  		fst(ut->data[i], &N, buf->data, &NN);
	}
  	
  	transposeMatrix(b, ut);
  	
  	for (i=0;i<ut->cols;++i)
  	{
		fstinv(b->data[i], &N, buf->data, &NN);
  	}

  	freeMatrix(ut);
  	freeVector(buf);
}

int main(int argc, char** argv)
{
	int i, j, N, flag, local;
	Matrix A=NULL, Q=NULL;
	Matrix b, e;
	Vector grid, lambda=NULL;
	double time, sum, h;
	if (argc < 3) 
	{
		printf("need two parameters, N and flag\n");
		printf(" - N is the problem size (in each direction\n");
		printf(" - flag = 12 -> Diagonalization, fst based\n");
		return 1;
	}
	N=atoi(argv[1]);
	flag=atoi(argv[2]);
	if (N < 0) 
	{
		printf("invalid problem size given\n");
		return 2;
	}
	if (flag < 0 || flag > 13) 
	{
		printf("invalid flag given\n");
		return 3;
	}
	grid = equidistantMesh(0.0, 1.0, N);
	b = createMatrix(N-1, N-1);
	e = createMatrix(N-1, N-1);
	local = 0;	// If local = 1, space is reserved for boundary in the matrix
	// evalMeshInternal fills the b vector with values from the source-function
	evalMeshInternal2(b, grid, source, local);
	h = grid->data[1]-grid->data[0];
	scaleVector(b->as_vec, pow(h, 2));

	if (flag >= 11 && flag < 13)
		lambda = generateEigenValuesP1D(N-1);

	time = WallTime();

	else if (flag == 12)
	{
		DiagonalizationPoisson2Dfst(b, lambda);
	}

	printf("elapsed: %f\n", WallTime()-time);

	evalMeshInternal2(e, grid, exact, local);
	axpy(b->as_vec,e->as_vec,-1.0);

	printf("max error: %e\n", maxNorm(b->as_vec));


	freeMatrix(b);
	freeMatrix(e);
	freeVector(grid);
	freeVector(lambda);
	
	return 0;
}