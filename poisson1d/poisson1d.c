#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common.h"
// #include "solvers.h"
#include "poissoncommon.h"

double alpha=0.1;

double exact(double x)
{
  return x*(pow(x,5.0)-1.0);
}

double source(double x)
{
  return -30*pow(x,4.0);
}

void DiagonalizationPoisson1Dfst(Vector u, const Vector lambda)
{
	Vector btilde = createVector(u->len);
	Vector buf = createVector(4*(u->len+1));
	int i;
	int N=u->len+1;
	int NN=4*N;
	copyVector(btilde, u);
	fst(btilde->data, &N, buf->data, &NN);
	printf("after fst btidle vector: \n");
	printVector(btilde);
	printf("\n");
	for (i=0;i<btilde->len;++i)
		btilde->data[i] /= (lambda->data[i]+alpha);

	printf("after the for loop, btilde: \n");
	printVector(btilde);
	printf("\n");
	fstinv(btilde->data, &N, buf->data, &NN);
	
	printf("The btidle after fstinv\n");
	printVector(btilde);
	printf("\n");
	printf("The bud vector: \n");
	printVector(buf);
	printf("\n");
	copyVector(u, btilde);
	freeVector(btilde);
	freeVector(buf);
}

int main(int argc, char** argv)
{
	// will try to do this with fst based approach
	int N, flag;
	double time, h, tol=1e-8;
	Vector grid, b, e, lambda=NULL;
	Matrix Q;
	if (argc < 3) 
	{
		printf("need two parameters, N and flag [and tolerance]\n");
		printf(" - N is the problem size (in each direction\n");
		// printf(" - flag = 1  -> Dense LU\n");
		// printf(" - flag = 2  -> Dense Cholesky\n");
		// printf(" - flag = 3  -> Full Gauss-Jacobi iterations\n");
		// printf(" - flag = 4  -> Full Gauss-Jacobi iterations using BLAS\n");
		// printf(" - flag = 5  -> Full Gauss-Seidel iterations\n");
		// printf(" - flag = 6  -> Full Gauss-Seidel iterations using BLAS\n");
		// printf(" - flag = 7  -> Full CG iterations\n");
		// printf(" - flag = 8  -> Matrix-less Gauss-Jacobi iterations\n");
		// printf(" - flag = 9  -> Matrix-less Gauss-Seidel iterations\n");
		// printf(" - flag = 10 -> Matrix-less Red-Black Gauss-Seidel iterations\n");
		// printf(" - flag = 11 -> Diagonalization\n");
		printf(" - flag = 1 -> Diagonalization - FST\n");
		// printf(" - flag = 13 -> Matrix-less CG iterations\n");
		return 1;
	}
	N = atoi(argv[1]);
	flag = atoi(argv[2]);
	if (argc > 3)
		tol = atof(argv[3]);
	if (N < 0) 
	{
		printf("invalid problem size given\n");
		return 2;
	}

	if ( (N !=1) && (N & (N-1)) )
	{
		printf("problem size needs to be a power of 2: \n");
		return 3;
	}

	h = 1.0 / N;	// The mesh spacing
	// creates a equidistance grid, a vector for 1d case
	// we will use standard x0 = 0 and X1 = 1
	grid = equidistantMesh(0.0, 1.0, N);

	b = createVector(N-1);
	e = createVector(N-1);
	evalMeshInternal(b, grid, source);
	printf("the grid\n");
	printVector(grid);
	printf("the b vector\n");
	printVector(b);
	printf("\n");
	evalMeshInternal(e, grid, exact);
	scaleVector(b, pow(h, 2));
	printf("The b vector\n"); //this gives f*h^2, and here f = 1
	printVector(b);
	printf("The e vector\n");
	printVector(e);
	if (flag >= 0 && flag < 2)
		lambda = generateEigenValuesP1D(N-1);

	// printf("Print the eigenmatrix\n");
	// printMatrix(Q);
	// printf("The eigenvalues\n");
	// printVector(lambda);
	time = WallTime();
	printf("the lambda vector\n");
	printVector(lambda);
	printf("\n");
	if (flag == 1)
	{
		DiagonalizationPoisson1Dfst(b,lambda);
	}
	printf("the b vector\n");
	printVector(b);
	printf("\n");

	axpy(b, e, alpha);
	printf("elapsed: %f\n", WallTime()-time);

	evalMeshInternal(e, grid, exact);
	axpy(b,e,-1.0);
	printf("max error: %e\n", maxNorm(b));
	freeVector(grid);
	freeVector(b);
	freeVector(e);
	if (lambda)
		freeVector(lambda);
	return 0;
}
