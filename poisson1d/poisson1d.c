#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common.h"
// #include "solvers.h"
#include "poissoncommon.h"


int main(int argc, char** argv)
{
	// will try to do this with fst based approach
	int N, flag;
	double time, h, tol=1e-4;
	Vector grid;
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
	printVector(grid);
}
