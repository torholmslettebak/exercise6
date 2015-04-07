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

typedef double Real;

/* function prototypes */
Real *createRealArray (int n);
Real **createReal2DArray (int m, int n);
void transpose (Real **bt, Real **b, int m);
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);

void printMatrix(Real **matrix, int m)
{
  int i,j;
	for(i = 0; i < m; i++)
	{
		for(j = 0; j < m; j++)
		{
			printf("	%lf	", matrix[i][j]);
		}
		printf("\n");
	}
}

void printVector(Real *vec, int m)
{
	int i = 0;
	for (i = 0; i < m; i++)
	{
		printf("	%lf	", vec[i]);
	}
	printf("\n");
}
int main(int argc, char **argv )
{
  Real *diag, **b, **bt, *z;
  Real pi, h, umax;
  int i, j, n, m, nn;

  /* the total number of grid points in each spatial direction is (n+1) */
  /* the total number of degrees-of-freedom in each spatial direction is (n-1) */
  /* this version requires n to be a power of 2 */

	if( argc < 2 )
	{
    printf("need a problem size\n");
    return 1;
	}

  n  = atoi(argv[1]);
  m  = n-1;
  nn = 4*n;

  diag = createRealArray (m);
  b    = createReal2DArray (m,m);
  bt   = createReal2DArray (m,m);
  z    = createRealArray (nn);

  h    = 1./(Real)n;
  pi   = 4.*atan(1.);
  // Generates the eigenvalues and stores it in diagonal
  for (i=0; i < m; i++) {
    diag[i] = 2.*(1.-cos((i+1)*pi/(Real)n));
  }
  printf("The diagonal vector of eigenvalues \n");
  printVector(diag, m);

  // Vector vec = createVector(m);
  // vec->data = diag;
  for (j=0; j < m; j++) {
    for (i=0; i < m; i++) {
      b[j][i] = h*h;
    }
  }
  // printing b matrix for testing
  printf("The B matrix filled with h*h\n");
  printMatrix(b, m);
  #pragma omp parallel for schedule(static)
  for (j=0; j < m; j++) {
    fst_(b[j], &n, z, &nn);
  }
  printf("The B matrix after fst_\n");
  printMatrix(b, m);

  transpose (bt,b,m);
  printf("The bt matrix after transpose\n");
  printMatrix(bt, m);
  printf("\n");
  printVector(bt[0], m);

  for (i=0; i < m; i++) {
    fstinv_(bt[i], &n, z, &nn);
  }
 	printf("The bt matrix after fstinv_\n");
  printMatrix(bt, m);
  
  for (j=0; j < m; j++) {
    for (i=0; i < m; i++) {
      // printf("The expression: bt[j][i] = bt[j][i]/(diag[i]+diag[j]),\n bt[j][i] = %lf, diag[i] = %lf ,diag[j] = %lf \n", bt[j][i], diag[i], diag[j]);
      bt[j][i] = bt[j][i]/(diag[i]+diag[j]);
      printf("bt[j][i] = %lf (j = %d, i = %d) \n", bt[j][i], j, i);
    }
  }
  printf("The bt matrix something\n");
  printMatrix(bt, m);
  
  for (i=0; i < m; i++) {
    fst_(bt[i], &n, z, &nn);
  }
  printf("The bt matrix after fst_\n");
  printMatrix(bt, m);

  transpose (b,bt,m);
  printf("The b matrix after another transpose \n");
  printMatrix(b, m);

  for (j=0; j < m; j++) {
    fstinv_(b[j], &n, z, &nn);
  }
	printf("The b matrix after fstinv_ \n");
	printMatrix(b, m);

  umax = 0.0;
  for (j=0; j < m; j++) {
    for (i=0; i < m; i++) {
      if (b[j][i] > umax) umax = b[j][i];
    }
  }
  printf (" umax = %e \n",umax);

  return 0;
}

void transpose (Real **bt, Real **b, int m)
{
  int i, j;
  for (j=0; j < m; j++) {
    for (i=0; i < m; i++) {
      bt[j][i] = b[i][j];
    }
  }
}

void local_transpose (Real **bt, Real **b, int m)
{
  int i, j;
  for (j=0; j < m; j++) {
    for (i=0; i < m; i++) {
      bt[j][i] = b[i][j];
    }
  }
}

void transposeMPI(Real **bt, Real **b, int m)
{

}

Real *createRealArray (int n)
{
  Real *a;
  int i;
  a = (Real *)malloc(n*sizeof(Real));
  for (i=0; i < n; i++) {
    a[i] = 0.0;
  }
  return (a);
}

Real **createReal2DArray (int n1, int n2)
{
  int i, n;
  Real **a;
  a    = (Real **)malloc(n1   *sizeof(Real *));
  a[0] = (Real  *)malloc(n1*n2*sizeof(Real));
  for (i=1; i < n1; i++) {
    a[i] = a[i-1] + n2;
  }
  n = n1*n2;
  memset(a[0],0,n*sizeof(Real));
  return (a);
}
