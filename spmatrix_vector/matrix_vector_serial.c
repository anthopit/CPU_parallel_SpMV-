// Computes matrix-vector product. Matrix A is in row-major order
// i.e. A[i, j] is stored in i * COLS + j element of the vector.

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "wtime.h"
#include "mmio.h"
#include "matrix_utils.h"


// Simple CPU implementation of matrix-vector product
void MatrixVector(int rows, int cols, const double* A, const double* x, double* restrict y) 
{
  int row,col, idx;
  double t;
  for (row = 0; row < rows; ++row) {
    t=0.0;
    for (col = 0; col < cols; ++col) {
      idx = row * cols + col;
      t = t + A[idx]*x[col];
    }
    y[row] = t;
  }
}

void MatrixVectorCSR(int rows, int cols, const int* maxnzr, const int* ja,const double* as, const double* x, double* restrict y) 
{
  int i,j;
  double t;
  for (i = 0; i < m; ++i) {
    t=0.0;
    for (j = irp(i); j < irp(i+1)-1; ++j) {
      t = t + as(j)*x(ja(j));
    }
    y(i) = t;
  }
}


void MatrixVectorELLPACK(int rows, int cols, const int* maxnz, const int** ja,const double** as, const double* x, double* restrict y) 
{
  int i,j;
  double t;
  for (i = 0; i < m; ++i) {
    t=0.0;
    for (j = 0; j < maxnz; ++j) {
      t = t + as[i][j]*x(ja[i][j]));
    }
    y(i) = t;
  }
}

int main(int argc, char** argv) 
{

  int ret_code;
  MM_typecode matcode;
  FILE *f;
  
  int m, n, nz;

  // CSR
  int *irp, *ja;
  double *as;

  // ELLPACK
  // int maxnz, **ja;
  // double **as;
  
  if (argc < 2)
  {
  	    fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
        exit(1);
  }
  else
  {
    if ((f = fopen(argv[1], "r")) == NULL){
       printf("Could not open file %s", argv[1]);
       exit(1);
     }
  }

  if (mm_read_banner(f, &matcode) != 0)
  {
      printf("Could not process Matrix Market banner.\n");
      exit(1);
  }

  /*  This is how one can screen matrix types if their application */
  /*  only supports a subset of the Matrix Market data types.      */

  if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
      mm_is_sparse(matcode) && mm_is_integer(matcode))
  {
      printf("Sorry, this application does not support ");
      printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
      exit(1);
  }

  /* find out size of sparse matrix .... */

  if ((ret_code = mm_read_mtx_crd_size(f, &m, &n, &nz)) !=0)
      exit(1);
      
      
  double* x = (double*) malloc(sizeof(double)*n );
  double* y = (double*) malloc(sizeof(double)*n );
  int row;
  
  srand(12345);
  for ( row = 0; row < n; ++row) {
    x[row] = 100.0f * ((double) rand()) / RAND_MAX;      
  }    
      

//   read_mtx_coo_ellpack(f, m, n, nz, &MAXNZ, &jaL, &as, mm_is_symmetric(matcode));
//   print_ellpack_mtx_ellpack(m, n, &MAXNZ, ja, as);
//   print_ellpack_mtx_2D(m, n, &MAXNZ, ja, as);

    read_mtx_coo_csr(f, m, N, nz, &irp, &ja, &as, mm_is_symmetric(matcode));
//    print_csr_mtx_csr(m, n, nz, irp, ja, as);
//    print_csr_mtx_2D(m, n, nz, irp, ja, as);
  
  double t1 = wtime();
  MatrixVectorCSR(m, n, &maxnzr, ja, as, x, y);
  double t2 = wtime();
  double tmlt = (t2-t1);
  double mflops = (2.0e-6)*n*m/tmlt;
  
  fprintf(stdout,"Matrix-Vector product of size %d x %d with %d non zero element with 1 thread: time %lf  MFLOPS %lf \n",
	  n,m,nz,tmlt,mflops);

  free(irp);
  free(ja);
  free(as);
  free(x);
  free(y);
  return 0;
}

