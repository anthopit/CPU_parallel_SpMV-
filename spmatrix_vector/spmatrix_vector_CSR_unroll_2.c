// Computes matrix-vector product. Matrix A is in row-major order
// i.e. A[i, j] is stored in i * COLS + j element of the vector.

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "wtime.h"
#include "mmio.h"
#include "matrix_utils.h"
const int ntimes=20;

inline double dmin ( double a, double b ) { return a < b ? a : b;}

inline double max ( double a, double b ) { return a > b ? a : b; }
inline double min ( double a, double b ) { return a < b ? a : b; }

void MatrixVectorCSRserial(int m, int n, const int* irp, const int* ja,const double* as, const double* x, double* restrict y)
{
    int i,j;
    for (i = 0; i < m; ++i) {
        double t=0.0;
        for (j = irp[i]; j < irp[i+1]; ++j) {
            t = t + as[j]*x[ja[j]];
        }
        y[i] = t;
    }
}

void MatrixVectorCSR_unrolled_2(int m, int n, const int* irp, const int* ja,const double* as, const double* x, double* restrict y)
{
  int i,j;
#pragma omp parallel for shared(x,y,irp,ja,as) private(i,j) schedule(static)
  for (i = 0; i < m; i++) {
    double t0 = 0.0;
    double t1 = 0.0;
    j = irp[i];
    for (; j < irp[i+1] - 1; j += 2) {
        t0 = t0 + as[j + 0]*x[ja[j + 0]];
        t1 = t1 + as[j + 1]*x[ja[j + 1]];
    }
    for (; j < irp[i+1]; j++) {
        t0 = t0 + as[j]*x[ja[j]];
    }
    y[i] = t0 + t1;
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
  //int maxnz, **ja;
  //double **as;
  
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
  double* y1 = (double*) malloc(sizeof(double)*n );
  int row;
  
  srand(12345);
  for ( row = 0; row < n; ++row) {
    x[row] = 100.0f * ((double) rand()) / RAND_MAX;
  }

    read_mtx_coo_csr(f, m, n, &nz, &irp, &ja, &as, mm_is_symmetric(matcode), mm_is_pattern(matcode));
//    print_csr_mtx_csr(m, n, nz, irp, ja, as);
//    print_csr_mtx_2D(m, n, nz, irp, ja, as);

    double tmlt = 1e100;
    for (int try=0; try < ntimes; try ++ ) {
        double t1 = wtime();
        MatrixVectorCSR_unrolled_2(m, n, irp, ja, as, x, y);
        double t2 = wtime();
        tmlt = dmin(tmlt,(t2-t1));
    }
    double mflops = (2.0e-6)*nz/tmlt;

    MatrixVectorCSRserial(m,n, irp,ja,as,x,y1);

    double error = 0.0;
    double diff = 0.0;
    for (int i = 0; i < m; ++i) {
        diff = y[i] - y1[i];
        if (diff < 0.0) diff = -diff;
        error = max(error, diff);
    }

#pragma omp parallel
    {
#pragma omp master
        {
            fprintf(stdout,"%d;%d;%d;%d;%lf;%lf;%lf\n",
                    m,n,nz,omp_get_num_threads(),tmlt,mflops,error);
        }
    }

  free(irp);
  free(ja);
  free(as);
  free(x);
  free(y);
  return 0;
}
