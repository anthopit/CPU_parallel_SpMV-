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

void MatrixVectorELLSerial(int m, int n, const int* maxnz, int** ja, double** as, const double* x, double* restrict y)
{
    int i,j;
    for (i = 0; i < m; ++i) {
        double t=0.0;
        for (j = 0; j < *maxnz; ++j) {
            t = t + as[i][j]*x[ja[i][j]];
        }
        y[i] = t;
    }
}

void innerMatrixVector(int rows, const int* maxnz, int** ja, double** as, int ncallc, const double* x,
                       double beta, double* restrict y)
{
    int row,col;
    double t0, t1, t2, t3, t4, t5, t6, t7;

#pragma omp parallel for shared(x,y,maxnz,ja,as,rows,ncallc) private(row,col,t0, t1, t2, t3, t4, t5, t6, t7 ) schedule(static)
    for (row = 0; row < rows - rows%8; row += 8) {
        if (beta == 0.0) {
            t0 = beta;
            t1 = beta;
            t2 = beta;
            t3 = beta;
            t4 = beta;
            t5 = beta;
            t6 = beta;
            t7 = beta;
        } else {
            t0 = beta *  y[row+0];
            t1 = beta *  y[row+1];
            t2 = beta *  y[row+2];
            t3 = beta *  y[row+3];
            t4 = beta *  y[row+4];
            t5 = beta *  y[row+5];
            t6 = beta *  y[row+6];
            t7 = beta *  y[row+7];
        }

        for (col = 0; col < *maxnz - *maxnz%8 ; col+=8) {
            t0 += as[row+0][col+0] * x[ja[row+0][col+0]] + as[row+0][col+1] * x[ja[row+0][col+1]]
                     +as[row+0][col+2] * x[ja[row+0][col+2]] + as[row+0][col+3] * x[ja[row+0][col+3]]
                    +as[row+0][col+4] * x[ja[row+0][col+4]] + as[row+0][col+5] * x[ja[row+0][col+5]]
                  +as[row+0][col+6] * x[ja[row+0][col+6]] + as[row+0][col+7] * x[ja[row+0][col+7]];
            t1 += as[row+1][col+0] * x[ja[row+1][col+0]] + as[row+1][col+1] * x[ja[row+1][col+1]]
                  +as[row+1][col+2] * x[ja[row+1][col+2]] + as[row+1][col+3] * x[ja[row+1][col+3]]
                   +as[row+1][col+4] * x[ja[row+1][col+4]] + as[row+1][col+5] * x[ja[row+1][col+5]]
                    +as[row+1][col+6] * x[ja[row+1][col+6]] + as[row+1][col+7] * x[ja[row+1][col+7]];
            t2 += as[row+2][col+0] * x[ja[row+2][col+0]] + as[row+2][col+1] * x[ja[row+2][col+1]]
                    +as[row+2][col+2] * x[ja[row+2][col+2]] + as[row+2][col+3] * x[ja[row+2][col+3]]
                  +as[row+2][col+4] * x[ja[row+2][col+4]] + as[row+2][col+5] * x[ja[row+2][col+5]]
                +as[row+2][col+6] * x[ja[row+2][col+6]] + as[row+2][col+7] * x[ja[row+2][col+7]];
            t3 += as[row+3][col+0] * x[ja[row+3][col+0]] + as[row+3][col+1] * x[ja[row+3][col+1]]
                  +as[row+3][col+2] * x[ja[row+3][col+2]] + as[row+3][col+3] * x[ja[row+3][col+3]]
                   +as[row+3][col+4] * x[ja[row+3][col+4]] + as[row+3][col+5] * x[ja[row+3][col+5]]
                    +as[row+3][col+6] * x[ja[row+3][col+6]] + as[row+3][col+7] * x[ja[row+3][col+7]];
            t4 += as[row+4][col+0] * x[ja[row+4][col+0]] + as[row+4][col+1] * x[ja[row+4][col+1]]
                    +as[row+4][col+2] * x[ja[row+4][col+2]] + as[row+4][col+3] * x[ja[row+4][col+3]]
                  +as[row+4][col+4] * x[ja[row+4][col+4]] + as[row+4][col+5] * x[ja[row+4][col+5]]
                +as[row+4][col+6] * x[ja[row+4][col+6]] + as[row+4][col+7] * x[ja[row+4][col+7]];
            t5 += as[row+5][col+0] * x[ja[row+5][col+0]] + as[row+5][col+1] * x[ja[row+5][col+1]]
                    +as[row+5][col+2] * x[ja[row+5][col+2]] + as[row+5][col+3] * x[ja[row+5][col+3]]
                     +as[row+5][col+4] * x[ja[row+5][col+4]] + as[row+5][col+5] * x[ja[row+5][col+5]]
                        +as[row+5][col+6] * x[ja[row+5][col+6]] + as[row+5][col+7] * x[ja[row+5][col+7]];
            t6 += as[row+6][col+0] * x[ja[row+6][col+0]] + as[row+6][col+1] * x[ja[row+6][col+1]]
                    +as[row+6][col+2] * x[ja[row+6][col+2]] + as[row+6][col+3] * x[ja[row+6][col+3]]
                     +as[row+6][col+4] * x[ja[row+6][col+4]] + as[row+6][col+5] * x[ja[row+6][col+5]]
                        +as[row+6][col+6] * x[ja[row+6][col+6]] + as[row+6][col+7] * x[ja[row+6][col+7]];
            t7 += as[row+7][col+0] * x[ja[row+7][col+0]] + as[row+7][col+1] * x[ja[row+7][col+1]]
                    +as[row+7][col+2] * x[ja[row+7][col+2]] + as[row+7][col+3] * x[ja[row+7][col+3]]
                     +as[row+7][col+4] * x[ja[row+7][col+4]] + as[row+7][col+5] * x[ja[row+7][col+5]]
                        +as[row+7][col+6] * x[ja[row+7][col+6]] + as[row+7][col+7] * x[ja[row+7][col+7]];
        }
        for (col = *maxnz - *maxnz%8; col < *maxnz; col++) {
            t0 += as[row+0][col] * x[ja[row+0][col]];
            t1 += as[row+1][col] * x[ja[row+1][col]];
            t2 += as[row+2][col] * x[ja[row+2][col]];
            t3 += as[row+3][col] * x[ja[row+3][col]];
            t4 += as[row+4][col] * x[ja[row+4][col]];
            t5 += as[row+5][col] * x[ja[row+5][col]];
            t6 += as[row+6][col] * x[ja[row+6][col]];
            t7 += as[row+7][col] * x[ja[row+7][col]];
        }
        y[row+0] = t0;
        y[row+1] = t1;
        y[row+2] = t2;
        y[row+3] = t3;
        y[row+4] = t4;
        y[row+5] = t5;
        y[row+6] = t6;
        y[row+7] = t7;
    }

    for (row = rows - rows%8; row < rows; row++) {
        double t=0.0;
        for (col = 0; col < *maxnz; col++) {
            t += as[row][col] * x[ja[row][col]];
        }
        y[row]=t;
    }
}


int main(int argc, char** argv) 
{

  int ret_code;
  MM_typecode matcode;
  FILE *f;
  
  int m, n, nz;

  // ELLPACK
  int maxnz, **ja;
  double **as;
  
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

   read_mtx_coo_ellpack(f, m, n, &nz, &maxnz, &ja, &as, mm_is_symmetric(matcode), mm_is_pattern(matcode));
//   print_ellpack_mtx_ellpack(m, n, &maxnz, ja, as);
//   print_ellpack_mtx_2D(m, n, &maxnz, ja, as);

    double tmlt = 1e100;
    for (int try=0; try < ntimes; try ++ ) {
        double t1 = wtime();
        innerMatrixVector(m, &maxnz, ja, as, n, x, 0.0, y);
        double t2 = wtime();
        tmlt = dmin(tmlt,(t2-t1));
    }
    double mflops = (2.0e-6)*nz/tmlt;


    MatrixVectorELLSerial(m, n, &maxnz, ja, as, x, y1);

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

  free(ja);
  free(as);
  free(x);
  free(y);
  return 0;
}
