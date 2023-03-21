/*
*   Matrix Market I/O example program
*
*   Read a real (non-complex) sparse matrix from a Matrix Market (v. 2.0) file.
*   and copies it to stdout.  This porgram does nothing useful, but
*   illustrates common usage of the Matrix Matrix I/O routines.
*   (See http://math.nist.gov/MatrixMarket for details.)
*
*   Usage:  a.out [filename] > output
*
*
*   NOTES:
*
*   1) Matrix Market files are always 1-based, i.e. the index of the first
*      element of a matrix is (1,1), not (0,0) as in C.  ADJUST THESE/
*      OFFSETS ACCORDINGLY offsets accordingly when reading and writing
*      to files.
*
*   2) ANSI C requires one to use the "l" format modifier when reading
*      double precision floating point numbers in scanf() and
*      its variants.  For example, use "%lf", "%lg", or "%le"
*      when reading doubles, otherwise errors will occur.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mmio.h"
#include "matrix_utils.h"


int main(int argc, char *argv[])
{

    clock_t start, end;
    double elapsed;
    
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;

    // CSR
    int *IRP, *JA;
    double *AS;

    // ELLPACK
    int MAXNZ, **JA_ELL;
    double **AS_ELL;
    
    start = clock();


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

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);

//   read_mtx_coo_ellpack(f, M, N, nz, &MAXNZ, &JA_ELL, &AS_ELL, mm_is_symmetric(matcode));
//   print_ellpack_mtx_ellpack(M, N, &MAXNZ, JA_ELL, AS_ELL);
//   print_ellpack_mtx_2D(M, N, &MAXNZ, JA_ELL, AS_ELL);

    read_mtx_coo_csr(f, M, N, nz, &IRP, &JA, &AS, mm_is_symmetric(matcode));
//    print_csr_mtx_csr(M, N, nz, IRP, JA, AS);
//    print_csr_mtx_2D(M, N, nz, IRP, JA, AS);


    if (f !=stdin) fclose(f);

    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    
    printf("The elapsed time is %f seconds.\n", elapsed);

    return 0;
}


