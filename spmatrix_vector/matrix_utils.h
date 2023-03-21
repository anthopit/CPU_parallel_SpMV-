//
// Created by anthony on 08/02/23.
//

#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include <stdbool.h>
int get_nz_mtx_symetric(int *I, int *J, int nz);
int read_mtx_coo_csr(FILE *f, int M, int N, int *nz, int **IRP_, int **JA_, double **AS_, bool is_symmetric, bool is_pattern);
int read_mtx_coo_ellpack(FILE *f, int M, int N, int *nz, int *MAXNZ, int ***JA_ELL_, double ***AS_ELL_, bool is_symmetric, bool is_pattern);
int read_mtx_coo_ellpack_t(FILE *f, int M, int N, int *nz, int *MAXNZ, int **JA_ELL_, double **AS_ELL_, bool is_symmetric, bool is_pattern);
void print_csr_mtx_csr(int M, int N, int nz, int *IRP, int *JA, double *AS);
void print_csr_mtx_2D(int M, int N, int nz, int *IRP, int *JA, double *AS);
void print_ellpack_mtx_ellpack(int M, int N, int *MAXNZ, int **JA_ELL, double **AS_ELL);
void print_ellpack_mtx_2D(int M, int N, int *MAXNZ, int **JA_ELL, double **AS_ELL);
#endif //MATRIX_UTILS_H

