//
// Created by anthony on 08/02/23.
//

#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "mmio.h"
#include "matrix_utils.h"

//------------------------------------- Read Matrix Market file -------------------------------------//

int get_nz_mtx_symetric(int *I, int *J, int nz) {
    int nz_sym = 0;
    for (int k=0; k<nz; k++)
    {
        if (I[k] != J[k]) {
            nz_sym += 2;
        } else {
            nz_sym += 1;
        }
    }
     return nz = nz_sym;
}

int read_mtx_coo_csr(FILE *f, int M, int N, int *nz, int **IRP_, int **JA_, double **AS_, bool is_symmetric, bool is_pattern) {

    int *J, *I;
    double *val;

    int *IRP, *JA;
    double *AS;

    int *IRP_temp, *IRP_temp2;

    IRP = (int *) calloc(M+1 , sizeof(int));
    IRP_temp = (int *) calloc(M+1 ,sizeof(int));
    if (is_symmetric) {
        IRP_temp2 = (int *) calloc(M+1 ,sizeof(int));
    }

    I = (int *) malloc(*nz * sizeof(int));
    J = (int *) malloc(*nz * sizeof(int));
    val = (double *) malloc(*nz * sizeof(double));

    if (is_symmetric) {
        if (is_pattern) {
            for (int i=0; i< *nz; i++)
            {
                fscanf(f, "%d %d\n", &I[i], &J[i]);
                val[i] = 1.0;
                if (I[i] != J[i]) {
                    IRP[I[i]]++;
                    IRP[J[i]]++;
                } else {
                    IRP[I[i]]++;
                }
                I[i]--;  /* adjust from 1-based to 0-based */
                J[i]--;
            }
        } else {
            for (int i=0; i< *nz; i++)
            {
                fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
                if (I[i] != J[i]) {
                    IRP[I[i]]++;
                    IRP[J[i]]++;
                } else {
                    IRP[I[i]]++;
                }
                I[i]--;  /* adjust from 1-based to 0-based */
                J[i]--;
            }
        }

    } else {
        if (is_pattern) {
            for (int i=0; i< *nz; i++)
            {
                fscanf(f, "%d %d\n", &I[i], &J[i]);
                val[i] = 1.0;
                IRP[I[i]]++;
                I[i]--;  /* adjust from 1-based to 0-based */
                J[i]--;
            }
        } else {
            for (int i=0; i< *nz; i++)
            {
                fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
                IRP[I[i]]++;
                I[i]--;  /* adjust from 1-based to 0-based */
                J[i]--;
            }
        }
    }

    int nz_sym = get_nz_mtx_symetric(I, J, *nz);

    if (is_symmetric) {
        IRP[M] = nz_sym;
    } else {
        IRP[M] = *nz;
    }

    for (int i = 1; i < M; ++i) {
        IRP[i] += IRP[i-1];
    }

    memcpy(IRP_temp, IRP, sizeof(int) * M+1);
    if (is_symmetric) {
        memcpy(IRP_temp2, IRP, sizeof(int) * M+1);
    }

    if (is_symmetric) {
        JA = (int *) calloc(nz_sym , sizeof(int));
        AS = (double *) calloc(nz_sym , sizeof(double));
    } else {
        JA = (int *) calloc(*nz , sizeof(int));
        AS = (double *) calloc(*nz , sizeof(double));
    }


    *IRP_ = IRP;
    *JA_ = JA;
    *AS_ = AS;

    if (is_symmetric) {
        for (int j = *nz - 1; j >= 0; j--) {
            if (I[j] != J[j]) {
                JA[IRP_temp2[J[j]+1]-1] = I[j];
                AS[IRP_temp2[J[j]+1]-1] = val[j];
                IRP_temp2[J[j]+1] --;
            }
        }
        for (int j = 0; j < *nz; ++j) {
            JA[IRP_temp[I[j]]] = J[j];
            AS[IRP_temp[I[j]]] = val[j];
            IRP_temp[I[j]] ++;
        }
    } else {
        for (int j = 0; j < *nz; ++j) {
            JA[IRP_temp[I[j]]] = J[j];
            AS[IRP_temp[I[j]]] = val[j];
            IRP_temp[I[j]] ++;
        }
    }
    
    if (is_symmetric) {
    	*nz = nz_sym;
    }

    return 0;
}

int read_mtx_coo_ellpack_t(FILE *f, int M, int N, int *nz, int *MAXNZ, int **JA_ELL_, double **AS_ELL_, bool is_symmetric, bool is_pattern){

    int *J, *I;
    double *val;

    int *JA_ELL;
    double *AS_ELL;

    int *MAXNZ_temp;

    MAXNZ_temp = (int *) calloc(M , sizeof(int));

    I = (int *) calloc(*nz, sizeof(int));
    J = (int *) calloc(*nz, sizeof(int));
    val = (double *) calloc(*nz, sizeof(double));

    if (is_symmetric) {
        if (is_pattern) {
            for (int i=0; i< *nz; i++)
            {
                fscanf(f, "%d %d\n", &I[i], &J[i]);
                val[i] = 1.0;
                if (I[i] != J[i]) {
                    MAXNZ_temp[I[i]-1] += 2;
                } else {
                    MAXNZ_temp[I[i]-1]++;
                }
                I[i]--;  /* adjust from 1-based to 0-based */
                J[i]--;
            }
        } else {
            for (int i=0; i< *nz; i++)
            {
                fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
                if (I[i] != J[i]) {
                    MAXNZ_temp[I[i]-1] += 2;
                } else {
                    MAXNZ_temp[I[i]-1]++;
                }
                I[i]--;  /* adjust from 1-based to 0-based */
                J[i]--;
            }
        }
    } else {
        if (is_pattern) {
            for (int i=0; i< *nz; i++)
            {
                fscanf(f, "%d %d\n", &I[i], &J[i]);
                val[i] = 1.0;
                MAXNZ_temp[I[i]-1]++;
                I[i]--;  /* adjust from 1-based to 0-based */
                J[i]--;
            }
        } else {
            for (int i = 0; i < *nz; i++) {
                fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
                MAXNZ_temp[I[i] - 1]++;
                I[i]--;  /* adjust from 1-based to 0-based */
                J[i]--;
            }
        }
    }

    *MAXNZ = MAXNZ_temp[0];
    for (int i = 1; i < M; i++) {
        if (MAXNZ_temp[i] > *MAXNZ) {
            *MAXNZ = MAXNZ_temp[i];
        }
    }

    for (int i = 0; i < M; ++i) {
        MAXNZ_temp[i] = *MAXNZ;
    }

    JA_ELL = (int *) calloc(M * *MAXNZ, sizeof(int *));
    AS_ELL = (double *) calloc(M * *MAXNZ, sizeof(double *));

    *JA_ELL_ = JA_ELL;
    *AS_ELL_ = AS_ELL;

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < *MAXNZ; j++) {
            JA_ELL[i * *MAXNZ + j] = 0;
            AS_ELL[i * *MAXNZ + j] = 0.0;
        }
    }


    if (is_symmetric) {
        for (int i=0; i< *nz; i++)
        {
            if (I[i] != J[i]) {
                JA_ELL[I[i] * *MAXNZ + MAXNZ_temp[I[i]]-*MAXNZ] = J[i];
                AS_ELL[I[i] * *MAXNZ + MAXNZ_temp[I[i]]-*MAXNZ] = val[i];
                MAXNZ_temp[I[i]]++;
                JA_ELL[J[i] * *MAXNZ + MAXNZ_temp[J[i]]-*MAXNZ] = I[i];
                AS_ELL[J[i] * *MAXNZ + MAXNZ_temp[J[i]]-*MAXNZ] = val[i];
                MAXNZ_temp[J[i]]++;
            } else {
                JA_ELL[I[i] * *MAXNZ + MAXNZ_temp[I[i]]-*MAXNZ] = J[i];
                AS_ELL[I[i] * *MAXNZ +MAXNZ_temp[I[i]]-*MAXNZ] = val[i];
                MAXNZ_temp[I[i]]++;
            }
        }
    } else {
        for (int i=0; i< *nz; i++)
        {
            JA_ELL[I[i] * *MAXNZ + MAXNZ_temp[I[i]]-*MAXNZ] = J[i];
            AS_ELL[I[i] * *MAXNZ + MAXNZ_temp[I[i]]-*MAXNZ] = val[i];
            MAXNZ_temp[I[i]]++;
        }
    }

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < *MAXNZ; j++) {
            if (JA_ELL[i * *MAXNZ + j] == 0) {
                JA_ELL[i * *MAXNZ + j] = JA_ELL[i * *MAXNZ + j -1];
            }
        }
    }

    if (is_symmetric) {
    	*nz = get_nz_mtx_symetric(I, J, *nz);
    }

    return 0;
}

int read_mtx_coo_ellpack(FILE *f, int M, int N, int *nz, int *MAXNZ, int ***JA_ELL_, double ***AS_ELL_, bool is_symmetric, bool is_pattern){

    int *J, *I;
    double *val;

    int **JA_ELL;
    double **AS_ELL;

    int *MAXNZ_temp;

    MAXNZ_temp = (int *) calloc(M , sizeof(int));

    I = (int *) calloc(*nz, sizeof(int));
    J = (int *) calloc(*nz, sizeof(int));
    val = (double *) calloc(*nz, sizeof(double));

    if (is_symmetric) {
        if (is_pattern) {
            for (int i=0; i< *nz; i++)
            {
                fscanf(f, "%d %d\n", &I[i], &J[i]);
                val[i] = 1.0;
                if (I[i] != J[i]) {
                    MAXNZ_temp[I[i]-1] += 2;
                } else {
                    MAXNZ_temp[I[i]-1]++;
                }
                I[i]--;  /* adjust from 1-based to 0-based */
                J[i]--;
            }
        } else {
            for (int i=0; i< *nz; i++)
            {
                fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
                if (I[i] != J[i]) {
                    MAXNZ_temp[I[i]-1] += 2;
                } else {
                    MAXNZ_temp[I[i]-1]++;
                }
                I[i]--;  /* adjust from 1-based to 0-based */
                J[i]--;
            }
        }
    } else {
        if (is_pattern) {
            for (int i=0; i< *nz; i++)
            {
                fscanf(f, "%d %d\n", &I[i], &J[i]);
                val[i] = 1.0;
                MAXNZ_temp[I[i]-1]++;
                I[i]--;  /* adjust from 1-based to 0-based */
                J[i]--;
            }
        } else {
            for (int i = 0; i < *nz; i++) {
                fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
                MAXNZ_temp[I[i] - 1]++;
                I[i]--;  /* adjust from 1-based to 0-based */
                J[i]--;
            }
        }
    }

    *MAXNZ = MAXNZ_temp[0];
    for (int i = 1; i < M; i++) {
        if (MAXNZ_temp[i] > *MAXNZ) {
            *MAXNZ = MAXNZ_temp[i];
        }
    }

    for (int i = 0; i < M; ++i) {
        MAXNZ_temp[i] = *MAXNZ;
    }

    JA_ELL = (int **) calloc(M , sizeof(int *));
    AS_ELL = (double **) calloc(M , sizeof(double *));
    for (int i=0; i<M; i++) {
        JA_ELL[i] = (int *) calloc(*MAXNZ , sizeof(int));
        AS_ELL[i] = (double *) calloc(*MAXNZ , sizeof(double));
    }

    *JA_ELL_ = JA_ELL;
    *AS_ELL_ = AS_ELL;

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < *MAXNZ; j++) {
            JA_ELL[i][j] = 0;
            AS_ELL[i][j] = 0.0;
        }
    }


    if (is_symmetric) {
        for (int i=0; i< *nz; i++)
        {
            if (I[i] != J[i]) {
                JA_ELL[I[i]][MAXNZ_temp[I[i]]-*MAXNZ] = J[i];
                AS_ELL[I[i]][MAXNZ_temp[I[i]]-*MAXNZ] = val[i];
                MAXNZ_temp[I[i]]++;
                JA_ELL[J[i]][MAXNZ_temp[J[i]]-*MAXNZ] = I[i];
                AS_ELL[J[i]][MAXNZ_temp[J[i]]-*MAXNZ] = val[i];
                MAXNZ_temp[J[i]]++;
            } else {
                JA_ELL[I[i]][MAXNZ_temp[I[i]]-*MAXNZ] = J[i];
                AS_ELL[I[i]][MAXNZ_temp[I[i]]-*MAXNZ] = val[i];
                MAXNZ_temp[I[i]]++;
            }
        }
    } else {
        for (int i=0; i< *nz; i++)
        {
            JA_ELL[I[i]][MAXNZ_temp[I[i]]-*MAXNZ] = J[i];
            AS_ELL[I[i]][MAXNZ_temp[I[i]]-*MAXNZ] = val[i];
            MAXNZ_temp[I[i]]++;
        }
    }

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < *MAXNZ; j++) {
            if (JA_ELL[i][j] == 0) {
                JA_ELL[i][j] = JA_ELL[i][j-1];
            }
        }
    }

    if (is_symmetric) {
        *nz = get_nz_mtx_symetric(I, J, *nz);
    }

    return 0;
}

//------------------------------------- Print matrix -------------------------------------//

 void print_csr_mtx_csr(int M, int N, int nz, int *IRP, int *JA, double *AS) {
    printf("M: %d, N: %d, nz: %d\n", M, N, nz);
     printf("IRP:\n");
    for (int j = 0; j < M+1; ++j) {
        printf("%d ", IRP[j]);
    }
    printf("\n");
    printf("JA:\n");
    for (int i = 0; i < 46; ++i) {
        printf("%d ", JA[i]);
    }
    printf("\n");
    printf("AS:\n");
    for (int j = 0; j < 46; ++j) {
        printf("%lg ", AS[j]);
    }
    printf("\n");
}

void print_csr_mtx_2D(int M, int N, int nz, int *IRP, int *JA, double *AS) {
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            int found = 0;
            for (int k = IRP[i]; k < IRP[i + 1]; k++) {
                if (JA[k] == j) {
//                    printf("%.2lg ", AS[k]);
                    printf("x ");
                    found = 1;
                    break;
                }
            }
            if (!found) {
                printf("- ");
            }
        }
        printf("\n");
    }
}


void print_ellpack_mtx_ellpack(int M, int N, int *MAXNZ, int **JA_ELL, double **AS_ELL) {
    printf("M: %d, N: %d \n", M, N);
    printf("MAXNZ: %d\n", *MAXNZ);
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < *MAXNZ; j++) {
            printf("%d ", JA_ELL[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < *MAXNZ; j++) {
            printf("%lf ", AS_ELL[i][j]);
        }
        printf("\n");
    }
}

void print_ellpack_mtx_2D(int M, int N, int *MAXNZ, int **JA_ELL, double **AS_ELL) {
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            int found = 0;
            for (int k = 0; k < *MAXNZ; k++) {
                if (JA_ELL[i][k] == j) {
                    printf("x ", AS_ELL[i][k]);
                    found = 1;
                    break;
                }
            }
            if (!found) {
                printf("- ");
            }
        }
        printf("\n");
    }
}
