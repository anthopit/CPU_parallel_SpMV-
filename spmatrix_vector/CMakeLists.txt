# CSR

add_executable(spmatrix_vector_CSR spmatrix_vector_CSR.c mmio.h mmio.c matrix_utils.h matrix_utils.c wtime.c)
add_executable(spmatrix_vector_CSR_unroll_2 spmatrix_vector_CSR_unroll_2.c mmio.h mmio.c matrix_utils.h matrix_utils.c wtime.c)
add_executable(spmatrix_vector_CSR_unroll_4 spmatrix_vector_CSR_unroll_4.c mmio.h mmio.c matrix_utils.h matrix_utils.c wtime.c)
add_executable(spmatrix_vector_CSR_unroll_8 spmatrix_vector_CSR_unroll_8.c mmio.h mmio.c matrix_utils.h matrix_utils.c wtime.c) 

# ELLPACK

add_executable(spmatrix_vector_ELL spmatrix_vector_ELL.c mmio.h mmio.c matrix_utils.h matrix_utils.c wtime.c)
add_executable(spmatrix_vector_ELL_unroll_2 spmatrix_vector_ELL_unroll_2.c mmio.h mmio.c matrix_utils.h matrix_utils.c wtime.c)
add_executable(spmatrix_vector_ELL_unroll_8 spmatrix_vector_ELL_unroll_8.c mmio.h mmio.c matrix_utils.h matrix_utils.c wtime.c)
add_executable(spmatrix_vector_ELL_unroll_block_2 spmatrix_vector_ELL_unroll_block_2.c mmio.h mmio.c matrix_utils.h matrix_utils.c wtime.c)
add_executable(spmatrix_vector_ELL_unroll_block_2_ro spmatrix_vector_ELL_unroll_block_2_ro.c mmio.h mmio.c matrix_utils.h matrix_utils.c wtime.c)
add_executable(spmatrix_vector_ELL_unroll_block_8 spmatrix_vector_ELL_unroll_block_8.c mmio.h mmio.c matrix_utils.h matrix_utils.c wtime.c)

