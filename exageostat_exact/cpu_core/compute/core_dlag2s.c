/**
 *
 * Copyright (c) 2017-2023, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file core_dlag2s.c
 *
 * Covert double descriptor to single descriptor with transpose.
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#include <lapacke_utils.h>
#include "coreblas/lapacke.h"

/***************************************************************************//**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 **/

void core_dlag2s(int m, int n,
                 const double* A, int lda,
                 float *B, int ldb) {
    double* C;
    C = (double* ) malloc(m * n * sizeof(double));
    LAPACKE_dge_trans(LAPACK_COL_MAJOR, m, n, A, lda, C, ldb);
    LAPACKE_dlag2s(LAPACK_COL_MAJOR, n, m, C, ldb, B, ldb);
    free(C);
}
/***************************************************************************/