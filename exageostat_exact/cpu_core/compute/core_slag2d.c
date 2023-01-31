/**
 *
 * Copyright (c) 2017-2023, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file core_slag2d.c
 *
 * Covert single descriptor to double descriptor with transpose.
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#include <lapacke_utils.h>
#include "../include/exageostatcore.h"

/***************************************************************************//**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 **/

void core_slag2d(int m, int n,
                 const float *A, int lda,
                 double* B, int ldb) {

    float *C;
    C = (float *) malloc(m * n * sizeof(float));
    LAPACKE_sge_trans(LAPACK_COL_MAJOR, m, n, A, lda, C, ldb);
    LAPACKE_slag2d(LAPACK_COL_MAJOR, n, m, C, ldb, B, ldb);
    free(C);

}
/***************************************************************************/


