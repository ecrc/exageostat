/**
 *
 * Copyright (c) 2017-2020, King Abdullah University of Science and Technology
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
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2019-01-16
 *
 **/
#include "coreblas/lapacke.h"
#include "coreblas.h"
#include "../include/exageostatcore.h"
/***************************************************************************//**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 **/

void core_slag2d(int m, int n,
        const float *A, int lda,
        double *B, int ldb)
{

    float *C;
    C = (float*) malloc(m * n * sizeof(float));
    LAPACKE_sge_trans(LAPACK_COL_MAJOR, m, n, A, lda, C, ldb);
    LAPACKE_slag2d(LAPACK_COL_MAJOR, n, m, C, ldb, B, ldb);
    free(C);

}
/***************************************************************************/


