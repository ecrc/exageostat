/**
 *
 * @file codelet_strsm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon strsm StarPU codelet
 *
 * @version 1.2.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @generated s Sun Nov 25 16:59:45 2018
 *
 */
#include "chameleon_starpu.h"
//#include <chameleon/runtime/starpu/include/runtime_codelet_z.h>
#include <runtime/starpu/runtime_codelet_s.h>
#include "../include/starpu_exageostat.h"

/**
 *
 * @ingroup CORE_float
 *
 */
void EXAGEOSTAT_TASK_sexageostat_trsm(const RUNTIME_option_t *options,
                                     CHAM_enum side, CHAM_enum uplo, CHAM_enum transA, CHAM_enum diag,
                                     int m, int n, int nb,
                                     float alpha, const CHAM_desc_t *A, int Am, int An, int lda,
                                     const CHAM_desc_t *B, int Bm, int Bn, int ldb) {
    (void) nb;
    struct starpu_codelet *codelet = &cl_strsm;
    void (*callback)(void *) = options->profiling ? cl_strsm_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
        CHAMELEON_ACCESS_R(A, Am, An);
        CHAMELEON_ACCESS_RW(B, Bm, Bn);CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE, &side, sizeof(CHAM_enum),
            STARPU_VALUE, &uplo, sizeof(CHAM_enum),
            STARPU_VALUE, &transA, sizeof(CHAM_enum),
            STARPU_VALUE, &diag, sizeof(CHAM_enum),
            STARPU_VALUE, &m, sizeof(int),
            STARPU_VALUE, &n, sizeof(int),
            STARPU_VALUE, &alpha, sizeof(float),
            STARPU_R, EXAGEOSTAT_RTBLKADDR(A, ChamRealFloat, Am, An),
            STARPU_VALUE, &lda, sizeof(int),
            STARPU_RW, EXAGEOSTAT_RTBLKADDR(B, ChamRealFloat, Bm, Bn),
            STARPU_VALUE, &ldb, sizeof(int),
            STARPU_PRIORITY, options->priority,
            STARPU_CALLBACK, callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "strsm",
#endif
            0);
}


#if !defined(CHAMELEON_SIMULATION)

static void cl_strsm_cpu_func(void *descr[], void *cl_arg) {
    CHAM_enum side;
    CHAM_enum uplo;
    CHAM_enum transA;
    CHAM_enum diag;
    int m;
    int n;
    float alpha;
    float *A;
    int lda;
    float *B;
    int ldb;

    A = (float *) STARPU_MATRIX_GET_PTR(descr[0]);
    B = (float *) STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &side, &uplo, &transA, &diag, &m, &n, &alpha, &lda, &ldb);
    // printf("%s, %f, %f, %f, %f, %f, %f, %f, %f, \n", __func__,  A[0], A[1], A[2], A[3], B[0], B[1], B[2], B[3]);
    CORE_dtrsm(side, uplo,
               transA, diag,
               m, n,
               alpha, A, lda,
               B, ldb);
}

#ifdef CHAMELEON_USE_CUDA
static void cl_strsm_cuda_func(void *descr[], void *cl_arg)
{
    CHAM_enum side;
    CHAM_enum uplo;
    CHAM_enum transA;
    CHAM_enum diag;
    int m;
    int n;
    float alpha;
    const float *A;
    int lda;
    float *B;
    int ldb;

    A = (const float *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (float *)STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &side, &uplo, &transA, &diag, &m, &n, &alpha, &lda, &ldb);

    RUNTIME_getStream(stream);

    CUDA_strsm(
        side, uplo, transA, diag,
        m, n,
        &alpha, A, lda,
        B, ldb,
        stream);

#ifndef STARPU_CUDA_ASYNC
    cudaStreamSynchronize( stream );
#endif

    return;
}
#endif /* CHAMELEON_USE_CUDA */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS(sexageostat_trsm, cl_strsm_cpu_func, cl_strsm_cuda_func, STARPU_CUDA_ASYNC)
