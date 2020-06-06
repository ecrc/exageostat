/**
 *
 * @file codelet_sgemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon sgemm StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @generated s Sun Nov 25 16:59:44 2018
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_s.h"
#include "../include/starpu_exageostat.h"
/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 */
//SCODELETS_HEADER(exageostat_gemm)
void MORSE_TASK_sexageostat_gemm(const MORSE_option_t *options,
                      MORSE_enum transA, int transB,
                      int m, int n, int k, int nb,
                      float alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                                               const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                      float beta,  const MORSE_desc_t *C, int Cm, int Cn, int ldc)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_sgemm;
    void (*callback)(void*) = options->profiling ? cl_sgemm_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_R(A, Am, An);
    MORSE_ACCESS_R(B, Bm, Bn);
    MORSE_ACCESS_RW(C, Cm, Cn);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &transA,            sizeof(MORSE_enum),
        STARPU_VALUE,    &transB,            sizeof(MORSE_enum),
        STARPU_VALUE,    &m,                 sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_VALUE,    &k,                 sizeof(int),
        STARPU_VALUE,    &alpha,             sizeof(float),
        STARPU_R,        EXAGEOSTAT_RTBLKADDR(A, MorseRealFloat, Am, An),
        STARPU_VALUE,    &lda,               sizeof(int),
        STARPU_R,        EXAGEOSTAT_RTBLKADDR(B, MorseRealFloat, Bm, Bn),
        STARPU_VALUE,    &ldb,               sizeof(int),
        STARPU_VALUE,    &beta,              sizeof(float),
        STARPU_RW,       EXAGEOSTAT_RTBLKADDR(C, MorseRealFloat, Cm, Cn),
        STARPU_VALUE,    &ldc,               sizeof(int),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "sgemm",
#endif
        0);
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_sgemm_cpu_func(void *descr[], void *cl_arg)
{
    MORSE_enum transA;
    MORSE_enum transB;
    int m;
    int n;
    int k;
    float alpha;
    float *A;
    int lda;
    float *B;
    int ldb;
    float beta;
    float *C;
    int ldc;

    A = (float *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (float *)STARPU_MATRIX_GET_PTR(descr[1]);
    C = (float *)STARPU_MATRIX_GET_PTR(descr[2]);
    starpu_codelet_unpack_args(cl_arg, &transA, &transB, &m, &n, &k, &alpha, &lda, &ldb, &beta, &ldc);
    CORE_sgemm(transA, transB,
        m, n, k,
        alpha, A, lda,
        B, ldb,
        beta, C, ldc);
}

#ifdef CHAMELEON_USE_CUDA
static void cl_sgemm_cuda_func(void *descr[], void *cl_arg)
{
    MORSE_enum transA;
    MORSE_enum transB;
    int m;
    int n;
    int k;
    float alpha;
    const float *A;
    int lda;
    const float *B;
    int ldb;
    float beta;
    float *C;
    int ldc;

    A = (const float *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (const float *)STARPU_MATRIX_GET_PTR(descr[1]);
    C = (float *)STARPU_MATRIX_GET_PTR(descr[2]);
    starpu_codelet_unpack_args(cl_arg, &transA, &transB, &m, &n, &k, &alpha, &lda, &ldb, &beta, &ldc);

    RUNTIME_getStream( stream );

    CUDA_sgemm(
        transA, transB,
        m, n, k,
        &alpha, A, lda,
                B, ldb,
        &beta,  C, ldc,
        stream);

#ifndef STARPU_CUDA_ASYNC
    cudaStreamSynchronize( stream );
#endif

    return;
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS(sexageostat_gemm, 3, cl_sgemm_cpu_func, cl_sgemm_cuda_func, STARPU_CUDA_ASYNC)
