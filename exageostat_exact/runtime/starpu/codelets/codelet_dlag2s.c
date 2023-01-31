/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014, 2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file codelet_dlag2s.c
 *
 *  CHAMELEON codelets kernel
 *  CHAMELEON is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions mixed zc -> ds
 *
 **/
#include "chameleon_starpu.h"
#include "runtime_codelet_d.h"
#include "../include/starpu_exageostat.h"
/**
 *
 * @ingroup CORE_ChamRealDouble
 *
 **/
void CHAMELEON_TASK_dlag2s(const RUNTIME_option_t *options,
                       int m, int n, int nb,
                       const CHAM_desc_t *A, int Am, int An, int lda,
                       const CHAM_desc_t *B, int Bm, int Bn, int ldb)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_lag2s;
    //void (*callback)(void*) = options->profiling ? cl_dlag2s_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_W(B, Bm, Bn);
    CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &m,                 sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_R,         EXAGEOSTAT_RTBLKADDR(A, ChamRealDouble, Am, An),
        STARPU_VALUE,    &lda,               sizeof(int),
        STARPU_W,         EXAGEOSTAT_RTBLKADDR(B, ChamRealFloat, Bm, Bn),
        STARPU_VALUE,    &ldb,               sizeof(int),
        STARPU_PRIORITY,  options->priority,
     //   STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "dlag2s",
#endif
        0);
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_dlag2s_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    ChamRealDouble *A;
    int lda;
    ChamRealFloat *B;
    int ldb;

    A = (ChamRealDouble *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (ChamRealFloat *)STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &lda, &ldb);
    CORE_dlag2s( m, n, A, lda, B, ldb);
}
#endif /* !defined(CHAMELEON_SIMULATION) */
/*
void CHAMELEON_TASK_dlag2s(const RUNTIME_option_t *options,
                       int m, int n, int nb,
                       const CHAM_desc_t *A, int Am, int An, int lda,
                       const CHAM_desc_t *B, int Bm, int Bn, int ldb)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_dlag2s;
    void (*callback)(void*) = options->profiling ? cl_dlag2s_callback : NULL;

    if ( chameleon_desc_islocal( A, Am, An ) ||
         chameleon_desc_islocal( B, Bm, Bn ) )
    {
        starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE,    &m,                 sizeof(int),
            STARPU_VALUE,    &n,                 sizeof(int),
            STARPU_R,         RTBLKADDR(A, ChamRealFloat, Am, An),
            STARPU_VALUE,    &lda,               sizeof(int),
            STARPU_W,         RTBLKADDR(B, ChamRealDouble, Bm, Bn),
            STARPU_VALUE,    &ldb,               sizeof(int),
            STARPU_PRIORITY,  options->priority,
            STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "dlag2s",
#endif
            0);
    }
}
*/


/*
 * Codelet definition
 */
//CODELETS_CPU(dlag2s, 1, cl_dlag2s_cpu_func)
