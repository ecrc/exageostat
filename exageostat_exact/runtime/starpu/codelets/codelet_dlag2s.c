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
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
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
 * @ingroup CORE_MorseRealDouble
 *
 **/
void MORSE_TASK_dlag2s(const MORSE_option_t *options,
                       int m, int n, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_dlag2s;
    //void (*callback)(void*) = options->profiling ? cl_dlag2s_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_R(A, Am, An);
    MORSE_ACCESS_W(B, Bm, Bn);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &m,                 sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_R,         EXAGEOSTAT_RTBLKADDR(A, MorseRealDouble, Am, An),
        STARPU_VALUE,    &lda,               sizeof(int),
        STARPU_W,         EXAGEOSTAT_RTBLKADDR(B, MorseRealFloat, Bm, Bn),
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
    MorseRealDouble *A;
    int lda;
    MorseRealFloat *B;
    int ldb;

    A = (MorseRealDouble *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (MorseRealFloat *)STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &lda, &ldb);
    CORE_dlag2s( m, n, A, lda, B, ldb);
}
#endif /* !defined(CHAMELEON_SIMULATION) */
/*
void MORSE_TASK_dlag2s(const MORSE_option_t *options,
                       int m, int n, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_dlag2s;
    void (*callback)(void*) = options->profiling ? cl_dlag2s_callback : NULL;

    if ( morse_desc_islocal( A, Am, An ) ||
         morse_desc_islocal( B, Bm, Bn ) )
    {
        starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE,    &m,                 sizeof(int),
            STARPU_VALUE,    &n,                 sizeof(int),
            STARPU_R,         RTBLKADDR(A, MorseRealFloat, Am, An),
            STARPU_VALUE,    &lda,               sizeof(int),
            STARPU_W,         RTBLKADDR(B, MorseRealDouble, Bm, Bn),
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
