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
 * @file codelet_dpotrf.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 1.1.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2018-11-11
 * @generated d Fri Dec  1 14:38:46 2017
 *
 **/
#include "chameleon_starpu.h"
#include "runtime_codelet_d.h"
#include "../include/starpu_exageostat_approx.h"
/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 **/

void MORSE_TASK_dpotrf_diag(const MORSE_option_t *options,
                       MORSE_enum uplo, int n, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       int iinfo)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_dpotrf;
    void (*callback)(void*) = options->profiling ? cl_dpotrf_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_RW(A, Am, An);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &uplo,                      sizeof(MORSE_enum),
        STARPU_VALUE,    &n,                         sizeof(int),
        STARPU_RW,       EXAGEOSTAT_RTBLKADDR(A, MorseRealDouble, Am, An),
        STARPU_VALUE,    &lda,                       sizeof(int),
        STARPU_VALUE,    &iinfo,                     sizeof(int),
        /* STARPU_SCRATCH,   options->ws_worker, */
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "dpotrf_diag",
#endif
        0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_dpotrf_diag_cpu_func(void *descr[], void *cl_arg)
{
    MORSE_enum uplo;
    int n;
    double *A;
    int lda;
    int iinfo;
    int info = 0;

    A = (double *)STARPU_MATRIX_GET_PTR(descr[0]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &n, &lda, &iinfo);
    //keep the core function from Chameleon
    CORE_dpotrf(uplo, n, A, lda, &info);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(dpotrf_diag, 1, cl_dpotrf_diag_cpu_func)

