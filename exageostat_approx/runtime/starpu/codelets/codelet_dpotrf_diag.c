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
 *  CHAMELEON codelets kernel
 *  CHAMELEON is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 1.2.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2022-11-09
 * @generated d Fri Dec  1 14:38:46 2017
 *
 **/
#include "chameleon_starpu.h"
//#include "chameleon/runtime/starpu/include/runtime_codelet_z.h"
#include <runtime/starpu/runtime_codelet_d.h>
#include "../include/starpu_exageostat_approx.h"

/**
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 **/

void CHAM_TASK_dpotrf_diag(const RUNTIME_option_t *options,
                           CHAM_enum uplo, int n, int nb,
                           const CHAM_desc_t *A, int Am, int An, int lda,
                           int iinfo) {
    (void) nb;
    struct starpu_codelet *codelet = &cl_dpotrf;
    void (*callback)(void *) = options->profiling ? cl_dpotrf_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
        CHAMELEON_ACCESS_RW(A, Am, An);CHAMELEON_END_ACCESS_DECLARATION;

    starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE, &uplo, sizeof(CHAM_enum),
            STARPU_VALUE, &n, sizeof(int),
            STARPU_RW, EXAGEOSTAT_RTBLKADDR(A, ChamRealDouble, Am, An),
            STARPU_VALUE, &lda, sizeof(int),
            STARPU_VALUE, &iinfo, sizeof(int),
            /* STARPU_SCRATCH,   options->ws_worker, */
            STARPU_PRIORITY, options->priority,
            STARPU_CALLBACK, callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "dpotrf_diag",
#endif
            0);
}


#if !defined(CHAMELEON_SIMULATION)

static void cl_dpotrf_diag_cpu_func(void *descr[], void *cl_arg) {
    CHAM_enum uplo;
    int n;
    double* A;
    int lda;
    int iinfo;
    int info = 0;

    A = (double* ) STARPU_MATRIX_GET_PTR(descr[0]);

    starpu_codelet_unpack_args(cl_arg, &uplo, &n, &lda, &iinfo);
    //keep the core function from Chameleon
    CORE_dpotrf(uplo, n, A, lda, &info);
}

#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(dpotrf_diag, cl_dpotrf_diag_cpu_func)