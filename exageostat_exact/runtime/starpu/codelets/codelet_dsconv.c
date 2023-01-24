/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 **/

/**
 *
 * @file codelet_dsconv.c
 *
 *  CHAM codelets kernel
 *
 * @version 1.2.0
 * @author Sameh Abdulah
 *
 * @date 2019-01-20
 * @precisions mixed zc -> ds
 *
 **/
#include "chameleon_starpu.h"
#include "../include/starpu_exageostat.h"

/**
 *
 * @ingroup CORE_ChamRealDouble
 *
 **/

static void cl_dsconv_cpu_func(void *descr[], void *cl_arg) {
    int m;
    int n;
    double* A;
    int lda;
    float *B;
    int ldb;
    int *info;

    A = (double* ) STARPU_MATRIX_GET_PTR(descr[0]);
    B = (float *) STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &lda, &ldb);
    core_dlag2s(m, n, A, lda, B, ldb);
}

#if defined(EXAGEOSTAT_USE_CUDA)
static void cl_dsconv_cuda_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    double* A;
    int lda;
    float  *B;
    int ldb;
    int *info;
    A = (double* )STARPU_MATRIX_GET_PTR(descr[0]);
    B = (float *)STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &lda, &ldb);
    //RUNTIME_getStream(stream);
    cudaStream_t stream = starpu_cuda_get_local_stream();

    double2float_array( m, n, A, lda, B, ldb, CUBLAS_OP_T, stream);
    cudaStreamSynchronize( stream );
}
#endif /* defined(EXAGEOSTAT_USE_CUDA) */

static struct starpu_codelet cl_dsconv =
        {
                .where          =  STARPU_CPU | STARPU_CUDA,
                .cpu_func      =  cl_dsconv_cpu_func,
#if defined(EXAGEOSTAT_USE_CUDA)
                .cuda_func      = cl_dsconv_cuda_func,
#endif
                .nbuffers       = 2,
                .modes          = {STARPU_R, STARPU_W},
                .name           = "dsconv"
        };

void EXAGEOSTAT_TASK_dsconv(const RUNTIME_option_t *options,
                           int m, int n, int nb,
                           const CHAM_desc_t *A, int Am, int An, int lda,
                           const CHAM_desc_t *B, int Bm, int Bn, int ldb) {
    (void) nb;
    struct starpu_codelet *codelet = &cl_dsconv;

    int execution_rankA = A->get_rankof(A, Am, An);
    int execution_rankB = A->get_rankof(A, Am, An);

    starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE, &m, sizeof(int),
            STARPU_VALUE, &n, sizeof(int),
            STARPU_R, EXAGEOSTAT_RTBLKADDR(A, ChamRealDouble, Am, An),
            STARPU_VALUE, &lda, sizeof(int),
            STARPU_W, EXAGEOSTAT_RTBLKADDR(B, ChamRealFloat, Bm, Bn),
            STARPU_VALUE, &ldb, sizeof(int),
            STARPU_PRIORITY, options->priority,

#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "dsconv",
#endif
            0);
}