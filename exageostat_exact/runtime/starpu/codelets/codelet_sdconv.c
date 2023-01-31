/**
 *
 * @copyright (c) 2017-2023 King Abdullah University of
 *                          Sience and technelogy (KAUST).
 *                          All rights reserved.
 *
 **/

/**
 *
 * @file codelet_sdconv.c
 *
 *  CHAMELEON codelets kernel
 *
 * @version 1.2.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Sameh Abdulah
 * @date 2019-01-27
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

static void cl_sdconv_cpu_func(void *descr[], void *cl_arg) {
    int m;
    int n;
    float *A;
    int lda;
    double* B;
    int ldb;

    A = (float *) STARPU_MATRIX_GET_PTR(descr[0]);
    B = (double* ) STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &lda, &ldb);

    core_slag2d(m, n, A, lda, B, ldb);
}

#if defined(EXAGEOSTAT_USE_CUDA)
static void cl_sdconv_cuda_func(void *descr[], void *cl_arg)
{

    int m;
    int n;
    float *A;
    int lda;
    double* B;
    int ldb;
    int *info;

    A = (float *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (double* )STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &lda, &ldb);
    cudaStream_t stream = starpu_cuda_get_local_stream();
    float2double_array( m, n, A, lda, B, ldb, CUBLAS_OP_T, stream);
    cudaStreamSynchronize( stream );
}
#endif /* defined(EXAGEOSTAT_USE_CUDA) */

static struct starpu_codelet cl_sdconv =
        {
                .where          = STARPU_CPU | STARPU_CUDA,
                .cpu_func      = {cl_sdconv_cpu_func},
#if defined(EXAGEOSTAT_USE_CUDA)
                .cuda_func      = {cl_sdconv_cuda_func},
#endif
                .nbuffers       = 2,
                .modes          = {STARPU_R, STARPU_W},
                .name           = "sdconv"
        };

void EXAGEOSTAT_TASK_sdconv(const RUNTIME_option_t *options,
                           int m, int n, int nb,
                           const CHAM_desc_t *A, int Am, int An, int lda,
                           const CHAM_desc_t *B, int Bm, int Bn, int ldb) {
    (void) nb;
    struct starpu_codelet *codelet = &cl_sdconv;

    starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE, &m, sizeof(int),
            STARPU_VALUE, &n, sizeof(int),
            STARPU_R, EXAGEOSTAT_RTBLKADDR(A, ChamRealFloat, Am, An),
            STARPU_VALUE, &lda, sizeof(int),
            STARPU_W, EXAGEOSTAT_RTBLKADDR(B, ChamRealDouble, Bm, Bn),
            STARPU_VALUE, &ldb, sizeof(int),
            STARPU_PRIORITY, options->priority,

#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "sdconv",
#endif
            0);
}