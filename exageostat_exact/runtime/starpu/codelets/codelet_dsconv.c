/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 **/

/**
 *
 * @file codelet_dsconv.c
 *
 *  MORSE codelets kernel
 *
 * @version 1.1.0
 * @author Sameh Abdulah
 *
 * @date 2019-01-20
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

static void cl_dsconv_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    double *A;
    int lda;
    float *B;
    int ldb;
    int *info;

    A = (double *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (float *)STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &lda, &ldb);
    //printf("CPU_CORE_dlag2s========================\n");
    core_dlag2s( m, n, A, lda, B, ldb );
}


#if defined(EXAGEOSTAT_USE_CUDA)
static void cl_dsconv_cuda_func(void *descr[], void *cl_arg)
{

    int m;
    int n;
    double *A;
    int lda;
    float  *B;
    int ldb;
    int *info;
    A = (double *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (float *)STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &lda, &ldb);
    //RUNTIME_getStream(stream);
    cudaStream_t stream = starpu_cuda_get_local_stream();

    double2float_array( m, n, A, lda, B, ldb, CUBLAS_OP_T, stream);

    //#ifndef STARPU_CUDA_ASYNC
    cudaStreamSynchronize( stream );
    //#endif
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



void MORSE_TASK_dsconv(const MORSE_option_t *options,
        int m, int n, int nb,
        const MORSE_desc_t *A, int Am, int An, int lda,
        const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_dsconv;
    //void (*callback)(void*) = options->profiling ? cl_dsconv_callback : NULL;


    int execution_rankA = A->get_rankof( A, Am, An );
    int execution_rankB = A->get_rankof( A, Am, An );

    //printf ("%d -%d\n", execution_rankA, execution_rankB );

    //MORSE_BEGIN_ACCESS_DECLARATION;
    //MORSE_ACCESS_R(A, Am, An);
    //MORSE_ACCESS_W(B, Bm, Bn);
    //MORSE_END_ACCESS_DECLARATION;

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
            //#if defined(CHAMELEON_USE_MPI)
            //			STARPU_EXECUTE_ON_NODE, execution_rank,
            //#endif
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "dsconv",
#endif
            0);
}


/*
 * Codelet definition
 */
//CODELETS(dsconv, 4, cl_dsconv_cuda_func, cl_dsconv_cpu_func, STARPU_CUDA_ASYNC)
