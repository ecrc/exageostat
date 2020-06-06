/**
 *
 * @copyright (c) 2017-2020 King Abdullah University of
 *                          Sience and technelogy (KAUST).
 *                          All rights reserved.
 *
 **/

/**
 *
 * @file codelet_sdconv.c
 *
 *  MORSE codelets kernel
 *
 * @version 1.1.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Sameh Abdulah
 * @date 2019-01-27
 * @precisions mixed zc -> ds
 *
 **/
//#include "../../../../include/coreblas_ds.h"
#include "chameleon_starpu.h"
#include "runtime_codelet_d.h"
#include "../include/starpu_exageostat.h"
/**
 *
 * @ingroup CORE_MorseRealDouble
 *
 **/

static void cl_sdconv_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    float *A;
    int lda;
    double *B;
    int ldb;

    A = (float *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (double *)STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &lda, &ldb);

    //printf("CPUCORE_slag2d**********\n");
    core_slag2d( m, n, A, lda, B, ldb);
}

#if defined(EXAGEOSTAT_USE_CUDA)
static void cl_sdconv_cuda_func(void *descr[], void *cl_arg)
{

    int m;
    int n;
    float *A;
    int lda;
    double *B;
    int ldb;
    int *info;

    A = (float *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (double *)STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &lda, &ldb);
    //      RUNTIME_getStream(stream);
    cudaStream_t stream = starpu_cuda_get_local_stream();
    float2double_array( m, n, A, lda, B, ldb, CUBLAS_OP_T, stream);

    //#ifndef STARPU_CUDA_ASYNC
    cudaStreamSynchronize( stream );
    //#endif

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
void MORSE_TASK_sdconv(const MORSE_option_t *options,
        int m, int n, int nb,
        const MORSE_desc_t *A, int Am, int An, int lda,
        const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_sdconv;
    //void (*callback)(void*) = options->profiling ? cl_sdconv_callback : NULL;
    //	int execution_rank = B->get_rankof( B, Bm, Bn );

    //	MORSE_BEGIN_ACCESS_DECLARATION;
    //	MORSE_ACCESS_R(A, Am, An);
    //	MORSE_ACCESS_W(B, Bm, Bn);
    //	MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE,    &m,                 sizeof(int),
            STARPU_VALUE,    &n,                 sizeof(int),
            STARPU_R,        EXAGEOSTAT_RTBLKADDR(A, MorseRealFloat, Am, An),
            STARPU_VALUE,    &lda,               sizeof(int),
            STARPU_W,        EXAGEOSTAT_RTBLKADDR(B, MorseRealDouble, Bm, Bn),
            STARPU_VALUE,    &ldb,               sizeof(int),
            STARPU_PRIORITY,  options->priority,
            //   STARPU_CALLBACK,  callback,
            //#if defined(CHAMELEON_USE_MPI)
            //			STARPU_EXECUTE_ON_NODE, execution_rank,
            //#endif
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "sdconv",
#endif
            0);




}


/*
   void MORSE_TASK_sdconv(const MORSE_option_t *options,
   int m, int n, int nb,
   const MORSE_desc_t *A, int Am, int An, int lda,
   const MORSE_desc_t *B, int Bm, int Bn, int ldb)
   {
   (void)nb;
   struct starpu_codelet *codelet = &cl_sdconv;
   void (*callback)(void*) = options->profiling ? cl_sdconv_callback : NULL;

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
STARPU_NAME, "sdconv",
#endif
0);
}
}
*/


/*
 * Codelet definition
 */
//CODELETS_CPU(sdconv, 1, cl_sdconv_cpu_func)
