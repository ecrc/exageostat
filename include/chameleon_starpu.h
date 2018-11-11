/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/
/**
 *
 * @file morse_starpu.h
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2018-11-11
 *
 **/
#ifndef _MORSE_STARPU_H_
#define _MORSE_STARPU_H_

#include "chameleon/chameleon_config.h"

/* StarPU options */
/* #undef HAVE_STARPU_FXT_PROFILING */
#define HAVE_STARPU_IDLE_PREFETCH
#define HAVE_STARPU_ITERATION_PUSH
#define HAVE_STARPU_DATA_WONT_USE
#define HAVE_STARPU_DATA_SET_COORDINATES
#define HAVE_STARPU_MALLOC_ON_NODE_SET_DEFAULT_FLAGS
/* #undef HAVE_STARPU_MPI_DATA_REGISTER */
/* #undef HAVE_STARPU_MPI_COMM_RANK */
/* #undef HAVE_STARPU_MPI_CACHED_RECEIVE */

#if defined(CHAMELEON_USE_MPI)
#include <starpu_mpi.h>
#else
#include <starpu.h>
#endif

#include <starpu_profiling.h>

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
#include <starpu_scheduler.h>
#include <starpu_cuda.h>

#include <cublas.h>
#include <starpu_cublas.h>
#if defined(CHAMELEON_USE_CUBLAS_V2)
#include <cublas_v2.h>
#include <starpu_cublas_v2.h>
#endif
#endif

#include "common.h"
#include "runtime/starpu/runtime_codelets.h"
#include "runtime/starpu/runtime_profiling.h"
#include "runtime/starpu/runtime_codelet_profile.h"
#include "runtime/starpu/runtime_workspace.h"

typedef struct starpu_conf starpu_conf_t;

/******************************************************************************/

/*
 * MPI Redefinitions
 */
#if defined(CHAMELEON_USE_MPI)
#undef STARPU_REDUX
//#define starpu_insert_task(...) starpu_mpi_insert_task(MPI_COMM_WORLD, __VA_ARGS__)
#define starpu_insert_task starpu_mpi_insert_task
#define starpu_mpi_codelet(_codelet_) MPI_COMM_WORLD, _codelet_

#else

#define starpu_mpi_codelet(_codelet_) _codelet_

#endif

/*
 * cuBlasAPI v2 - StarPU enable the support for cublas handle
 */
#if defined(CHAMELEON_USE_CUDA) && defined(CHAMELEON_USE_CUBLAS_V2)
#define RUNTIME_getStream(_stream_)                             \
    cublasHandle_t _stream_ = starpu_cublas_get_local_handle();
#else
#define RUNTIME_getStream(_stream_)                             \
    cudaStream_t _stream_ = starpu_cuda_get_local_stream();     \
    cublasSetKernelStream( stream );

#endif

/*
 * Enable codelets names
 */
#if (STARPU_MAJOR_VERSION > 1) || ((STARPU_MAJOR_VERSION == 1) && (STARPU_MINOR_VERSION > 1))
#define CHAMELEON_CODELETS_HAVE_NAME
#endif

/**
 * Access to block pointer and leading dimension
 */
#define RTBLKADDR( desc, type, m, n ) ( (starpu_data_handle_t)RUNTIME_data_getaddr( desc, m, n ) )

void RUNTIME_set_reduction_methods(starpu_data_handle_t handle, MORSE_enum dtyp);
#ifdef CHAMELEON_USE_MPI
#ifdef HAVE_STARPU_MPI_CACHED_RECEIVE
int RUNTIME_desc_iscached(const MORSE_desc_t *A, int Am, int An);
#endif
#endif

#if defined(CHAMELEON_USE_MPI)
#  if defined(HAVE_STARPU_MPI_CACHED_RECEIVE)
#    define RUNTIME_ACCESS_WRITE_CACHED(A, Am, An) do { if (RUNTIME_desc_iscached(A, Am, An)) __morse_need_submit = 1; } while(0)
#  else
#    warning "WAR dependencies need starpu_mpi_cached_receive support from StarPU 1.2.1 or greater"
#    define RUNTIME_ACCESS_WRITE_CACHED(A, Am, An)
#  endif
#else
#define RUNTIME_ACCESS_WRITE_CACHED(A, Am, An)
#endif

#ifdef CHAMELEON_ENABLE_PRUNING_STATS

#define RUNTIME_PRUNING_STATS_BEGIN_ACCESS_DECLARATION \
    int __morse_exec = 0; \
    int __morse_changed = 0;

#define RUNTIME_PRUNING_STATS_ACCESS_W(A, Am, An) \
    if (morse_desc_islocal(A, Am, An)) \
        __morse_exec = 1;

#define RUNTIME_PRUNING_STATS_END_ACCESS_DECLARATION \
    RUNTIME_total_tasks++; \
    if (__morse_exec) \
        RUNTIME_exec_tasks++; \
    else if (__morse_need_submit) \
        RUNTIME_comm_tasks++; \
    else if (__morse_changed) \
        RUNTIME_changed_tasks++;

#define RUNTIME_PRUNING_STATS_RANK_CHANGED(rank) \
    int __morse_myrank; \
    RUNTIME_comm_rank(&__morse_myrank); \
    __morse_exec = (rank) == __morse_myrank; \
    __morse_changed = 1; \

#else
#define RUNTIME_PRUNING_STATS_BEGIN_ACCESS_DECLARATION
#define RUNTIME_PRUNING_STATS_ACCESS_W(A, Am, An)
#define RUNTIME_PRUNING_STATS_END_ACCESS_DECLARATION
#define RUNTIME_PRUNING_STATS_RANK_CHANGED(rank)
#endif

#define RUNTIME_BEGIN_ACCESS_DECLARATION        \
    RUNTIME_PRUNING_STATS_BEGIN_ACCESS_DECLARATION

#define RUNTIME_ACCESS_R(A, Am, An)

#define RUNTIME_ACCESS_W(A, Am, An)             \
    RUNTIME_PRUNING_STATS_ACCESS_W(A, Am, An);  \
    RUNTIME_ACCESS_WRITE_CACHED(A, Am, An)

#define RUNTIME_ACCESS_RW(A, Am, An)            \
    RUNTIME_PRUNING_STATS_ACCESS_W(A, Am, An);  \
    RUNTIME_ACCESS_WRITE_CACHED(A, Am, An)

#define RUNTIME_RANK_CHANGED(rank)              \
    RUNTIME_PRUNING_STATS_RANK_CHANGED(rank)

#define RUNTIME_END_ACCESS_DECLARATION          \
    RUNTIME_PRUNING_STATS_END_ACCESS_DECLARATION;

#endif /* _MORSE_STARPU_H_ */
