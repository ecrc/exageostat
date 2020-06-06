/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2015 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file common.h
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2018-11-11
 *
 **/

/*******************************************************************************
 *  MORSE facilities of interest to both MORSE core developer
 *  and also of interest to MORSE community contributor.
 **/
#ifndef _MORSE_COMMON_H_
#define _MORSE_COMMON_H_


#if defined( _WIN32 ) || defined( _WIN64 )
#include <io.h>
#else
#include <unistd.h>
#endif

/** ****************************************************************************
 * Implementation headers
 **/
#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#if defined(CHAMELEON_USE_CUBLAS_V2)
#include <cublas.h>
#include <cublas_v2.h>
#else
#include <cublas.h>
#endif
#endif

#if defined(CHAMELEON_USE_OPENCL) && !defined(CHAMELEON_SIMULATION)
#include <OpenCL/cl.h>
#endif

#if defined(CHAMELEON_USE_MPI)
#include <mpi.h>
#endif

/** ****************************************************************************
 *  Line to avoid conflict with other linear algebra libraries, because, we
 *  don't know why but lapacke provide a wrong interface of lapack in fortran
 **/
#ifndef LAPACK_NAME
#define LAPACK_NAME(a, b) lapackef77_##a
#endif

/** ****************************************************************************
 *  Chameleon header files
 **/
#include "morse.h"

#include "global.h"
#include "auxiliary.h"
#include "context.h"
#include "descriptor.h"
#include "tile.h"
#include "async.h"

/*******************************************************************************
 *  Global shortcuts
 **/
#define MORSE_RANK        morse_rank(morse)
#define MORSE_SIZE        morse->world_size
#define MORSE_GRPSIZE     morse->group_size
#define MORSE_NB          morse->nb
#define MORSE_IB          morse->ib
#define MORSE_NBNBSIZE    morse->nbnbsize
#define MORSE_IBNBSIZE    morse->ibnbsize
#define MORSE_SCHEDULING  morse->scheduling
#define MORSE_RHBLK       morse->rhblock
#define MORSE_TRANSLATION morse->translation
#define MORSE_PARALLEL    morse->parallel_enabled
#define MORSE_PROFILING   morse->profiling_enabled
#if defined(CHAMELEON_USE_MPI)
#define MORSE_MPI_RANK    morse->my_mpi_rank
#define MORSE_MPI_SIZE    morse->mpi_comm_size
#endif

/*******************************************************************************
 *  IPT internal define
 **/
#define MorseIPT_NoDep   0
#define MorseIPT_Panel   1
#define MorseIPT_All     2


/*******************************************************************************
 *  Global array of LAPACK constants
 **/
extern char *morse_lapack_constants[];
#define morse_lapack_const(morse_const) morse_lapack_constants[morse_const][0]

#ifdef __cplusplus
extern "C" {
#endif

#include "compute_s.h"
#include "compute_d.h"
//#include "compute_ds.h"
#define COMPLEX
#include "compute_c.h"
#include "compute_z.h"
#undef COMPLEX

/*
void morse_pdlag2s(MORSE_context_t *morse);
void morse_pzlag2c(MORSE_context_t *morse);
void morse_pslag2d(MORSE_context_t *morse);
void morse_pclag2z(MORSE_context_t *morse);
*/

#ifdef __cplusplus
}
#endif

#endif
