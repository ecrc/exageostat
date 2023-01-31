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
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon common header file
 *
 * @version 1.2.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2022-11-09
 *
 */
/**
 *  CHAMELEON facilities of interest to both CHAMELEON core developer
 *  and also of interest to CHAMELEON community contributor.
 */

#ifndef _CHAMELEON_COMMON_H_
#define _CHAMELEON_COMMON_H_

#define _GNU_SOURCE 1

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

/**
 *  Chameleon header files
 */
#include "chameleon.h"


#if defined( _WIN32 ) || defined( _WIN64 )
#include <io.h>
#else

#include <unistd.h>

#endif

/**
 * Implementation headers
 */
#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <cublas_v2.h>
#endif


#if defined(CHAMELEON_USE_OPENCL) && !defined(CHAMELEON_SIMULATION)
#include <OpenCL/cl.h>
#endif

#if defined(CHAMELEON_USE_MPI)
#include <mpi.h>
#endif

/**
 *  Line to avoid conflict with other linear algebra libraries, because, we
 *  don't know why but lapacke provide a wrong interface of lapack in fortran
 */
#ifndef LAPACK_NAME
#define LAPACK_NAME(a, b) lapackef77_##a
#endif

/** ****************************************************************************
 *  Chameleon header files
 **/
#include "global.h"
#include "auxiliary.h"
#include "context.h"
#include "descriptor.h"
#include "async.h"

/**
 *  Global shortcuts
 */
#define CHAMELEON_RANK        chameleon_rank(chamctxt)
#define CHAMELEON_NB          chamctxt->nb
#define CHAMELEON_IB          chamctxt->ib
#define CHAMELEON_RHBLK       chamctxt->rhblock
#define CHAMELEON_TRANSLATION chamctxt->translation
#define CHAMELEON_PARALLEL    chamctxt->parallel_enabled
#define CHAMELEON_STATISTICS  chamctxt->statistics_enabled

/**
 *  IPT internal define
 */
#define ChamIPT_NoDep   0
#define ChamIPT_Panel   1
#define ChamIPT_All     2

/**
 *  Global array of LAPACK constants
 */
extern char *chameleon_lapack_constants[];
#define chameleon_lapack_const(chameleon_const) chameleon_lapack_constants[chameleon_const][0]

#ifdef __cplusplus
extern "C" {
#endif

void chameleon_pmap(cham_uplo_t uplo, CHAM_desc_t *A,
                    cham_unary_operator_t operator, void *op_args,
                    RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

#if defined(__GNUC__)

static inline int chameleon_asprintf(char **strp, const char *fmt, ...) __attribute__((format(printf, 2, 3)));

#endif

static inline int chameleon_asprintf(char **strp, const char *fmt, ...) {
    va_list ap;
    int rc;

    va_start(ap, fmt);
    rc = vsprintf(strp, fmt, ap);
    va_end(ap);

    assert(rc != -1);
    return rc;
}

#ifdef __cplusplus
}
#endif

#endif
