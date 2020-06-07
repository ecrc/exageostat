/**
 *
 * @file runtime_codelet_s.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU float codelets header
 *
 * @version 1.0.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2011-06-01
 * @generated s Sun Nov 25 16:59:42 2018
 *
 */
#ifndef _CODELETS_S_H_
#define _CODELETS_S_H_

#include <stdio.h>
#include "runtime/starpu/runtime_codelets.h"

#include "chameleon/morse_tasks_d.h"
#include "chameleon/morse_tasks_ds.h"
#if !defined(CHAMELEON_SIMULATION)
#include "coreblas/coreblas_d.h"
#if defined(CHAMELEON_USE_CUDA)
#include "cudablas.h"
#endif
#endif
/*
 * Management functions
 */
SCODELETS_HEADER(tile_zero)

/*
 * BLAS 1 functions
 */
SCODELETS_HEADER(axpy)

/*
 * BLAS 3 functions
 */
SCODELETS_HEADER(gemm)
SCODELETS_HEADER(symm)
SCODELETS_HEADER(syr2k)
SCODELETS_HEADER(syrk)
SCODELETS_HEADER(symm)
SCODELETS_HEADER(syr2k)
SCODELETS_HEADER(syrk)
SCODELETS_HEADER(trmm)
SCODELETS_HEADER(trsm)

/*
 * LAPACK functions
 */
SCODELETS_HEADER(gelqt)
SCODELETS_HEADER(geqrt)
SCODELETS_HEADER(gessm)
SCODELETS_HEADER(gessq)
SCODELETS_HEADER(getrf)
SCODELETS_HEADER(getrf_incpiv)
SCODELETS_HEADER(getrf_nopiv)
SCODELETS_HEADER(syrfb)
SCODELETS_HEADER(lauum)
SCODELETS_HEADER(potrf)
SCODELETS_HEADER(ssssm)
SCODELETS_HEADER(syssq)
SCODELETS_HEADER(trasm)
SCODELETS_HEADER(trssq)
SCODELETS_HEADER(trtri)
SCODELETS_HEADER(tplqt)
SCODELETS_HEADER(tpqrt)
SCODELETS_HEADER(tpmlqt)
SCODELETS_HEADER(tpmqrt)
SCODELETS_HEADER(tslqt)
SCODELETS_HEADER(tsmlq)
SCODELETS_HEADER(tsmqr)
SCODELETS_HEADER(tsmlq_hetra1)
SCODELETS_HEADER(tsmqr_hetra1)
SCODELETS_HEADER(tsqrt)
SCODELETS_HEADER(tstrf)
SCODELETS_HEADER(ttlqt)
SCODELETS_HEADER(ttmlq)
SCODELETS_HEADER(ttmqr)
SCODELETS_HEADER(ttqrt)
SCODELETS_HEADER(ormlq)
SCODELETS_HEADER(ormqr)

/*
 * Auxiliary functions
 */
SCODELETS_HEADER(geadd)
SCODELETS_HEADER(he2ge)
SCODELETS_HEADER(lascal)
SCODELETS_HEADER(tradd)
SCODELETS_HEADER(lacpy)
SCODELETS_HEADER(lange)
SCODELETS_HEADER(lange_max)
SCODELETS_HEADER(lansy)
SCODELETS_HEADER(lantr)
SCODELETS_HEADER(laset)
SCODELETS_HEADER(laset2)
SCODELETS_HEADER(latro)
SCODELETS_HEADER(plssq)
SCODELETS_HEADER(plssq2)

/*
 * MIXED PRECISION functions
 */
SCODELETS_HEADER(lag2d)

/*
 * DZ functions
 */
SCODELETS_HEADER(asum)

/*
 * CPU only functions
 */
SCODELETS_HEADER(plrnt)
SCODELETS_HEADER(build)

#if defined(PRECISION_z) || defined(PRECISION_c)
SCODELETS_HEADER(syssq)
SCODELETS_HEADER(lansy)
SCODELETS_HEADER(plgsy)
SCODELETS_HEADER(sytrf_nopiv)
#endif
SCODELETS_HEADER(plgsy)

#endif /* _CODELETS_S_H_ */
