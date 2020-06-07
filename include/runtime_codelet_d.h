/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2014 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file runtime_codelet_d.h
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 1.0.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2018-11-11
 * @generated d Fri Dec  1 14:38:51 2017
 *
 **/

#ifndef _CODELETS_D_H_
#define _CODELETS_D_H_

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
DCODELETS_HEADER(tile_zero)

/*
 * BLAS 1 functions
 */
DCODELETS_HEADER(axpy)

/*
 * BLAS 3 functions
 */
DCODELETS_HEADER(gemm)
DCODELETS_HEADER(symm)
DCODELETS_HEADER(syr2k)
DCODELETS_HEADER(syrk)
DCODELETS_HEADER(symm)
DCODELETS_HEADER(syr2k)
DCODELETS_HEADER(syrk)
DCODELETS_HEADER(trmm)
DCODELETS_HEADER(trsm)

/*
 * LAPACK functions
 */
DCODELETS_HEADER(gelqt)
DCODELETS_HEADER(geqrt)
DCODELETS_HEADER(gessm)
DCODELETS_HEADER(gessq)
DCODELETS_HEADER(getrf)
DCODELETS_HEADER(getrf_incpiv)
DCODELETS_HEADER(getrf_nopiv)
DCODELETS_HEADER(syrfb)
DCODELETS_HEADER(lauum)
DCODELETS_HEADER(potrf)
DCODELETS_HEADER(ssssm)
DCODELETS_HEADER(syssq)
DCODELETS_HEADER(trasm)
DCODELETS_HEADER(trssq)
DCODELETS_HEADER(trtri)
DCODELETS_HEADER(tpqrt)
DCODELETS_HEADER(tpmqrt)
DCODELETS_HEADER(tslqt)
DCODELETS_HEADER(tsmlq)
DCODELETS_HEADER(tsmqr)
DCODELETS_HEADER(tsmlq_hetra1)
DCODELETS_HEADER(tsmqr_hetra1)
DCODELETS_HEADER(tsqrt)
DCODELETS_HEADER(tstrf)
DCODELETS_HEADER(ttlqt)
DCODELETS_HEADER(ttmlq)
DCODELETS_HEADER(ttmqr)
DCODELETS_HEADER(ttqrt)
DCODELETS_HEADER(ormlq)
DCODELETS_HEADER(ormqr)

/*
 * Auxiliary functions
 */
DCODELETS_HEADER(geadd)
DCODELETS_HEADER(he2ge)
DCODELETS_HEADER(lascal)
DCODELETS_HEADER(tradd)
DCODELETS_HEADER(lacpy)
DCODELETS_HEADER(lange)
DCODELETS_HEADER(lange_max)
DCODELETS_HEADER(lansy)
DCODELETS_HEADER(lantr)
DCODELETS_HEADER(laset)
DCODELETS_HEADER(laset2)
DCODELETS_HEADER(latro)
DCODELETS_HEADER(plssq)
DCODELETS_HEADER(plssq2)

/*
 * MIXED PRECISION functions
 */
DCODELETS_HEADER(lag2s)
DCODELETS_HEADER(dsconv)
DCODELETS_HEADER(sdconv)
/*
 * DZ functions
 */
DCODELETS_HEADER(asum)

/*
 * CPU only functions
 */
DCODELETS_HEADER(plrnt)
DCODELETS_HEADER(build)

#if defined(PRECISION_z) || defined(PRECISION_c)
DCODELETS_HEADER(syssq)
DCODELETS_HEADER(lansy)
DCODELETS_HEADER(plgsy)
DCODELETS_HEADER(sytrf_nopiv)
#endif
DCODELETS_HEADER(plgsy)

#endif /* _CODELETS_D_H_ */
