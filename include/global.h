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
 * @file global.h
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Piotr Luszczek
 * @author Cedric Castagnede
 * @date 2018-11-11
 *
 **/

/*******************************************************************************
 *  MORSE internals of interest to MORSE core developers, but not necessarily
 *  of interest to MORSE community contributors.
 **/
#ifndef _MORSE_GLOBAL_H_
#define _MORSE_GLOBAL_H_

#if defined( _WIN32 ) || defined( _WIN64 )
#include "control/morsewinthread.h"
#else
#include <pthread.h>
#endif

/*******************************************************************************
 *  Numerical operations
 **/
#define MORSE_FUNC_SGELS    1
#define MORSE_FUNC_SPOSV    2
#define MORSE_FUNC_SGESV    3
#define MORSE_FUNC_DGELS    4
#define MORSE_FUNC_DPOSV    5
#define MORSE_FUNC_DGESV    6
#define MORSE_FUNC_CGELS    7
#define MORSE_FUNC_CPOSV    8
#define MORSE_FUNC_CGESV    9
#define MORSE_FUNC_ZGELS   10
#define MORSE_FUNC_ZPOSV   11
#define MORSE_FUNC_ZGESV   12
#define MORSE_FUNC_ZCGESV  13
#define MORSE_FUNC_DSGESV  14
#define MORSE_FUNC_ZCPOSV  15
#define MORSE_FUNC_DSPOSV  16
#define MORSE_FUNC_DSGELS  17
#define MORSE_FUNC_ZCGELS  18
#define MORSE_FUNC_SGEMM   19
#define MORSE_FUNC_DGEMM   20
#define MORSE_FUNC_CGEMM   21
#define MORSE_FUNC_ZGEMM   22
#define MORSE_FUNC_SSYMM   23
#define MORSE_FUNC_DSYMM   24
#define MORSE_FUNC_CSYMM   25
#define MORSE_FUNC_ZSYMM   26
#define MORSE_FUNC_CHERK   27
#define MORSE_FUNC_ZHERK   28
#define MORSE_FUNC_SSYRK   29
#define MORSE_FUNC_DSYRK   30
#define MORSE_FUNC_CSYRK   31
#define MORSE_FUNC_ZSYRK   32
#define MORSE_FUNC_CHEMM   33
#define MORSE_FUNC_ZHEMM   34
#define MORSE_FUNC_ZHEEV   35
#define MORSE_FUNC_CHEEV   36
#define MORSE_FUNC_DSYEV   37
#define MORSE_FUNC_SSYEV   38
#define MORSE_FUNC_ZHEEVD  39
#define MORSE_FUNC_CHEEVD  40
#define MORSE_FUNC_DSYEVD  41
#define MORSE_FUNC_SSYEVD  42
#define MORSE_FUNC_ZHEGST  43
#define MORSE_FUNC_CHEGST  44
#define MORSE_FUNC_DSYGST  45
#define MORSE_FUNC_SSYGST  46
#define MORSE_FUNC_ZHEGV   47
#define MORSE_FUNC_CHEGV   48
#define MORSE_FUNC_DSYGV   49
#define MORSE_FUNC_SSYGV   50
#define MORSE_FUNC_ZHEGVD  51
#define MORSE_FUNC_CHEGVD  52
#define MORSE_FUNC_DSYGVD  53
#define MORSE_FUNC_SSYGVD  54
#define MORSE_FUNC_ZHETRD  55
#define MORSE_FUNC_CHETRD  56
#define MORSE_FUNC_DSYTRD  57
#define MORSE_FUNC_SSYTRD  58
#define MORSE_FUNC_ZGESVD  59
#define MORSE_FUNC_CGESVD  60
#define MORSE_FUNC_DGESVD  61
#define MORSE_FUNC_SGESVD  62
#define MORSE_FUNC_ZGEEV   63
#define MORSE_FUNC_CGEEV   64
#define MORSE_FUNC_DGEEV   65
#define MORSE_FUNC_SGEEV   66
#define MORSE_FUNC_ZGEHRD  67
#define MORSE_FUNC_CGEHRD  68
#define MORSE_FUNC_DGEHRD  69
#define MORSE_FUNC_SGEHRD  70
#define MORSE_FUNC_ZGEBRD  71
#define MORSE_FUNC_CGEBRD  72
#define MORSE_FUNC_DGEBRD  73
#define MORSE_FUNC_SGEBRD  74
#define MORSE_FUNC_ZSYSV  75
#define MORSE_FUNC_CSYSV  76

#endif
