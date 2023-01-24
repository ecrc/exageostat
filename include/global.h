/**
 *
 * @file global.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon global variables header
 *
 * @version 1.2.0
 * @author Jakub Kurzak
 * @author Piotr Luszczek
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2022-11-09
 *
 */
/**
 *  CHAMELEON internals of interest to CHAMELEON core developers, but not necessarily
 *  of interest to CHAMELEON community contributors.
 */
#ifndef _chameleon_global_h_
#define _chameleon_global_h_

/**
 *  Numerical operations
 */
#define CHAMELEON_FUNC_SGELS    1
#define CHAMELEON_FUNC_SPOSV    2
#define CHAMELEON_FUNC_SGESV    3
#define CHAMELEON_FUNC_DGELS    4
#define CHAMELEON_FUNC_DPOSV    5
#define CHAMELEON_FUNC_DGESV    6
#define CHAMELEON_FUNC_CGELS    7
#define CHAMELEON_FUNC_CPOSV    8
#define CHAMELEON_FUNC_CGESV    9
#define CHAMELEON_FUNC_ZGELS   10
#define CHAMELEON_FUNC_ZPOSV   11
#define CHAMELEON_FUNC_ZGESV   12
#define CHAMELEON_FUNC_ZCGESV  13
#define CHAMELEON_FUNC_DSGESV  14
#define CHAMELEON_FUNC_ZCPOSV  15
#define CHAMELEON_FUNC_DSPOSV  16
#define CHAMELEON_FUNC_DSGELS  17
#define CHAMELEON_FUNC_ZCGELS  18
#define CHAMELEON_FUNC_SGEMM   19
#define CHAMELEON_FUNC_DGEMM   20
#define CHAMELEON_FUNC_CGEMM   21
#define CHAMELEON_FUNC_ZGEMM   22
#define CHAMELEON_FUNC_SSYMM   23
#define CHAMELEON_FUNC_DSYMM   24
#define CHAMELEON_FUNC_CSYMM   25
#define CHAMELEON_FUNC_ZSYMM   26
#define CHAMELEON_FUNC_CHERK   27
#define CHAMELEON_FUNC_ZHERK   28
#define CHAMELEON_FUNC_SSYRK   29
#define CHAMELEON_FUNC_DSYRK   30
#define CHAMELEON_FUNC_CSYRK   31
#define CHAMELEON_FUNC_ZSYRK   32
#define CHAMELEON_FUNC_CHEMM   33
#define CHAMELEON_FUNC_ZHEMM   34
#define CHAMELEON_FUNC_ZHEEV   35
#define CHAMELEON_FUNC_CHEEV   36
#define CHAMELEON_FUNC_DSYEV   37
#define CHAMELEON_FUNC_SSYEV   38
#define CHAMELEON_FUNC_ZHEEVD  39
#define CHAMELEON_FUNC_CHEEVD  40
#define CHAMELEON_FUNC_DSYEVD  41
#define CHAMELEON_FUNC_SSYEVD  42
#define CHAMELEON_FUNC_ZHEGST  43
#define CHAMELEON_FUNC_CHEGST  44
#define CHAMELEON_FUNC_DSYGST  45
#define CHAMELEON_FUNC_SSYGST  46
#define CHAMELEON_FUNC_ZHEGV   47
#define CHAMELEON_FUNC_CHEGV   48
#define CHAMELEON_FUNC_DSYGV   49
#define CHAMELEON_FUNC_SSYGV   50
#define CHAMELEON_FUNC_ZHEGVD  51
#define CHAMELEON_FUNC_CHEGVD  52
#define CHAMELEON_FUNC_DSYGVD  53
#define CHAMELEON_FUNC_SSYGVD  54
#define CHAMELEON_FUNC_ZHETRD  55
#define CHAMELEON_FUNC_CHETRD  56
#define CHAMELEON_FUNC_DSYTRD  57
#define CHAMELEON_FUNC_SSYTRD  58
#define CHAMELEON_FUNC_ZGESVD  59
#define CHAMELEON_FUNC_CGESVD  60
#define CHAMELEON_FUNC_DGESVD  61
#define CHAMELEON_FUNC_SGESVD  62
#define CHAMELEON_FUNC_ZGEEV   63
#define CHAMELEON_FUNC_CGEEV   64
#define CHAMELEON_FUNC_DGEEV   65
#define CHAMELEON_FUNC_SGEEV   66
#define CHAMELEON_FUNC_ZGEHRD  67
#define CHAMELEON_FUNC_CGEHRD  68
#define CHAMELEON_FUNC_DGEHRD  69
#define CHAMELEON_FUNC_SGEHRD  70
#define CHAMELEON_FUNC_ZGEBRD  71
#define CHAMELEON_FUNC_CGEBRD  72
#define CHAMELEON_FUNC_DGEBRD  73
#define CHAMELEON_FUNC_SGEBRD  74
#define CHAMELEON_FUNC_ZSYSV  75
#define CHAMELEON_FUNC_CSYSV  76

typedef int CHAM_enum;

#endif /* _chameleon_global_h_ */
