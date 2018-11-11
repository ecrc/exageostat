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
 * @file control/descriptor.h
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2018-11-11
 *
 **/
#ifndef _MORSE_DESCRIPTOR_H_
#define _MORSE_DESCRIPTOR_H_

#include <assert.h>
#include "chameleon/chameleon_config.h"
#include "chameleon/morse_struct.h"
#include "auxiliary.h"
#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************************
 *  Internal routines
 **/
inline static void* morse_geteltaddr(const MORSE_desc_t *A, int m, int n, int eltsize);
inline static void* morse_getaddr_cm    (const MORSE_desc_t *A, int m, int n);
inline static void* morse_getaddr_ccrb  (const MORSE_desc_t *A, int m, int n);
inline static void* morse_getaddr_null  (const MORSE_desc_t *A, int m, int n);
inline static int   morse_getblkldd_cm  (const MORSE_desc_t *A, int m);
inline static int   morse_getblkldd_ccrb(const MORSE_desc_t *A, int m);

/*****************************************************************
 *  Data distributions
 */
inline static int   morse_getrankof_2d(const MORSE_desc_t *desc, int m, int n);
inline static int   morse_getrankof_2d_diag(const MORSE_desc_t *desc, int m, int n);

MORSE_desc_t morse_desc_init(MORSE_enum dtyp, int mb, int nb, int bsiz,
                             int lm, int ln, int i, int j, int m, int n, int p, int q);
MORSE_desc_t morse_desc_init_diag(MORSE_enum dtyp, int mb, int nb, int bsiz,
                                  int lm, int ln, int i, int j, int m, int n, int p, int q);
MORSE_desc_t morse_desc_init_user(MORSE_enum dtyp, int mb, int nb, int bsiz,
                                  int lm, int ln, int i, int j,
                                  int m,  int n,  int p, int q,
                                  void* (*get_blkaddr)( const MORSE_desc_t*, int, int ),
                                  int (*get_blkldd)( const MORSE_desc_t*, int ),
                                  int (*get_rankof)( const MORSE_desc_t*, int, int ));
MORSE_desc_t* morse_desc_submatrix(MORSE_desc_t *descA, int i, int j, int m, int n);

int morse_desc_check    (MORSE_desc_t *desc);
int morse_desc_mat_alloc(MORSE_desc_t *desc);
int morse_desc_mat_free (MORSE_desc_t *desc);

#define BLKLDD(A, k) A->get_blkldd( A,k )

/*******************************************************************************
 *  Internal function to return address of block (m,n) with m,n = block indices
 **/
inline static void* morse_getaddr_ccrb(const MORSE_desc_t *A, int m, int n)
{
    size_t mm = m + A->i / A->mb;
    size_t nn = n + A->j / A->nb;
    size_t eltsize = MORSE_Element_Size(A->dtyp);
    size_t offset = 0;

#if defined(CHAMELEON_USE_MPI)
    assert( A->myrank == A->get_rankof( A, mm, nn) );
    mm = mm / A->p;
    nn = nn / A->q;
#endif

    if (mm < (size_t)(A->llm1)) {
        if (nn < (size_t)(A->lln1))
            offset = (size_t)(A->bsiz) * (mm + (size_t)(A->llm1) * nn);
        else
            offset = A->A12 + ((size_t)(A->mb * (A->lln%A->nb)) * mm);
    }
    else {
        if (nn < (size_t)(A->lln1))
            offset = A->A21 + ((size_t)((A->llm%A->mb) * A->nb) * nn);
        else
            offset = A->A22;
    }

    return (void*)((intptr_t)A->mat + (offset*eltsize) );
}

/*******************************************************************************
 *  Internal function to return address of block (m,n) with m,n = block indices
 **/
inline static void *morse_getaddr_cm(const MORSE_desc_t *A, int m, int n)
{
    size_t mm = m + A->i / A->mb;
    size_t nn = n + A->j / A->nb;
    size_t eltsize = MORSE_Element_Size(A->dtyp);
    size_t offset = 0;

#if defined(CHAMELEON_USE_MPI)
    assert( A->myrank == A->get_rankof( A, mm, nn) );
    mm = mm / A->p;
    nn = nn / A->q;
#endif

    offset = (size_t)(A->llm * A->nb) * nn + (size_t)(A->mb) * mm;
    return (void*)((intptr_t)A->mat + (offset*eltsize) );
}

/*******************************************************************************
 *  Internal function to return address of block (m,n) with m,n = block indices
 *  This version lets the runtime allocate on-demand.
 **/
inline static void *morse_getaddr_null(const MORSE_desc_t *A, int m, int n)
{
    (void)A; (void)m; (void)n;
    return NULL;
}

/*******************************************************************************
 *  Internal function to return address of element A(m,n) with m,n = matrix indices
 **/
inline static void* morse_geteltaddr(const MORSE_desc_t *A, int m, int n, int eltsize) // Not used anywhere ?!
{
    size_t mm = (m + A->i)/A->mb;
    size_t nn = (n + A->j)/A->nb;
    size_t offset = 0;

#if defined(CHAMELEON_USE_MPI)
    assert( A->myrank == A->get_rankof( A, mm, nn) );
    mm = mm / A->p;
    nn = nn / A->q;
#endif

    if (mm < (size_t)(A->llm1)) {
        if (nn < (size_t)(A->lln1))
            offset = A->bsiz*(mm+A->llm1*nn) + m%A->mb + A->mb*(n%A->nb);
        else
            offset = A->A12 + (A->mb*(A->lln%A->nb)*mm) + m%A->mb + A->mb*(n%A->nb);
    }
    else {
        if (nn < (size_t)(A->lln1))
            offset = A->A21 + ((A->llm%A->mb)*A->nb*nn) + m%A->mb + (A->llm%A->mb)*(n%A->nb);
        else
            offset = A->A22 + m%A->mb  + (A->llm%A->mb)*(n%A->nb);
    }
    return (void*)((intptr_t)A->mat + (offset*eltsize) );
}

/*******************************************************************************
 *  Internal function to return the leading dimension of element A(m,*) with m,n = block indices
 **/
inline static int morse_getblkldd_ccrb(const MORSE_desc_t *A, int m)
{
    int mm = m + A->i / A->mb;
    return ( ((mm+1) == A->lmt) && ((A->lm % A->mb) != 0)) ? A->lm % A->mb : A->mb;
}

inline static int morse_getblkldd_cm(const MORSE_desc_t *A, int m) {
    (void)m;
    return A->llm;
}


/*******************************************************************************
 *  Internal function to return MPI rank of element A(m,n) with m,n = block indices
 **/
inline static int morse_getrankof_2d(const MORSE_desc_t *desc, int m, int n)
{
    return (m % desc->p) * desc->q + (n % desc->q);
}

/*******************************************************************************
 *  Internal function to return MPI rank of element DIAG(m,0) with m,n = block indices
 **/
inline static int morse_getrankof_2d_diag(const MORSE_desc_t *desc, int m, int n)
{
    assert( n == 0 );
    return (m % desc->p) * desc->q + (m % desc->q);
}


/*******************************************************************************
 * Detect if the tile is local or not
 **/
inline static int morse_desc_islocal( const MORSE_desc_t *A, int m, int n )
{
#if defined(CHAMELEON_USE_MPI)
    return (A->myrank == A->get_rankof(A, m, n));
#else
    (void)A; (void)m; (void)n;
    return 1;
#endif /* defined(CHAMELEON_USE_MPI) */
}

/*******************************************************************************
 * Declare data accesses of codelets using these macros, for instance:
 * MORSE_BEGIN_ACCESS_DECLARATION
 * MORSE_ACCESS_R(A, Am, An)
 * MORSE_ACCESS_R(B, Bm, Bn)
 * MORSE_ACCESS_RW(C, Cm, Cn)
 * MORSE_END_ACCESS_DECLARATION
 */
#define MORSE_BEGIN_ACCESS_DECLARATION { \
    unsigned __morse_need_submit = 0; \
    RUNTIME_BEGIN_ACCESS_DECLARATION

#define MORSE_ACCESS_R(A, Am, An) do { \
    if (morse_desc_islocal(A, Am, An)) __morse_need_submit = 1; \
    RUNTIME_ACCESS_R(A, Am, An); \
} while(0)

#define MORSE_ACCESS_W(A, Am, An) do { \
    if (morse_desc_islocal(A, Am, An)) __morse_need_submit = 1; \
    RUNTIME_ACCESS_W(A, Am, An); \
} while(0)

#define MORSE_ACCESS_RW(A, Am, An) do { \
    if (morse_desc_islocal(A, Am, An)) __morse_need_submit = 1; \
    RUNTIME_ACCESS_RW(A, Am, An); \
} while(0)

#define MORSE_RANK_CHANGED(rank) do {\
    __morse_need_submit = 1; \
    RUNTIME_RANK_CHANGED(rank); \
} while (0)

#define MORSE_END_ACCESS_DECLARATION \
    RUNTIME_END_ACCESS_DECLARATION; \
    if (!__morse_need_submit) return; \
}

#ifdef __cplusplus
}
#endif

#endif
