/**
 *
 * @file descriptor.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon descriptor header
 *
 * @version 1.2.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Guillaume Sylvand
 * @author Raphael Boucherie
 * @author Samuel Thibault
 * @date 2022-11-09
 *
 */
#ifndef _chameleon_descriptor_h_
#define _chameleon_descriptor_h_

#include <assert.h>
#include "chameleon/include/chameleon/struct.h"
#include "auxiliary.h"
#include "chameleon/include/chameleon/constants.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 *  Internal routines
 */
inline static void *chameleon_geteltaddr(const CHAM_desc_t *A, int m, int n, int eltsize);

inline static void *chameleon_getaddr_cm(const CHAM_desc_t *A, int m, int n);

inline static void *chameleon_getaddr_ccrb(const CHAM_desc_t *A, int m, int n);

inline static void *chameleon_getaddr_null(const CHAM_desc_t *A, int m, int n);

inline static void *chameleon_getaddr_diag(const CHAM_desc_t *A, int m, int n);

inline static int chameleon_getblkldd_cm(const CHAM_desc_t *A, int m);

inline static int chameleon_getblkldd_ccrb(const CHAM_desc_t *A, int m);

/**
 *  Data distributions
 */
int chameleon_getrankof_2d(const CHAM_desc_t *desc, int m, int n);

int chameleon_getrankof_2d_diag(const CHAM_desc_t *desc, int m, int n);

int chameleon_desc_init(CHAM_desc_t *desc, void *mat,
                        cham_flttype_t dtyp, int mb, int nb, int bsiz,
                        int lm, int ln, int i, int j,
                        int m, int n, int p, int q,
                        void *(*get_blkaddr)(const CHAM_desc_t *, int, int),
                        int   (*get_blkldd)(const CHAM_desc_t *, int),
                        int   (*get_rankof)(const CHAM_desc_t *, int, int));

CHAM_desc_t *chameleon_desc_submatrix(CHAM_desc_t *descA, int i, int j, int m, int n);

void chameleon_desc_destroy(CHAM_desc_t *desc);

int chameleon_desc_check(const CHAM_desc_t *desc);

#define BLKLDD(A, k) A->get_blkldd( A, k )

/**
 *  Internal function to return address of block (m,n) with m,n = block indices
 */
inline static void *chameleon_getaddr_ccrb(const CHAM_desc_t *A, int m, int n) {
    size_t mm = m + A->i / A->mb;
    size_t nn = n + A->j / A->nb;
    size_t eltsize = CHAMELEON_Element_Size(A->dtyp);
    size_t offset = 0;

#if defined(CHAMELEON_USE_MPI)
    assert( A->myrank == A->get_rankof( A, mm, nn) );
    mm = mm / A->p;
    nn = nn / A->q;
#endif

    if (mm < (size_t) (A->llm1)) {
        if (nn < (size_t) (A->lln1))
            offset = (size_t) (A->bsiz) * (mm + (size_t) (A->llm1) * nn);
        else
            offset = A->A12 + ((size_t) (A->mb * (A->lln % A->nb)) * mm);
    } else {
        if (nn < (size_t) (A->lln1))
            offset = A->A21 + ((size_t) ((A->llm % A->mb) * A->nb) * nn);
        else
            offset = A->A22;
    }

    return (void *) ((intptr_t) A->mat + (offset * eltsize));
}

/**
 *  Internal function to return address of block (m,n) with m,n = block indices
 */
inline static void *chameleon_getaddr_cm(const CHAM_desc_t *A, int m, int n) {
    size_t mm = m + A->i / A->mb;
    size_t nn = n + A->j / A->nb;
    size_t eltsize = CHAMELEON_Element_Size(A->dtyp);
    size_t offset = 0;

#if defined(CHAMELEON_USE_MPI)
    assert( A->myrank == A->get_rankof( A, mm, nn) );
    mm = mm / A->p;
    nn = nn / A->q;
#endif

    offset = (size_t) (A->llm * A->nb) * nn + (size_t) (A->mb) * mm;
    return (void *) ((intptr_t) A->mat + (offset * eltsize));
}

/**
 *  Internal function to return address of block (m,n) with m,n = block indices
 */
inline static void *chameleon_getaddr_diag(const CHAM_desc_t *A, int m, int n) {
    assert(m == n);
    return chameleon_getaddr_ccrb(A, m, 0);
}

/**
 *  Internal function to return address of block (m,n) with m,n = block indices
 *  This version lets the runtime allocate on-demand.
 */
inline static void *chameleon_getaddr_null(const CHAM_desc_t *A, int m, int n) {
    (void) A;
    (void) m;
    (void) n;
    return NULL;
}

/**
 *  Internal function to return address of block (m,n) with m,n = block indices
 */
inline static CHAM_tile_t *chameleon_desc_gettile(const CHAM_desc_t *A, int m, int n) {
    size_t mm = m + A->i / A->mb;
    size_t nn = n + A->j / A->nb;
    size_t offset = 0;

    assert(A->tiles != NULL);

    offset = A->lmt * nn + mm;
    return A->tiles + offset;
}

/**
 *  Internal function to return address of element A(m,n) with m,n = matrix indices
 */
inline static void *chameleon_geteltaddr(const CHAM_desc_t *A, int m, int n, int eltsize) // Not used anywhere ?!
{
    size_t mm = (m + A->i) / A->mb;
    size_t nn = (n + A->j) / A->nb;
    size_t offset = 0;

#if defined(CHAMELEON_USE_MPI)
    assert( A->myrank == A->get_rankof( A, mm, nn) );
    mm = mm / A->p;
    nn = nn / A->q;
#endif

    if (mm < (size_t) (A->llm1)) {
        if (nn < (size_t) (A->lln1))
            offset = A->bsiz * (mm + A->llm1 * nn) + m % A->mb + A->mb * (n % A->nb);
        else
            offset = A->A12 + (A->mb * (A->lln % A->nb) * mm) + m % A->mb + A->mb * (n % A->nb);
    } else {
        if (nn < (size_t) (A->lln1))
            offset = A->A21 + ((A->llm % A->mb) * A->nb * nn) + m % A->mb + (A->llm % A->mb) * (n % A->nb);
        else
            offset = A->A22 + m % A->mb + (A->llm % A->mb) * (n % A->nb);
    }
    return (void *) ((intptr_t) A->mat + (offset * eltsize));
}

/**
 *  Internal function to return the leading dimension of element A(m,*) with m,n = block indices
 */
inline static int chameleon_getblkldd_ccrb(const CHAM_desc_t *A, int m) {
    int mm = m + A->i / A->mb;
    return (((mm + 1) == A->lmt) && ((A->lm % A->mb) != 0)) ? A->lm % A->mb : A->mb;
}

inline static int chameleon_getblkldd_cm(const CHAM_desc_t *A, int m) {
    (void) m;
    return A->llm;
}

/**
 * Detect if the tile is local or not
 */
inline static int chameleon_desc_islocal(const CHAM_desc_t *A, int m, int n) {
#if defined(CHAMELEON_USE_MPI)
    return (A->myrank == A->get_rankof(A, m, n));
#else
    (void) A;
    (void) m;
    (void) n;
    return 1;
#endif /* defined(CHAMELEON_USE_MPI) */
}

/**
 * Declare data accesses of codelets using these macros, for instance:
 * CHAMELEON_BEGIN_ACCESS_DECLARATION
 * CHAMELEON_ACCESS_R(A, Am, An)
 * CHAMELEON_ACCESS_R(B, Bm, Bn)
 * CHAMELEON_ACCESS_RW(C, Cm, Cn)
 * CHAMELEON_END_ACCESS_DECLARATION
 */
#define CHAMELEON_BEGIN_ACCESS_DECLARATION { \
    unsigned __chameleon_need_submit = 0; \
    RUNTIME_BEGIN_ACCESS_DECLARATION

#define CHAMELEON_ACCESS_R(A, Am, An) do { \
    if (chameleon_desc_islocal(A, Am, An)) __chameleon_need_submit = 1; \
    RUNTIME_ACCESS_R(A, Am, An); \
} while(0)

#define CHAMELEON_ACCESS_W(A, Am, An) do { \
    if (chameleon_desc_islocal(A, Am, An)) __chameleon_need_submit = 1; \
    RUNTIME_ACCESS_W(A, Am, An); \
} while(0)

#define CHAMELEON_ACCESS_RW(A, Am, An) do { \
    if (chameleon_desc_islocal(A, Am, An)) __chameleon_need_submit = 1; \
    RUNTIME_ACCESS_RW(A, Am, An); \
} while(0)

#define CHAMELEON_RANK_CHANGED(rank) do {\
    __chameleon_need_submit = 1; \
    RUNTIME_RANK_CHANGED(rank); \
} while (0)

#define CHAMELEON_END_ACCESS_DECLARATION \
    RUNTIME_END_ACCESS_DECLARATION; \
    if (!__chameleon_need_submit) return; \
}

#ifdef __cplusplus
}
#endif

#endif /* _chameleon_descriptor_h_ */
