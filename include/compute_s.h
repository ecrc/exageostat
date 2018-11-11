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
 * @file compute_s.h
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2018-11-11
 * @generated s Fri Dec  1 14:38:09 2017
 *
 **/

/***************************************************************************//**
 *  Macro for matrix conversion / Lapack interface
 **/
#define morse_sdesc_alloc_diag(descA, mb, nb, lm, ln, i, j, m, n, p, q) \
    descA = morse_desc_init_diag(                                       \
        MorseRealFloat, (mb), (nb), ((mb)*(nb)),                    \
        (m), (n), (i), (j), (m), (n), p, q);                            \
    morse_desc_mat_alloc( &(descA) );

#define morse_sdesc_alloc( descA, mb, nb, lm, ln, i, j, m, n, free)     \
    descA = morse_desc_init(                                            \
        MorseRealFloat, (mb), (nb), ((mb)*(nb)),                    \
        (m), (n), (i), (j), (m), (n), 1, 1);                            \
    if ( morse_desc_mat_alloc( &(descA) ) ) {                           \
        morse_error( __func__, "morse_desc_mat_alloc() failed");        \
        {free;};                                                        \
        return MORSE_ERR_OUT_OF_RESOURCES;                              \
    }

#define morse_sooplap2tile( descA, A, mb, nb, lm, ln, i, j, m, n, seq, req, free) \
    descA = morse_desc_init(                                           \
        MorseRealFloat, (mb), (nb), ((mb)*(nb)),                   \
        (lm), (ln), (i), (j), (m), (n), 1, 1);                         \
    if ( morse_desc_mat_alloc( &(descA) ) ) {                          \
        morse_error( __func__, "morse_desc_mat_alloc() failed");        \
        {free;};                                                        \
        return MORSE_ERR_OUT_OF_RESOURCES;                             \
    }                                                                   \
    morse_pslapack_to_tile(A, lm, &descA, seq, req);

#define morse_siplap2tile( descA, A, mb, nb, lm, ln, i, j, m, n, seq, req) \
    descA = morse_desc_init(                                         \
        MorseRealFloat, (mb), (nb), ((mb)*(nb)),                 \
        (lm), (ln), (i), (j), (m), (n), 1, 1);                        \
    descA.mat = A;                                                    \
    MORSE_sgecfi_Async((lm), (ln), (A), MorseCM, (mb), (nb),        \
                        MorseCCRB, (mb), (nb), (seq), (req));


#define morse_sooptile2lap( descA, A, mb, nb, lm, ln, seq, req)    \
    morse_pstile_to_lapack(&descA, A, lm, seq, req);

#define morse_siptile2lap( descA, A, mb, nb, lm, ln, seq, req)         \
    MORSE_sgecfi_Async((lm), (ln), (A), MorseCCRB, (mb), (nb),        \
                        MorseCM, (mb), (nb), (seq), (req));

/***************************************************************************//**
 *  Declarations of internal sequential functions
 **/
int morse_sshift(MORSE_context_t *morse, int m, int n, float *A,
                  int nprob, int me, int ne, int L,
                  MORSE_sequence_t *sequence, MORSE_request_t *request);

/***************************************************************************//**
 *  Declarations of parallel functions (dynamic scheduling) - alphabetical order
 **/
void morse_psbarrier_pnl2tl(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psbarrier_row2tl(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psbarrier_tl2pnl(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psbarrier_tl2row(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psgebrd_gb2bd(MORSE_enum uplo, MORSE_desc_t *A, float *D, float *E, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psgebrd_ge2gb(MORSE_desc_t A, MORSE_desc_t T, MORSE_desc_t D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psgelqf(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psgelqfrh(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psgemm(MORSE_enum transA, MORSE_enum transB, float alpha, MORSE_desc_t *A, MORSE_desc_t *B, float beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psgeqrf(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psgeqrfrh(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psgetmi2(MORSE_enum idep, MORSE_enum odep, MORSE_enum storev, int m, int n, int mb, int nb, float *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psgetrf_incpiv(MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psgetrf_nopiv(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psgetrf_reclap(MORSE_desc_t *A, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psgetrf_rectil(MORSE_desc_t *A, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pssygst(MORSE_enum itype, MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
#ifdef COMPLEX
void morse_pssymm(MORSE_enum side, MORSE_enum uplo, float alpha, MORSE_desc_t *A, MORSE_desc_t *B, float beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pssyrk(MORSE_enum uplo, MORSE_enum trans, float alpha, MORSE_desc_t *A, float beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pssyr2k(MORSE_enum uplo, MORSE_enum trans, float alpha, MORSE_desc_t *A, MORSE_desc_t *B, float beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
#endif
void morse_pssytrd_sy2sb(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *E, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pslacpy(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pslag2d(MORSE_desc_t *A, MORSE_desc_t *SB, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pslange(MORSE_enum norm, MORSE_desc_t *A, float *result, MORSE_sequence_t *sequence, MORSE_request_t *request);
#ifdef COMPLEX
void morse_pslansy(MORSE_enum norm, MORSE_enum uplo, MORSE_desc_t *A, float *result, MORSE_sequence_t *sequence, MORSE_request_t *request);
#endif
void morse_pslansy(MORSE_enum norm, MORSE_enum uplo, MORSE_desc_t *A, float *result, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pslantr(MORSE_enum norm, MORSE_enum uplo, MORSE_enum diag, MORSE_desc_t *A, float *result, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pslascal(MORSE_enum uplo, float alpha, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pslaset( MORSE_enum uplo, float alpha, float beta, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pslaset2(MORSE_enum uplo, float alpha,                          MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pslaswp(MORSE_desc_t *B, int *IPIV, int inc, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pslaswpc(MORSE_desc_t *B, int *IPIV, int inc, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pslauum(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
#ifdef COMPLEX
void morse_psplgsy(float bump, MORSE_enum uplo, MORSE_desc_t *A, unsigned long long int seed, MORSE_sequence_t *sequence, MORSE_request_t *request );
#endif
void morse_psplgsy(float bump, MORSE_enum uplo, MORSE_desc_t *A, unsigned long long int seed, MORSE_sequence_t *sequence, MORSE_request_t *request );
void morse_psplrnt(MORSE_desc_t *A, unsigned long long int seed, MORSE_sequence_t *sequence, MORSE_request_t *request );
void morse_pspotrf(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pspotrimm(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psshift(int, int, int, float *, int *, int, int, int, MORSE_sequence_t*, MORSE_request_t*);
void morse_pssymm(MORSE_enum side, MORSE_enum uplo, float alpha, MORSE_desc_t *A, MORSE_desc_t *B, float beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pssyrk(MORSE_enum uplo, MORSE_enum trans, float alpha, MORSE_desc_t *A, float beta,  MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pssyr2k(MORSE_enum uplo, MORSE_enum trans, float alpha, MORSE_desc_t *A, MORSE_desc_t *B, float beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pssytrf(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pstile2band(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *descAB, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pstpgqrt( int L, MORSE_desc_t *V1, MORSE_desc_t *T1, MORSE_desc_t *V2, MORSE_desc_t *T2, MORSE_desc_t *Q1, MORSE_desc_t *Q2, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request );
void morse_pstpqrt( int L, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request );
void morse_pstradd(MORSE_enum uplo, MORSE_enum trans, float alpha, MORSE_desc_t *A, float beta, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pstrmm(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, float alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pstrsm(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, float alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pstrsmpl(MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *L, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pstrsmrv(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, float alpha, MORSE_desc_t *A, MORSE_desc_t *W, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pstrtri(MORSE_enum uplo, MORSE_enum diag, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psorgbr(MORSE_enum side, MORSE_desc_t *A, MORSE_desc_t *O, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psorgbrrh(MORSE_enum side, MORSE_desc_t *A, MORSE_desc_t *O, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psorgqr(MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psorgqrrh(MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_desc_t *D,int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psorglq(MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psorglqrh(MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psorgtr(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psormqr(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psormqrrh(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psormlq(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psormlqrh(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psbuild( MORSE_enum uplo, MORSE_desc_t *A, void *user_data, void* user_build_callback, MORSE_sequence_t *sequence, MORSE_request_t *request );

void morse_psgelqf_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psgeqrf_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psormlq_param(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans,
                         MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psormqr_param(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans,
                         MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psorglq_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *Q,
                         MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_psorgqr_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *Q,
                         MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
