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
 * @file compute_c.h
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
 * @generated c Fri Dec  1 14:38:09 2017
 *
 **/

/***************************************************************************//**
 *  Macro for matrix conversion / Lapack interface
 **/
#define morse_cdesc_alloc_diag(descA, mb, nb, lm, ln, i, j, m, n, p, q) \
    descA = morse_desc_init_diag(                                       \
        MorseComplexFloat, (mb), (nb), ((mb)*(nb)),                    \
        (m), (n), (i), (j), (m), (n), p, q);                            \
    morse_desc_mat_alloc( &(descA) );

#define morse_cdesc_alloc( descA, mb, nb, lm, ln, i, j, m, n, free)     \
    descA = morse_desc_init(                                            \
        MorseComplexFloat, (mb), (nb), ((mb)*(nb)),                    \
        (m), (n), (i), (j), (m), (n), 1, 1);                            \
    if ( morse_desc_mat_alloc( &(descA) ) ) {                           \
        morse_error( __func__, "morse_desc_mat_alloc() failed");        \
        {free;};                                                        \
        return MORSE_ERR_OUT_OF_RESOURCES;                              \
    }

#define morse_cooplap2tile( descA, A, mb, nb, lm, ln, i, j, m, n, seq, req, free) \
    descA = morse_desc_init(                                           \
        MorseComplexFloat, (mb), (nb), ((mb)*(nb)),                   \
        (lm), (ln), (i), (j), (m), (n), 1, 1);                         \
    if ( morse_desc_mat_alloc( &(descA) ) ) {                          \
        morse_error( __func__, "morse_desc_mat_alloc() failed");        \
        {free;};                                                        \
        return MORSE_ERR_OUT_OF_RESOURCES;                             \
    }                                                                   \
    morse_pclapack_to_tile(A, lm, &descA, seq, req);

#define morse_ciplap2tile( descA, A, mb, nb, lm, ln, i, j, m, n, seq, req) \
    descA = morse_desc_init(                                         \
        MorseComplexFloat, (mb), (nb), ((mb)*(nb)),                 \
        (lm), (ln), (i), (j), (m), (n), 1, 1);                        \
    descA.mat = A;                                                    \
    MORSE_cgecfi_Async((lm), (ln), (A), MorseCM, (mb), (nb),        \
                        MorseCCRB, (mb), (nb), (seq), (req));


#define morse_cooptile2lap( descA, A, mb, nb, lm, ln, seq, req)    \
    morse_pctile_to_lapack(&descA, A, lm, seq, req);

#define morse_ciptile2lap( descA, A, mb, nb, lm, ln, seq, req)         \
    MORSE_cgecfi_Async((lm), (ln), (A), MorseCCRB, (mb), (nb),        \
                        MorseCM, (mb), (nb), (seq), (req));

/***************************************************************************//**
 *  Declarations of internal sequential functions
 **/
int morse_cshift(MORSE_context_t *morse, int m, int n, MORSE_Complex32_t *A,
                  int nprob, int me, int ne, int L,
                  MORSE_sequence_t *sequence, MORSE_request_t *request);

/***************************************************************************//**
 *  Declarations of parallel functions (dynamic scheduling) - alphabetical order
 **/
void morse_pcbarrier_pnl2tl(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcbarrier_row2tl(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcbarrier_tl2pnl(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcbarrier_tl2row(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcgebrd_gb2bd(MORSE_enum uplo, MORSE_desc_t *A, float *D, float *E, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcgebrd_ge2gb(MORSE_desc_t A, MORSE_desc_t T, MORSE_desc_t D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcgelqf(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcgelqfrh(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcgemm(MORSE_enum transA, MORSE_enum transB, MORSE_Complex32_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_Complex32_t beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcgeqrf(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcgeqrfrh(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcgetmi2(MORSE_enum idep, MORSE_enum odep, MORSE_enum storev, int m, int n, int mb, int nb, MORSE_Complex32_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcgetrf_incpiv(MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcgetrf_nopiv(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcgetrf_reclap(MORSE_desc_t *A, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcgetrf_rectil(MORSE_desc_t *A, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pchegst(MORSE_enum itype, MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
#ifdef COMPLEX
void morse_pchemm(MORSE_enum side, MORSE_enum uplo, MORSE_Complex32_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_Complex32_t beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcherk(MORSE_enum uplo, MORSE_enum trans, float alpha, MORSE_desc_t *A, float beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcher2k(MORSE_enum uplo, MORSE_enum trans, MORSE_Complex32_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, float beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
#endif
void morse_pchetrd_he2hb(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *E, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pclacpy(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pclag2z(MORSE_desc_t *A, MORSE_desc_t *SB, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pclange(MORSE_enum norm, MORSE_desc_t *A, float *result, MORSE_sequence_t *sequence, MORSE_request_t *request);
#ifdef COMPLEX
void morse_pclanhe(MORSE_enum norm, MORSE_enum uplo, MORSE_desc_t *A, float *result, MORSE_sequence_t *sequence, MORSE_request_t *request);
#endif
void morse_pclansy(MORSE_enum norm, MORSE_enum uplo, MORSE_desc_t *A, float *result, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pclantr(MORSE_enum norm, MORSE_enum uplo, MORSE_enum diag, MORSE_desc_t *A, float *result, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pclascal(MORSE_enum uplo, MORSE_Complex32_t alpha, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pclaset( MORSE_enum uplo, MORSE_Complex32_t alpha, MORSE_Complex32_t beta, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pclaset2(MORSE_enum uplo, MORSE_Complex32_t alpha,                          MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pclaswp(MORSE_desc_t *B, int *IPIV, int inc, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pclaswpc(MORSE_desc_t *B, int *IPIV, int inc, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pclauum(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
#ifdef COMPLEX
void morse_pcplghe(float bump, MORSE_enum uplo, MORSE_desc_t *A, unsigned long long int seed, MORSE_sequence_t *sequence, MORSE_request_t *request );
#endif
void morse_pcplgsy(MORSE_Complex32_t bump, MORSE_enum uplo, MORSE_desc_t *A, unsigned long long int seed, MORSE_sequence_t *sequence, MORSE_request_t *request );
void morse_pcplrnt(MORSE_desc_t *A, unsigned long long int seed, MORSE_sequence_t *sequence, MORSE_request_t *request );
void morse_pcpotrf(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcpotrimm(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcshift(int, int, int, MORSE_Complex32_t *, int *, int, int, int, MORSE_sequence_t*, MORSE_request_t*);
void morse_pcsymm(MORSE_enum side, MORSE_enum uplo, MORSE_Complex32_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_Complex32_t beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcsyrk(MORSE_enum uplo, MORSE_enum trans, MORSE_Complex32_t alpha, MORSE_desc_t *A, MORSE_Complex32_t beta,  MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcsyr2k(MORSE_enum uplo, MORSE_enum trans, MORSE_Complex32_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_Complex32_t beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcsytrf(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pctile2band(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *descAB, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pctpgqrt( int L, MORSE_desc_t *V1, MORSE_desc_t *T1, MORSE_desc_t *V2, MORSE_desc_t *T2, MORSE_desc_t *Q1, MORSE_desc_t *Q2, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request );
void morse_pctpqrt( int L, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request );
void morse_pctradd(MORSE_enum uplo, MORSE_enum trans, MORSE_Complex32_t alpha, MORSE_desc_t *A, MORSE_Complex32_t beta, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pctrmm(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, MORSE_Complex32_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pctrsm(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, MORSE_Complex32_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pctrsmpl(MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *L, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pctrsmrv(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, MORSE_Complex32_t alpha, MORSE_desc_t *A, MORSE_desc_t *W, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pctrtri(MORSE_enum uplo, MORSE_enum diag, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcungbr(MORSE_enum side, MORSE_desc_t *A, MORSE_desc_t *O, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcungbrrh(MORSE_enum side, MORSE_desc_t *A, MORSE_desc_t *O, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcungqr(MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcungqrrh(MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_desc_t *D,int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcunglq(MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcunglqrh(MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcungtr(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcunmqr(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcunmqrrh(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcunmlq(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcunmlqrh(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcbuild( MORSE_enum uplo, MORSE_desc_t *A, void *user_data, void* user_build_callback, MORSE_sequence_t *sequence, MORSE_request_t *request );

void morse_pcgelqf_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcgeqrf_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcunmlq_param(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans,
                         MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcunmqr_param(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans,
                         MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcunglq_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *Q,
                         MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pcungqr_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *Q,
                         MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
