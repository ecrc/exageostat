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
 * @file compute_z.h
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
 * @precisions normal z -> c d s
 *
 **/

/***************************************************************************//**
 *  Macro for matrix conversion / Lapack interface
 **/
#define morse_zdesc_alloc_diag(descA, mb, nb, lm, ln, i, j, m, n, p, q) \
    descA = morse_desc_init_diag(                                       \
        MorseComplexDouble, (mb), (nb), ((mb)*(nb)),                    \
        (m), (n), (i), (j), (m), (n), p, q);                            \
    morse_desc_mat_alloc( &(descA) );

#define morse_zdesc_alloc( descA, mb, nb, lm, ln, i, j, m, n, free)     \
    descA = morse_desc_init(                                            \
        MorseComplexDouble, (mb), (nb), ((mb)*(nb)),                    \
        (m), (n), (i), (j), (m), (n), 1, 1);                            \
    if ( morse_desc_mat_alloc( &(descA) ) ) {                           \
        morse_error( __func__, "morse_desc_mat_alloc() failed");        \
        {free;};                                                        \
        return MORSE_ERR_OUT_OF_RESOURCES;                              \
    }

#define morse_zooplap2tile( descA, A, mb, nb, lm, ln, i, j, m, n, seq, req, free) \
    descA = morse_desc_init(                                           \
        MorseComplexDouble, (mb), (nb), ((mb)*(nb)),                   \
        (lm), (ln), (i), (j), (m), (n), 1, 1);                         \
    if ( morse_desc_mat_alloc( &(descA) ) ) {                          \
        morse_error( __func__, "morse_desc_mat_alloc() failed");        \
        {free;};                                                        \
        return MORSE_ERR_OUT_OF_RESOURCES;                             \
    }                                                                   \
    morse_pzlapack_to_tile(A, lm, &descA, seq, req);

#define morse_ziplap2tile( descA, A, mb, nb, lm, ln, i, j, m, n, seq, req) \
    descA = morse_desc_init(                                         \
        MorseComplexDouble, (mb), (nb), ((mb)*(nb)),                 \
        (lm), (ln), (i), (j), (m), (n), 1, 1);                        \
    descA.mat = A;                                                    \
    MORSE_zgecfi_Async((lm), (ln), (A), MorseCM, (mb), (nb),        \
                        MorseCCRB, (mb), (nb), (seq), (req));


#define morse_zooptile2lap( descA, A, mb, nb, lm, ln, seq, req)    \
    morse_pztile_to_lapack(&descA, A, lm, seq, req);

#define morse_ziptile2lap( descA, A, mb, nb, lm, ln, seq, req)         \
    MORSE_zgecfi_Async((lm), (ln), (A), MorseCCRB, (mb), (nb),        \
                        MorseCM, (mb), (nb), (seq), (req));

/***************************************************************************//**
 *  Declarations of internal sequential functions
 **/
int morse_zshift(MORSE_context_t *morse, int m, int n, MORSE_Complex64_t *A,
                  int nprob, int me, int ne, int L,
                  MORSE_sequence_t *sequence, MORSE_request_t *request);

/***************************************************************************//**
 *  Declarations of parallel functions (dynamic scheduling) - alphabetical order
 **/
void morse_pzbarrier_pnl2tl(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzbarrier_row2tl(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzbarrier_tl2pnl(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzbarrier_tl2row(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgebrd_gb2bd(MORSE_enum uplo, MORSE_desc_t *A, double *D, double *E, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgebrd_ge2gb(MORSE_desc_t A, MORSE_desc_t T, MORSE_desc_t D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgelqf(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgelqfrh(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgemm(MORSE_enum transA, MORSE_enum transB, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_Complex64_t beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgeqrf(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgeqrfrh(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgetmi2(MORSE_enum idep, MORSE_enum odep, MORSE_enum storev, int m, int n, int mb, int nb, MORSE_Complex64_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgetrf_incpiv(MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgetrf_nopiv(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgetrf_reclap(MORSE_desc_t *A, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgetrf_rectil(MORSE_desc_t *A, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzhegst(MORSE_enum itype, MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
#ifdef COMPLEX
void morse_pzhemm(MORSE_enum side, MORSE_enum uplo, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_Complex64_t beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzherk(MORSE_enum uplo, MORSE_enum trans, double alpha, MORSE_desc_t *A, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzher2k(MORSE_enum uplo, MORSE_enum trans, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
#endif
void morse_pzhetrd_he2hb(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *E, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzlacpy(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzlag2c(MORSE_desc_t *A, MORSE_desc_t *SB, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzlange(MORSE_enum norm, MORSE_desc_t *A, double *result, MORSE_sequence_t *sequence, MORSE_request_t *request);
#ifdef COMPLEX
void morse_pzlanhe(MORSE_enum norm, MORSE_enum uplo, MORSE_desc_t *A, double *result, MORSE_sequence_t *sequence, MORSE_request_t *request);
#endif
void morse_pzlansy(MORSE_enum norm, MORSE_enum uplo, MORSE_desc_t *A, double *result, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzlantr(MORSE_enum norm, MORSE_enum uplo, MORSE_enum diag, MORSE_desc_t *A, double *result, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzlascal(MORSE_enum uplo, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzlaset( MORSE_enum uplo, MORSE_Complex64_t alpha, MORSE_Complex64_t beta, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzlaset2(MORSE_enum uplo, MORSE_Complex64_t alpha,                          MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzlaswp(MORSE_desc_t *B, int *IPIV, int inc, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzlaswpc(MORSE_desc_t *B, int *IPIV, int inc, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzlauum(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
#ifdef COMPLEX
void morse_pzplghe(double bump, MORSE_enum uplo, MORSE_desc_t *A, unsigned long long int seed, MORSE_sequence_t *sequence, MORSE_request_t *request );
#endif
void morse_pzplgsy(MORSE_Complex64_t bump, MORSE_enum uplo, MORSE_desc_t *A, unsigned long long int seed, MORSE_sequence_t *sequence, MORSE_request_t *request );
void morse_pzplrnt(MORSE_desc_t *A, unsigned long long int seed, MORSE_sequence_t *sequence, MORSE_request_t *request );
void morse_pzpotrf(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzpotrimm(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzshift(int, int, int, MORSE_Complex64_t *, int *, int, int, int, MORSE_sequence_t*, MORSE_request_t*);
void morse_pzsymm(MORSE_enum side, MORSE_enum uplo, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_Complex64_t beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzsyrk(MORSE_enum uplo, MORSE_enum trans, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_Complex64_t beta,  MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzsyr2k(MORSE_enum uplo, MORSE_enum trans, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_Complex64_t beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzsytrf(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pztile2band(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *descAB, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pztpgqrt( int L, MORSE_desc_t *V1, MORSE_desc_t *T1, MORSE_desc_t *V2, MORSE_desc_t *T2, MORSE_desc_t *Q1, MORSE_desc_t *Q2, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request );
void morse_pztpqrt( int L, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request );
void morse_pztradd(MORSE_enum uplo, MORSE_enum trans, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_Complex64_t beta, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pztrmm(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pztrsm(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pztrsmpl(MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *L, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pztrsmrv(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *W, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pztrtri(MORSE_enum uplo, MORSE_enum diag, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzungbr(MORSE_enum side, MORSE_desc_t *A, MORSE_desc_t *O, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzungbrrh(MORSE_enum side, MORSE_desc_t *A, MORSE_desc_t *O, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzungqr(MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzungqrrh(MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_desc_t *D,int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzunglq(MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzunglqrh(MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzungtr(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzunmqr(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzunmqrrh(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzunmlq(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzunmlqrh(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzbuild( MORSE_enum uplo, MORSE_desc_t *A, void *user_data, void* user_build_callback, MORSE_sequence_t *sequence, MORSE_request_t *request );

void morse_pzgelqf_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzgeqrf_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzunmlq_param(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans,
                         MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzunmqr_param(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans,
                         MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzunglq_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *Q,
                         MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pzungqr_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *Q,
                         MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
