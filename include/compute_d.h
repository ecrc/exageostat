/**
 *
 * @file compute_d.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon computational functions header
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2018-11-11
 * @generated d Wed Oct 31 11:23:28 2018
 *
 */
/**
 *  LAPACK/Tile Descriptor accesses
 */
#define MorseDescInput  1
#define MorseDescOutput 2
#define MorseDescInout  (MorseDescInput | MorseDescOutput)

/**
 *  Macro for matrix conversion / Lapack interface
 */
#define morse_ddesc_alloc_diag( descA, mb, nb, lm, ln, i, j, m, n, p, q) \
    descA = morse_desc_init_diag(                                       \
        MorseRealDouble, (mb), (nb), ((mb)*(nb)),                    \
        (m), (n), (i), (j), (m), (n), p, q);                            \
    morse_desc_mat_alloc( &(descA) );                                   \
    RUNTIME_desc_create( &(descA) );

#define morse_ddesc_alloc( descA, mb, nb, lm, ln, i, j, m, n, free)     \
    descA = morse_desc_init(                                            \
        MorseRealDouble, (mb), (nb), ((mb)*(nb)),                    \
        (m), (n), (i), (j), (m), (n), 1, 1);                            \
    if ( morse_desc_mat_alloc( &(descA) ) ) {                           \
        morse_error( __func__, "morse_desc_mat_alloc() failed");        \
        {free;};                                                        \
        return MORSE_ERR_OUT_OF_RESOURCES;                              \
    }                                                                   \
    RUNTIME_desc_create( &(descA) );

/**
 *  Declarations of internal sequential functions
 */
int morse_dshift(MORSE_context_t *morse, int m, int n, double *A,
                  int nprob, int me, int ne, int L,
                  MORSE_sequence_t *sequence, MORSE_request_t *request);

/**
 *  Declarations of parallel functions (dynamic scheduling) - alphabetical order
 */
void morse_pdbarrier_pnl2tl(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdbarrier_row2tl(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdbarrier_tl2pnl(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdbarrier_tl2row(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdgebrd_gb2bd(MORSE_enum uplo, MORSE_desc_t *A, double *D, double *E, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdgebrd_ge2gb(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdgelqf(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdgelqfrh(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdgemm(MORSE_enum transA, MORSE_enum transB, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdgeqrf(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdgeqrfrh(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdgetmi2(MORSE_enum idep, MORSE_enum odep, MORSE_enum storev, int m, int n, int mb, int nb, double *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdgetrf_incpiv(MORSE_desc_t *A, MORSE_desc_t *L, MORSE_desc_t *D, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdgetrf_nopiv(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdgetrf_reclap(MORSE_desc_t *A, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdgetrf_rectil(MORSE_desc_t *A, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdsygst(MORSE_enum itype, MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdsymm(MORSE_enum side, MORSE_enum uplo, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdsyrk(MORSE_enum uplo, MORSE_enum trans, double alpha, MORSE_desc_t *A, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdsyr2k(MORSE_enum uplo, MORSE_enum trans, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdsytrd_sy2sb(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *E, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdlacpy(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdlag2s(MORSE_desc_t *A, MORSE_desc_t *SB, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdlange(MORSE_enum norm, MORSE_desc_t *A, double *result, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdlansy(MORSE_enum norm, MORSE_enum uplo, MORSE_desc_t *A, double *result, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdlansy(MORSE_enum norm, MORSE_enum uplo, MORSE_desc_t *A, double *result, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdlantr(MORSE_enum norm, MORSE_enum uplo, MORSE_enum diag, MORSE_desc_t *A, double *result, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdlascal(MORSE_enum uplo, double alpha, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdlaset( MORSE_enum uplo, double alpha, double beta, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdlaset2(MORSE_enum uplo, double alpha,                          MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdlaswp(MORSE_desc_t *B, int *IPIV, int inc, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdlaswpc(MORSE_desc_t *B, int *IPIV, int inc, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdlauum(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdplgsy(double bump, MORSE_enum uplo, MORSE_desc_t *A, unsigned long long int seed, MORSE_sequence_t *sequence, MORSE_request_t *request );
void morse_pdplgsy(double bump, MORSE_enum uplo, MORSE_desc_t *A, unsigned long long int seed, MORSE_sequence_t *sequence, MORSE_request_t *request );
void morse_pdplrnt(MORSE_desc_t *A, unsigned long long int seed, MORSE_sequence_t *sequence, MORSE_request_t *request );
void morse_pdpotrf(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdpotrimm(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdshift(int, int, int, double *, int *, int, int, int, MORSE_sequence_t*, MORSE_request_t*);
void morse_pdsymm(MORSE_enum side, MORSE_enum uplo, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdsyrk(MORSE_enum uplo, MORSE_enum trans, double alpha, MORSE_desc_t *A, double beta,  MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdsyr2k(MORSE_enum uplo, MORSE_enum trans, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdsytrf(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdtile2band(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *descAB, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdtpgqrt( int L, MORSE_desc_t *V1, MORSE_desc_t *T1, MORSE_desc_t *V2, MORSE_desc_t *T2, MORSE_desc_t *Q1, MORSE_desc_t *Q2, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request );
void morse_pdtpqrt( int L, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request );
void morse_pdtradd(MORSE_enum uplo, MORSE_enum trans, double alpha, MORSE_desc_t *A, double beta, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdtrmm(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdtrsm(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdtrsmpl(MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *L, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdtrsmrv(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, double alpha, MORSE_desc_t *A, MORSE_desc_t *W, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdtrtri(MORSE_enum uplo, MORSE_enum diag, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdorgbr(MORSE_enum side, MORSE_desc_t *A, MORSE_desc_t *O, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdorgbrrh(MORSE_enum side, MORSE_desc_t *A, MORSE_desc_t *O, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdorgqr(MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdorgqrrh(MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_desc_t *D,int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdorglq(MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdorglqrh(MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdorgtr(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdormqr(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdormqrrh(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdormlq(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_desc_t *D, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdormlqrh(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_desc_t *D, int BS, MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdbuild( MORSE_enum uplo, MORSE_desc_t *A, void *user_data, void* user_build_callback, MORSE_sequence_t *sequence, MORSE_request_t *request );

void morse_pdgelqf_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdgeqrf_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdormlq_param(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans,
                         MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdormqr_param(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans,
                         MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdorglq_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *Q,
                         MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);
void morse_pdorgqr_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *Q,
                         MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request);


/**
 * @brief Internal function to convert the lapack format to tile format in
 * LAPACK interface calls
 */
static inline int
morse_dlap2tile( MORSE_context_t *morse,
                 MORSE_desc_t *descAl, MORSE_desc_t *descAt,
                 MORSE_enum mode, MORSE_enum uplo,
                 double *A, int mb, int nb, int lm, int ln, int m, int n,
                 MORSE_sequence_t *seq, MORSE_request_t *req )
{
    /* Initialize the Lapack descriptor */
    *descAl = morse_desc_init_user( MorseRealDouble, mb, nb, (mb)*(nb),
                                    lm, ln, 0, 0, m, n, 1, 1,
                                    morse_getaddr_cm, morse_getblkldd_cm, NULL  );
    descAl->mat = A;
    descAl->styp = MorseCM;

    /* Initialize the tile descriptor */
    *descAt = morse_desc_init( MorseRealDouble, mb, nb, (mb)*(nb),
                               lm, ln, 0, 0, m, n, 1, 1 );

    if ( MORSE_TRANSLATION == MORSE_OUTOFPLACE ) {
        if ( morse_desc_mat_alloc( descAt ) ) {
            morse_error( "morse_dlap2tile", "morse_desc_mat_alloc() failed");
            return MORSE_ERR_OUT_OF_RESOURCES;
        }

        RUNTIME_desc_create( descAl );
        RUNTIME_desc_create( descAt );

        if ( mode & MorseDescInput ) {
            morse_pdlacpy( uplo, descAl, descAt, seq, req );
        }
    }
    else {
        morse_fatal_error( "morse_dlap2tile", "INPLACE translation not supported yet");
        descAt->mat = A;

        RUNTIME_desc_create( descAl );
        RUNTIME_desc_create( descAt );

        if ( mode & MorseDescInput ) {
            /* MORSE_dgecfi_Async( lm, ln, A, MorseCM, mb, nb, */
            /*                     MorseCCRB, mb, nb, seq, req ); */
        }
        return MORSE_ERR_NOT_SUPPORTED;
    }

    return MORSE_SUCCESS;
}

/**
 * @brief Internal function to convert back the tile format to the lapack format
 * in LAPACK interface calls
 */
static inline int
morse_dtile2lap( MORSE_context_t *morse, MORSE_desc_t *descAl, MORSE_desc_t *descAt,
                 MORSE_enum mode, MORSE_enum uplo, MORSE_sequence_t *seq, MORSE_request_t *req )
{
    if ( MORSE_TRANSLATION == MORSE_OUTOFPLACE ) {
        if ( mode & MorseDescOutput ) {
            morse_pdlacpy( uplo, descAt, descAl, seq, req );
        }
    }
    else {
        morse_fatal_error( "morse_dtile2lap", "INPLACE translation not supported yet");
        if ( mode & MorseDescOutput ) {
            /* MORSE_dgecfi_Async( descAl->lm, descAl->ln, descAl->mat, */
            /*                     MorseCCRB, descAl->mb, descAl->nb,   */
            /*                     MorseCM, descAl->mb, descAl->nb, seq, req ); */
        }
        return MORSE_ERR_NOT_SUPPORTED;
    }
    RUNTIME_desc_flush( descAl, seq );
    RUNTIME_desc_flush( descAt, seq );

    return MORSE_SUCCESS;
}

/**
 * @brief Internal function to cleanup the temporary data from the layout
 * conversions in LAPACK interface calls
 */
static inline void
morse_dtile2lap_cleanup( MORSE_context_t *morse, MORSE_desc_t *descAl, MORSE_desc_t *descAt )
{
    if ( MORSE_TRANSLATION == MORSE_OUTOFPLACE ) {
        morse_desc_mat_free( descAt );
    }
    RUNTIME_desc_destroy( descAl );
    RUNTIME_desc_destroy( descAt );
}
