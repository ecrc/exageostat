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
 * @version 1.2.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2022-11-09
 * @generated d Wed Oct 31 11:23:28 2018
 *
 */
/**
 *  LAPACK/Tile Descriptor accesses
 */
#include <chameleon/constants.h>
#include <global.h>
#include <common.h>

#define ChamDescInput  1
#define ChamDescOutput 2
#define ChamDescInout  (ChamDescInput | ChamDescOutput)

/**
 *  Macro for matrix conversion / Lapack interface
 */
#define chameleon_ddesc_alloc_diag(descA, mb, nb, lm, ln, i, j, m, n, p, q) \
    descA = chameleon_desc_init_diag(                                       \
        ChamRealDouble, (mb), (nb), ((mb)*(nb)),                    \
        (m), (n), (i), (j), (m), (n), p, q);                            \
    chameleon_desc_mat_alloc( &(descA) );                                   \
    RUNTIME_desc_create( &(descA) );

#define chameleon_ddesc_alloc(descA, mb, nb, lm, ln, i, j, m, n, free)     \
    descA = chameleon_desc_init(                                            \
        ChamRealDouble, (mb), (nb), ((mb)*(nb)),                    \
        (m), (n), (i), (j), (m), (n), 1, 1);                            \
    if ( chameleon_desc_mat_alloc( &(descA) ) ) {                           \
        chameleon_error( __func__, "chameleon_desc_mat_alloc() failed");        \
        {free;};                                                        \
        return CHAMELEON_ERR_OUT_OF_RESOURCES;                              \
    }                                                                   \
    RUNTIME_desc_create( &(descA) );

/**
 *  Declarations of internal sequential functions
 */
int chameleon_dshift(CHAM_context_t *CHAM, int m, int n, double* A,
                     int nprob, int me, int ne, int L,
                     RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

/**
 *  Declarations of parallel functions (dynamic scheduling) - alphabetical order
 */
void chameleon_pdbarrier_pnl2tl(CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdbarrier_row2tl(CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdbarrier_tl2pnl(CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdbarrier_tl2row(CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdgebrd_gb2bd(CHAM_enum uplo, CHAM_desc_t *A, double* D, double* E, CHAM_desc_t *T,
                             RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdgebrd_ge2gb(CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence,
                             RUNTIME_request_t *request);

void chameleon_pdgelqf(CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence,
                       RUNTIME_request_t *request);

void chameleon_pdgelqfrh(CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D, int BS, RUNTIME_sequence_t *sequence,
                         RUNTIME_request_t *request);

void chameleon_pdgemm(CHAM_enum transA, CHAM_enum transB, double alpha, CHAM_desc_t *A, CHAM_desc_t *B, double beta,
                      CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdgeqrf(CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence,
                       RUNTIME_request_t *request);

void chameleon_pdgeqrfrh(CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D, int BS, RUNTIME_sequence_t *sequence,
                         RUNTIME_request_t *request);

void chameleon_pdgetmi2(CHAM_enum idep, CHAM_enum odep, CHAM_enum storev, int m, int n, int mb, int nb, double* A,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdgetrf_incpiv(CHAM_desc_t *A, CHAM_desc_t *L, CHAM_desc_t *D, int *IPIV, RUNTIME_sequence_t *sequence,
                              RUNTIME_request_t *request);

void chameleon_pdgetrf_nopiv(CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdgetrf_reclap(CHAM_desc_t *A, int *IPIV, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdgetrf_rectil(CHAM_desc_t *A, int *IPIV, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdsygst(CHAM_enum itype, CHAM_enum uplo, CHAM_desc_t *A, CHAM_desc_t *B, RUNTIME_sequence_t *sequence,
                       RUNTIME_request_t *request);

void chameleon_pdsymm(CHAM_enum side, CHAM_enum uplo, double alpha, CHAM_desc_t *A, CHAM_desc_t *B, double beta,
                      CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdsyrk(CHAM_enum uplo, CHAM_enum trans, double alpha, CHAM_desc_t *A, double beta, CHAM_desc_t *C,
                      RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdsyr2k(CHAM_enum uplo, CHAM_enum trans, double alpha, CHAM_desc_t *A, CHAM_desc_t *B, double beta,
                       CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void
chameleon_pdsytrd_sy2sb(CHAM_enum uplo, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *E, RUNTIME_sequence_t *sequence,
                        RUNTIME_request_t *request);

void chameleon_pdlacpy(CHAM_enum uplo, CHAM_desc_t *A, CHAM_desc_t *B, RUNTIME_sequence_t *sequence,
                       RUNTIME_request_t *request);

void chameleon_pdlag2s(CHAM_desc_t *A, CHAM_desc_t *SB, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdlange(CHAM_enum norm, CHAM_desc_t *A, double* result, RUNTIME_sequence_t *sequence,
                       RUNTIME_request_t *request);

void chameleon_pdlansy(CHAM_enum norm, CHAM_enum uplo, CHAM_desc_t *A, double* result, RUNTIME_sequence_t *sequence,
                       RUNTIME_request_t *request);

void chameleon_pdlansy(CHAM_enum norm, CHAM_enum uplo, CHAM_desc_t *A, double* result, RUNTIME_sequence_t *sequence,
                       RUNTIME_request_t *request);

void chameleon_pdlantr(CHAM_enum norm, CHAM_enum uplo, CHAM_enum diag, CHAM_desc_t *A, double* result,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdlascal(CHAM_enum uplo, double alpha, CHAM_desc_t *A, RUNTIME_sequence_t *sequence,
                        RUNTIME_request_t *request);

void chameleon_pdlaset(CHAM_enum uplo, double alpha, double beta, CHAM_desc_t *A, RUNTIME_sequence_t *sequence,
                       RUNTIME_request_t *request);

void chameleon_pdlaset2(CHAM_enum uplo, double alpha, CHAM_desc_t *A, RUNTIME_sequence_t *sequence,
                        RUNTIME_request_t *request);

void chameleon_pdlaswp(CHAM_desc_t *B, int *IPIV, int inc, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdlaswpc(CHAM_desc_t *B, int *IPIV, int inc, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdlauum(CHAM_enum uplo, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdplgsy(double bump, CHAM_enum uplo, CHAM_desc_t *A, unsigned long long int seed,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdplgsy(double bump, CHAM_enum uplo, CHAM_desc_t *A, unsigned long long int seed,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdplrnt(CHAM_desc_t *A, unsigned long long int seed, RUNTIME_sequence_t *sequence,
                       RUNTIME_request_t *request);

void chameleon_pdpotrf(CHAM_enum uplo, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdpotrimm(CHAM_enum uplo, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *C, RUNTIME_sequence_t *sequence,
                         RUNTIME_request_t *request);

void chameleon_pdshift(int, int, int, double* , int *, int, int, int, RUNTIME_sequence_t *, RUNTIME_request_t *);

void chameleon_pdsymm(CHAM_enum side, CHAM_enum uplo, double alpha, CHAM_desc_t *A, CHAM_desc_t *B, double beta,
                      CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdsyrk(CHAM_enum uplo, CHAM_enum trans, double alpha, CHAM_desc_t *A, double beta, CHAM_desc_t *C,
                      RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdsyr2k(CHAM_enum uplo, CHAM_enum trans, double alpha, CHAM_desc_t *A, CHAM_desc_t *B, double beta,
                       CHAM_desc_t *C, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdsytrf(CHAM_enum uplo, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdtile2band(CHAM_enum uplo, CHAM_desc_t *A, CHAM_desc_t *descAB, RUNTIME_sequence_t *sequence,
                           RUNTIME_request_t *request);

void chameleon_pdtpgqrt(int L, CHAM_desc_t *V1, CHAM_desc_t *T1, CHAM_desc_t *V2, CHAM_desc_t *T2, CHAM_desc_t *Q1,
                        CHAM_desc_t *Q2, CHAM_desc_t *D, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdtpqrt(int L, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T, RUNTIME_sequence_t *sequence,
                       RUNTIME_request_t *request);

void chameleon_pdtradd(CHAM_enum uplo, CHAM_enum trans, double alpha, CHAM_desc_t *A, double beta, CHAM_desc_t *B,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdtrmm(CHAM_enum side, CHAM_enum uplo, CHAM_enum transA, CHAM_enum diag, double alpha, CHAM_desc_t *A,
                      CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdtrsm(CHAM_enum side, CHAM_enum uplo, CHAM_enum transA, CHAM_enum diag, double alpha, CHAM_desc_t *A,
                      CHAM_desc_t *B, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdtrsmpl(CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *L, int *IPIV, RUNTIME_sequence_t *sequence,
                        RUNTIME_request_t *request);

void chameleon_pdtrsmrv(CHAM_enum side, CHAM_enum uplo, CHAM_enum transA, CHAM_enum diag, double alpha, CHAM_desc_t *A,
                        CHAM_desc_t *W, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdtrtri(CHAM_enum uplo, CHAM_enum diag, CHAM_desc_t *A, RUNTIME_sequence_t *sequence,
                       RUNTIME_request_t *request);

void chameleon_pdorgbr(CHAM_enum side, CHAM_desc_t *A, CHAM_desc_t *O, CHAM_desc_t *T, RUNTIME_sequence_t *sequence,
                       RUNTIME_request_t *request);

void chameleon_pdorgbrrh(CHAM_enum side, CHAM_desc_t *A, CHAM_desc_t *O, CHAM_desc_t *T, RUNTIME_sequence_t *sequence,
                         RUNTIME_request_t *request);

void chameleon_pdorgqr(CHAM_desc_t *A, CHAM_desc_t *Q, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence,
                       RUNTIME_request_t *request);

void chameleon_pdorgqrrh(CHAM_desc_t *A, CHAM_desc_t *Q, CHAM_desc_t *T, CHAM_desc_t *D, int BS,
                         RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdorglq(CHAM_desc_t *A, CHAM_desc_t *Q, CHAM_desc_t *T, CHAM_desc_t *D, RUNTIME_sequence_t *sequence,
                       RUNTIME_request_t *request);

void chameleon_pdorglqrh(CHAM_desc_t *A, CHAM_desc_t *Q, CHAM_desc_t *T, CHAM_desc_t *D, int BS,
                         RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdorgtr(CHAM_enum uplo, CHAM_desc_t *A, CHAM_desc_t *Q, CHAM_desc_t *T, RUNTIME_sequence_t *sequence,
                       RUNTIME_request_t *request);

void chameleon_pdormqr(CHAM_enum side, CHAM_enum trans, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T, CHAM_desc_t *D,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void
chameleon_pdormqrrh(CHAM_enum side, CHAM_enum trans, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T, CHAM_desc_t *D,
                    int BS, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdormlq(CHAM_enum side, CHAM_enum trans, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T, CHAM_desc_t *D,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void
chameleon_pdormlqrh(CHAM_enum side, CHAM_enum trans, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T, CHAM_desc_t *D,
                    int BS, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdbuild(CHAM_enum uplo, CHAM_desc_t *A, void *user_data, void *user_build_callback,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void
chameleon_pdgelqf_param(const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void
chameleon_pdgeqrf_param(const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdormlq_param(const libhqr_tree_t *qrtree, CHAM_enum side, CHAM_enum trans,
                             CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                             RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdormqr_param(const libhqr_tree_t *qrtree, CHAM_enum side, CHAM_enum trans,
                             CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                             RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdorglq_param(const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *Q,
                             CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                             RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

void chameleon_pdorgqr_param(const libhqr_tree_t *qrtree, CHAM_desc_t *A, CHAM_desc_t *Q,
                             CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                             RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);


/**
 * @brief Internal function to convert the lapack format to tile format in
 * LAPACK interface calls
 */
static inline int
chameleon_dlap2tile(CHAM_context_t *CHAM,
                    CHAM_desc_t *descAl, CHAM_desc_t *descAt,
                    CHAM_enum mode, CHAM_enum uplo,
                    double* A, int mb, int nb, int lm, int ln, int m, int n,
                    RUNTIME_sequence_t *seq, RUNTIME_request_t *req) {
    /* Initialize the Lapack descriptor */
    void *mat;
    int result = chameleon_desc_init(descAt, mat,
                                     ChamRealDouble, mb, nb, (mb) * (nb),
                                     lm, ln, 0, 0, m, n, 1, 1,
                                     chameleon_getaddr_cm, chameleon_getblkldd_cm, NULL);
    descAl->mat = A;
    descAl->styp = ChamCM;

    /* Initialize the tile descriptor */
    result = chameleon_desc_init(descAt, mat, ChamRealDouble, mb, nb, (mb) * (nb),
                                 lm, ln, 0, 0, m, n, 1, 1, chameleon_getaddr_cm, chameleon_getblkldd_cm, NULL);

    CHAM_context_t *chamctxt;
    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAM_diag_dpotrf", "CHAM not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    if (CHAMELEON_TRANSLATION == ChamOutOfPlace) {


        RUNTIME_desc_create(descAl);
        RUNTIME_desc_create(descAt);

        if (mode & ChamDescInput) {
            chameleon_pdlacpy(uplo, descAl, descAt, seq, req);
        }
    } else {
        chameleon_fatal_error("chameleon_dlap2tile", "INPLACE translation not supported yet");
        descAt->mat = A;

        RUNTIME_desc_create(descAl);
        RUNTIME_desc_create(descAt);

        return CHAMELEON_ERR_NOT_SUPPORTED;
    }

    return CHAMELEON_SUCCESS;
}

/**
 * @brief Internal function to convert back the tile format to the lapack format
 * in LAPACK interface calls
 */
static inline int
chameleon_dtile2lap(CHAM_context_t *CHAM, CHAM_desc_t *descAl, CHAM_desc_t *descAt,
                    CHAM_enum mode, CHAM_enum uplo, RUNTIME_sequence_t *seq, RUNTIME_request_t *req) {
    CHAM_context_t *chamctxt;
    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAM_diag_dpotrf", "CHAM not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    if (CHAMELEON_TRANSLATION == ChamOutOfPlace) {
        if (mode & ChamDescOutput) {
            chameleon_pdlacpy(uplo, descAt, descAl, seq, req);
        }
    } else {
        chameleon_fatal_error("chameleon_dtile2lap", "INPLACE translation not supported yet");
        return CHAMELEON_ERR_NOT_SUPPORTED;
    }
    RUNTIME_desc_flush(descAl, seq);
    RUNTIME_desc_flush(descAt, seq);

    return CHAMELEON_SUCCESS;
}

/**
 * @brief Internal function to cleanup the temporary data from the layout
 * conversions in LAPACK interface calls
 */
static inline void
chameleon_dtile2lap_cleanup(CHAM_context_t *CHAM, CHAM_desc_t *descAl, CHAM_desc_t *descAt) {
    CHAM_context_t *chamctxt;
    chamctxt = chameleon_context_self();


    if (CHAMELEON_TRANSLATION == ChamOutOfPlace) {
        chameleon_desc_destroy(descAt);
    }
    RUNTIME_desc_destroy(descAl);
    RUNTIME_desc_destroy(descAt);
}
