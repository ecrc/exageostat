/**
 *
 * Copyright (c) 2017-2023, King Abdullah University of Science and Technology
 * 			    All rights reserved.
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file psdportf.c
 *
 *  CHAMELEON auxiliary routines
 *  CHAMELEON is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 1.2.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 1.0.0
 * @author Sameh Abdulah
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2019-01-09
 * @generated d Fri Dec  1 14:38:19 2017
 *
 **/
#include "../include/mixed_prec.h"

#define A(m, n) A,  m,  n
#define temp(m) temp,  m, 0

/***************************************************************************//**
 *  Parallel tile Cholesky factorization - dynamic scheduling
 **/
void chameleon_psdpotrf(CHAM_enum uplo, CHAM_desc_t *A, CHAM_desc_t *temp, int diag_thick,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request) {

    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int k, m, n;
    int ldak, ldam, ldan, ldtempk;
    int tempkm, tempkn, tempmm, tempnn;
    size_t ws_host = 0;

    double zone = (double) 1.0;
    double mzone = (double) -1.0;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return;
    RUNTIME_options_init(&options, chamctxt, sequence, request);
    RUNTIME_options_ws_alloc(&options, 0, ws_host);

    /*
     *  ChamLower
     */
    if (uplo == ChamLower) {
        //Convert double to single precision for the lower triangular matrix.
        for (n = 0; n < A->nt; n++) {
            tempnn = n == A->nt - 1 ? A->n - n * A->nb : A->nb;
            ldan = BLKLDD(A, n);
            for (m = n + diag_thick; m < A->mt; m++) {

                tempmm = m == A->mt - 1 ? A->m - m * A->mb : A->mb;
                ldam = BLKLDD(A, m);
                EXAGEOSTAT_TASK_dsconv(
                        &options,
                        tempmm, tempnn, A->mb,
                        A(m, n), ldam,
                        A(n, m), ldan);

            }
        }
        //**********************************************
        for (k = 0; k < A->mt; k++) {
            RUNTIME_iteration_push(chamctxt, k);

            tempkm = k == A->mt - 1 ? A->m - k * A->mb : A->mb;
            tempkn = k == A->nt - 1 ? A->n - k * A->nb : A->nb;
            ldak = BLKLDD(A, k);
            ldtempk = BLKLDD(temp, k);

            options.priority = 2 * A->mt - 2 * k;
            INSERT_TASK_dpotrf(
                    &options,
                    ChamLower, tempkm, A->mb,
                    A(k, k), A->nb * k);

            EXAGEOSTAT_TASK_dsconv(
                    &options,
                    tempkm, tempkm, A->mb,
                    A(k, k), ldak,
                    temp(k), ldtempk);

            for (m = k + 1; m < A->mt; m++) {
                tempmm = m == A->mt - 1 ? A->m - m * A->mb : A->mb;
                ldam = BLKLDD(A, m);

                options.priority = 2 * A->mt - 2 * k - m;
                if (abs(m - k) < diag_thick) {
                    INSERT_TASK_dtrsm(
                            &options,
                            ChamRight, ChamLower, ChamTrans, ChamNonUnit,
                            tempmm, A->mb, A->mb,
                            zone, A(k, k),
                            A(m, k));
                } else {
                    EXAGEOSTAT_TASK_sexageostat_trsm(
                            &options,
                            ChamLeft, ChamUpper, ChamTrans, ChamNonUnit,
                            A->mb, tempmm, A->mb,
                            zone, temp(k), ldak,
                            A(k, m), ldak);

                    EXAGEOSTAT_TASK_sdconv(
                            &options,
                            tempkm, tempmm, A->mb,
                            A(k, m), ldak,
                            A(m, k), ldam);

                }
            }

            RUNTIME_data_flush(sequence, A(k, k));
            RUNTIME_data_flush(sequence, temp(k));  //helping in accuracy issue with MPI
            for (n = k + 1; n < A->nt; n++) {
                tempnn = n == A->nt - 1 ? A->n - n * A->nb : A->nb;
                ldan = BLKLDD(A, n);
                ldak = BLKLDD(A, k);
                options.priority = 2 * A->mt - 2 * k - n;

                INSERT_TASK_dsyrk(
                        &options,
                        ChamLower, ChamNoTrans,
                        tempnn, A->nb, A->mb,
                        -1.0, A(n, k),
                        1.0, A(n, n));

                if (abs(n - k) < diag_thick) {

                    EXAGEOSTAT_TASK_dsconv(
                            &options,
                            tempnn, tempkn, A->mb,
                            A(n, k), ldan,
                            A(k, n), ldak);
                }
                for (m = n + 1; m < A->mt; m++) {
                    tempmm = m == A->mt - 1 ? A->m - m * A->mb : A->mb;
                    ldam = BLKLDD(A, m);

                    options.priority = 2 * A->mt - 2 * k - n - m;

                    if (abs(m - n) < diag_thick) {
                        INSERT_TASK_dgemm(
                                &options,
                                ChamNoTrans, ChamTrans,
                                tempmm, tempnn, A->mb, A->mb,
                                mzone, A(m, k),
                                A(n, k),
                                zone, A(m, n));
                    } else {
                        EXAGEOSTAT_TASK_sexageostat_gemm(
                                &options,
                                ChamTrans, ChamNoTrans,
                                tempnn, tempmm, A->mb, A->mb,
                                mzone, A(k, n), ldak,
                                A(k, m), ldak,
                                zone, A(n, m), ldan);

                    }
                }
                RUNTIME_data_flush(sequence, A(n, k));
                RUNTIME_data_flush(sequence, A(k, n)); //fix pending problem with MPI
            }
            RUNTIME_iteration_pop(chamctxt);
        }
    }
        /*
         *  ChamUpper
         */
    else {
        for (k = 0; k < A->nt; k++) {
            RUNTIME_iteration_push(chamctxt, k);

            tempkm = k == A->nt - 1 ? A->n - k * A->nb : A->nb;
            ldak = BLKLDD(A, k);

            options.priority = 2 * A->nt - 2 * k;
            INSERT_TASK_dpotrf(
                    &options,
                    ChamUpper,
                    tempkm, A->mb,
                    A(k, k), A->nb * k);

            for (n = k + 1; n < A->nt; n++) {
                tempnn = n == A->nt - 1 ? A->n - n * A->nb : A->nb;

                options.priority = 2 * A->nt - 2 * k - n;
                INSERT_TASK_dtrsm(
                        &options,
                        ChamLeft, ChamUpper, ChamTrans, ChamNonUnit,
                        A->mb, tempnn, A->mb,
                        zone, A(k, k),
                        A(k, n));
            }
            RUNTIME_data_flush(sequence, A(k, k));

            for (m = k + 1; m < A->mt; m++) {
                tempmm = m == A->mt - 1 ? A->m - m * A->mb : A->mb;
                ldam = BLKLDD(A, m);

                options.priority = 2 * A->nt - 2 * k - m;
                INSERT_TASK_dsyrk(
                        &options,
                        ChamUpper, ChamTrans,
                        tempmm, A->mb, A->mb,
                        -1.0, A(k, m),
                        1.0, A(m, m));

                for (n = m + 1; n < A->nt; n++) {
                    tempnn = n == A->nt - 1 ? A->n - n * A->nb : A->nb;

                    options.priority = 2 * A->nt - 2 * k - n - m;
                    INSERT_TASK_dgemm(
                            &options,
                            ChamTrans, ChamNoTrans,
                            tempmm, tempnn, A->mb, A->mb,
                            mzone, A(k, m),
                            A(k, n),
                            zone, A(m, n));
                }
                RUNTIME_data_flush(sequence, A(k, m));
            }

            RUNTIME_iteration_pop(chamctxt);
        }
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}