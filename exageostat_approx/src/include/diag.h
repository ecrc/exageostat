/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file diag.h
 *
 * ExaGeoStat approx computation main functions header file.
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#ifndef _DIAG_H_
#define _DIAG_H_

#include "common.h"


void CHAM_pdpotrf_diag(CHAM_enum uplo, CHAM_desc_t *A, int diag_thick,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

int CHAM_dpotrf_diag(CHAM_enum uplo, int N,
                     double* A, int LDA, int diag_thick);

int CHAM_dpotrf_diag_Tile(CHAM_enum uplo, CHAM_desc_t *A, int diag_thick);

int CHAM_dpotrf_diag_Tile_Async(CHAM_enum uplo, CHAM_desc_t *A, int diag_thick,
                                RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

#endif
