/**
 *
 * Copyright (c) 2017-2018  King Abdullah University of Science and Technology
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
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2018-11-11
 *
 **/
#ifndef _DIAG_H_
#define _DIAG_H_
#include "common.h"


void morse_pdpotrf_diag(MORSE_enum uplo, MORSE_desc_t *A, int diag_thick,
                   MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dpotrf_diag(MORSE_enum uplo, int N,
                  double *A, int LDA, int diag_thick);
int MORSE_dpotrf_diag_Tile(MORSE_enum uplo, MORSE_desc_t *A, int diag_thick);
int MORSE_dpotrf_diag_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, int diag_thick,
                             MORSE_sequence_t *sequence, MORSE_request_t *request);

#endif
