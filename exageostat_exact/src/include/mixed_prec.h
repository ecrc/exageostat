/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file mixed_pre.h
 *
 * ExaGeoStat approx computation main functions header file.
 *
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2018-11-11
 *
 **/
#ifndef _MIXED_PREC_H_
#define _MIXED_PREC_H_
#include "starpu_exageostat.h"
#include "common.h"


int MORSE_sdpotrf_Tile(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *temp,
		int diag_thick);

int MORSE_sdpotrf_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *temp, 
		int diag_thick, MORSE_sequence_t *sequence, MORSE_request_t *request);

#endif
