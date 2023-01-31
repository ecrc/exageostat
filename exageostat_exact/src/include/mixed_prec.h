/**
 *
 * Copyright (c) 2017-2023  King Abdullah University of Science and Technology
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
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#ifndef _MIXED_PREC_H_
#define _MIXED_PREC_H_

#include <starpu_exageostat.h>
#include <common.h>


int EXAGEOSTAT_sdpotrf_Tile(CHAM_enum uplo, CHAM_desc_t *A, CHAM_desc_t *temp, int diag_thick);

int EXAGEOSTAT_sdpotrf_Tile_Async(CHAM_enum uplo, CHAM_desc_t *A, CHAM_desc_t *temp,
                                  int diag_thick, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request);

#endif
