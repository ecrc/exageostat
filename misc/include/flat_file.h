/**
 *
 * Copyright (c) 2017-2018  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file flat_file.h
 *
 * Heade file of auxiliary functions that are used to read and process flat files.
 *
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2018-11-11
 *
 **/
#ifndef _FLAT_FILE_H_
#define _FLAT_FILE_H_
#include "MLE_misc.h"

location* readLocsFile(char *locs_file, int n);
double * readObsFile(char *obsfile, int n);

#endif
