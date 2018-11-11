/**
 *
 * Copyright (c) 2017-2018  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file nc_file.h
 *
 * Heade file of auxiliary functions that are used to read and process flat files.
 *
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2018-11-11
 *
 **/
#ifndef _NC_FILE_H_
#define _NC_FILE_H_
#include <netcdf.h>
#include "MLE_misc.h"

int openFileNC(MLE_data* data, char *filename);
int countlinesNC(int ncid, char *dim1, char *dim2);
void readLocsNC_1d(MLE_data *data, int ncid);
void readLocsNC_2d(MLE_data *data, int ncid);
void readVarNCs(MLE_data *data, int ncid, char *varname, double *data_in, char *dim1, char *dim2);
void closeFileNC(MLE_data* data, int ncid);

#endif
