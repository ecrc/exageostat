/**
 *
 * Copyright (c) 2017-2020  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file exageostat-constants.h
 *
 * Auxiliary functions that are needed by ExaGeoStat.
 *
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2020-05-02
 *
 **/
#ifndef __EXAGEOSTAT_CONSTANTS_H__
#define __EXAGEOSTAT_CONSTANTS_H__

//! Enum for kernel types
enum EXAGEOSTAT_KERNEL
{
	STARSH_BACKEND_NOTSELECTED = -2,
	//!< Kernel has not been yet selected.
	STARSH_BACKEND_NOTSUPPORTED = -1,
	//!< Error, kernel is not supported.
	UNIVARIATE_MATERN_STATIONARY = 0,
        //!< Univariate matern stationary kernel.
	UNIVARIATE_POW_EXP_STATIONARY = 1,
        //!< Univariate powered exponential stationary kernel.
	UNIVARIATE_SQR_EXP_STATIONARY = 2,
        //!< Univariate square exponential stationary kernel.
	UNIVARIATE_MATERN_NON_STATIONARY = 3,
        //!< Univariate matern nonstationary kernel.
	BIVARIATE_MATERN_PARSIMONIOUS = 4,
        //!< Bivariate matern parsimonious kernel.
	BIVARIATE_MATERN_PARSIMONIOUS2 = 5,
        //!< Bivariate matern parsimonious2 kernel.
	BIVARIATE_MATERN_PARSIMONIOUS_PROFILE = 6,
        //!< Bivariate matern parsimonious kernel (profile).
	BIVARIATE_MATERN_PARSIMONIOUS2_PROFILE = 7,
        //!< Bivariate matern parsimonious2 kernel (profile).
	BIVARIATE_MATERN_FLEXIBLE = 8,
        //!< Bivariate matern flexiable kernel.
	BIVARIATE_MATERN_FLEXIBLE_PROFILE = 9
        //!< Bivariate matern flexiable kernel (profile).
};

#endif // __EXAGEOSTAT_CONSTANTS_H__
