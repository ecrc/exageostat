/**
 *
 * Copyright (c) 2017-2023, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file MLE.c
 *
 * ExaGeoStat main functions.
 *
 * @version 1.2.0
 *
 * @author Sameh Abdulah
 * @date 2022-11-09
 *
 **/
#include "../include/MLE.h"
#if defined(EXAGEOSTAT_USE_HICMA)
	#include "hicma/hicma_ext/control/hicma_context.h"
#endif
/** ****************************************************************************
 *  Current software version.
 **/
const char *argp_program_version = "Version 1.1.0";


void MLE_sdregister(MLE_data *data) {
    int precision = data->precision;
    if (precision == 2)      //single/double precision.
        EXAGEOSTAT_MLE_sdregister_Tile(data);
    else
        printf("MLE_sdregister is only needed with mixed precision mode (no actioin required)\n");
}

void MLE_zvg(MLE_data *data, double* Nrand, double* initial_theta, int n, int ts, int log, int p_grid, int q_grid)
//! Generate Z observation vector using generated or given locations
/*!  -- using only dense (single, double, or mixed) precision.
 * Returns Z observation vector
 * @param[in] data: MLE_data struct with different MLE inputs.
 * @param[in] Nrand: A uniform random vector with size n that is used to generate Z .
 * @param[in] initial_theta: Theta vector with three parameter (Variance, Range, Smoothness)
 *                           that is used to to generate the Covariance Matrix.
 * @param[in] n: Problem size (number spatial locations).
 * @param[in] ts: tile size (MB) is used only in the case of HiCMA not CHAMELEON.
 * @param[in] log: equals one if the user needs to generate log files for his problem.
 * @param[in] p_grid: p_grid in the case of distributed system.
 * @param[in] q_grid: q_grid in the case of distributed system.
 * */
{
    data->iter_count = 0;
    int async = data->async;
    char *computation = data->computation;
    int precision = data->precision;

    if (precision == 0)      //double precision.
    {
        if (strcmp(computation, "exact") == 0 && async == 0)
            EXAGEOSTAT_MLE_dzvg_Tile(data, Nrand, initial_theta, n, ts, log);
        else if (strcmp(computation, "exact") == 0 && async == 1)
            EXAGEOSTAT_MLE_dzvg_Tile_Async(data, Nrand, initial_theta, n, ts, log);
        else if (strcmp(computation, "diag_approx") == 0 && async == 0)
            EXAGEOSTAT_MLE_dzvg_Tile(data, Nrand, initial_theta, n, ts, log);
        else if (strcmp(computation, "diag_approx") == 0 && async == 1)
            EXAGEOSTAT_MLE_dzvg_Tile_Async(data, Nrand, initial_theta, n, ts, log);
#if defined( EXAGEOSTAT_USE_HICMA )
        else if (strcmp(computation, "lr_approx") == 0 && async == 0){
            EXAGEOSTAT_TLR_MLE_dzvg_Tile(data, Nrand, initial_theta, n, ts, log, p_grid, q_grid);
        }
        else if (strcmp(computation, "lr_approx") == 0 && async == 1)
            EXAGEOSTAT_TLR_MLE_dzvg_Tile_Async(data, Nrand, initial_theta, n, ts, log, p_grid, q_grid);
#endif
    } else if (precision == 1)    //single precision.
    {
        if (strcmp(computation, "exact") == 0 && async == 0)
            EXAGEOSTAT_MLE_szvg_Tile(data, Nrand, initial_theta, n, ts, log);
        else if (strcmp(computation, "exact") == 0 && async == 1)
            EXAGEOSTAT_MLE_szvg_Tile_Async(data, Nrand, initial_theta, n, ts, log);
        else if (strcmp(computation, "diag_approx") == 0 && async == 0)
            EXAGEOSTAT_MLE_szvg_Tile(data, Nrand, initial_theta, n, ts, log);
        else if (strcmp(computation, "diag_approx") == 0 && async == 1)
            EXAGEOSTAT_MLE_szvg_Tile_Async(data, Nrand, initial_theta, n, ts, log);
    } else if (precision == 2)     //Mixed-precision (single/double).
    {
        if (strcmp(computation, "exact") == 0 && async == 0)
            EXAGEOSTAT_MLE_sdzvg_Tile(data, Nrand, initial_theta, n, ts, log);
        else if (strcmp(computation, "exact") == 0 && async == 1)
            EXAGEOSTAT_MLE_sdzvg_Tile_Async(data, Nrand, initial_theta, n, ts, log);
    }
}

void MLE_ng_zvg(MLE_data *data, double* Nrand, double* initial_theta, int n, int ts, int log, int p_grid, int q_grid)
//! Generate Z observation vector using generated or given locations
/*!  -- using only dense (single, double, or mixed) precision.
 * Returns Z observation vector
 * @param[in] data: MLE_data struct with different MLE inputs.
 * @param[in] Nrand: A uniform random vector with size n that is used to generate Z .
 * @param[in] initial_theta: Theta vector with three parameter (Variance, Range, Smoothness)
 *                           that is used to to generate the Covariance Matrix.
 * @param[in] n: Problem size (number spatial locations).
 * @param[in] ts: tile size (MB) is used only in the case of HiCMA not CHAMELEON.
 * @param[in] log: equals one if the user needs to generate log files for his problem.
 * @param[in] p_grid: p_grid in the case of distributed system.
 * @param[in] q_grid: q_grid in the case of distributed system.
 * */
{
    data->iter_count = 0;
    int async = data->async;
    char *computation = data->computation;
    int precision = data->precision;

    if (precision == 0)      //double precision.
    {
        if (strcmp(computation, "exact") == 0 && async == 0)
            EXAGEOSTAT_MLE_ng_dzvg_Tile(data, Nrand, initial_theta, n, ts, log);

#if defined( EXAGEOSTAT_USE_HICMA )
        else if (strcmp(computation, "lr_approx") == 0 && async == 0)
            EXAGEOSTAT_TLR_MLE_dzvg_ng_Tile(data, Nrand, initial_theta, n, ts, log, p_grid, q_grid);
#endif
    } else if (precision == 1)    //single precision.
    {
        //// TODO: Implement
    } else if (precision == 2)     //Mixed-precision (single/double).
    {
        //// TODO: Implement
    }
}

//! Read Z observation vector using a give file.
/*!  -- using dense or approximate computation
 * Returns Z observation vector
 * @param[in] data: MLE_data struct with different MLE inputs.
 * @param[in] n: Problem size (number spatial locations).
 * */
double MLE_alg(unsigned n, const double* theta, double* grad, void *data)
//! Maximum Likelihood Evaluation (MLE)
/*!  -- using exact or approximation computation, and (single, double, or mixed) precision.
 * Returns the loglikelihhod value for the given theta.
 * @param[in] n: unsigned variable used by NLOPT package.
 * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
 *                           that is used to to generate the Covariance Matrix.
 * @param[in] grad: double variable used by NLOPT package.
 * @param[in] data: MLE_data struct with different MLE inputs.
 */
{
    char *computation = ((MLE_data *) data)->computation;
    int async = ((MLE_data *) data)->async;
    int precision = ((MLE_data *) data)->precision;

    if (precision == 0) {
        if (strcmp(computation, "exact") == 0) {
            if (async == 0)
                return EXAGEOSTAT_dmle_Tile(n, theta, grad, data);
            else
                return EXAGEOSTAT_dmle_Tile_Async(n, theta, grad, data);
        } else if (strcmp(computation, "diag_approx") == 0) {
            if (async == 0)
                return EXAGEOSTAT_dmle_diag_Tile(n, theta, grad, data);
            else
                return EXAGEOSTAT_dmle_diag_Tile_Async(n, theta, grad, data);
        }
#if defined(EXAGEOSTAT_USE_HICMA)
        else if (strcmp(computation, "lr_approx") == 0) {
            if (async == 0){
                return EXAGEOSTAT_TLR_dmle_Tile(n, theta, grad, data);
            }
            else{
                return EXAGEOSTAT_TLR_dmle_Tile_Async(n, theta, grad, data);
            }
        }
#endif
    } else if (precision == 1)    //single precision.
    {
        if (strcmp(computation, "exact") == 0) {
            if (async == 0) {
                return EXAGEOSTAT_smle_Tile(n, theta, grad, data);
            } else {
                return EXAGEOSTAT_smle_Tile_Async(n, theta, grad, data);
            }
        }
    } else if (precision == 2)   //mixed-precision (single/double).
    {

        if (strcmp(computation, "exact") == 0 && async == 0)
            return EXAGEOSTAT_sdmle_Tile(n, theta, grad, data);
        else if (strcmp(computation, "exact") == 0 && async == 1)
            return EXAGEOSTAT_sdmle_Tile_Async(n, theta, grad, data);
    }
    return 0;
}

double MLE_ng_alg(unsigned n, const double* theta, double* grad, void *data)
//! Maximum Likelihood Evaluation (MLE)
/*!  -- using exact or approximation computation, and (single, double, or mixed) precision.
 * Returns the loglikelihhod value for the given theta.
 * @param[in] n: unsigned variable used by NLOPT package.
 * @param[in] theta: theta Vector with three parameter (Variance, Range, Smoothness)
 *                           that is used to to generate the Covariance Matrix.
 * @param[in] grad: double variable used by NLOPT package.
 * @param[in] data: MLE_data struct with different MLE inputs.
 */
{
    char *computation = ((MLE_data *) data)->computation;
    int async = ((MLE_data *) data)->async;
    int precision = ((MLE_data *) data)->precision;

    if (precision == 0) {
        if (strcmp(computation, "exact") == 0) {
            if (async == 0)
                return EXAGEOSTAT_dmle_ng_Tile(n, theta, grad, data);
        }
	else if (strcmp(computation, "diag_approx") == 0) {
            if (async == 0)
                //return EXAGEOSTAT_dmle_diag_ng_Tile(n, theta, grad, data);  // TODO
       		return 0;
       	}
#if defined(EXAGEOSTAT_USE_HICMA)
        else if (strcmp(computation, "lr_approx") == 0) {
            if (async == 0)
                return EXAGEOSTAT_TLR_dmle_ng_Tile(n, theta, grad, data);
        }
#endif
    } else if (precision == 2)   //mixed-precision (single/double).
    {
        //// TODO: Implement
    }
    return 0;
}

void exageostat_init(int *ncores, int *gpus, int *dts, int *lts)
//! initialize exageostat (initiate underlying library)
/* @param[in] ncores: number of used CPU cores.
 * @param[in] gpus: number of used GPU units.
 * @param[in] ts: tile size.
 */
{
    CHAM_context_t *chamctxt;
    chamctxt = chameleon_context_self();
    if (chamctxt != NULL) {
        printf("Another instance of Chameleon is already running...!");
    } else {
        CHAMELEON_user_tag_size(31, 26);
        CHAMELEON_Init(*ncores, *gpus);
    }


#if defined(EXAGEOSTAT_USE_HICMA)
    HICMA_context_t *hicmatxt;
    hicmatxt = hicma_context_self();
    if (hicmatxt != NULL) {
        printf("Another instance of HiCMA is already running...!");
    } else {
        HICMA_user_tag_size(31, 26);
        HICMA_Init(*ncores, *gpus);
    }
#endif    
}


void exageostat_finalize()
//! finalize exageostat (initiate underlying library)
{
    CHAM_context_t *chamctxt;
    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        printf("No active instance oh Chameleon...please use exageostat_init() function to initiate a new instance!\n");
    } else
        CHAMELEON_Finalize();

#if defined(EXAGEOSTAT_USE_HICMA)    
    HICMA_context_t *hicmatxt;
    hicmatxt = hicma_context_self();
    if (hicmatxt == NULL) {
        printf("No active instance of HICMA...please use exageostat_init() function to initiate a new instance!\n");
    } else
        HICMA_Finalize();
#endif
}


void prediction_init(MLE_data *data, int nZmiss, int nZobs, int ts, int p_grid, int q_grid, int mse_flag)
//! initialize exageostat prediction allocation (allocate memory).
/* @param[in] data: MLE_data struct with different MLE inputs.
 * @param[in] nZmiss: number of missing measurements.
 * @param[in] nZobs: number of observed measurements.
 * @param[in] p_grid: p_grid in the case of distributed system.
 * @param[in] q_grid: q_grid in the case of distributed system.
 * @param[in] mse_flag: flag to enable or disable Mean Square Error (MSE) computing.
 */
{

    if (data->precision == 0) {
        if (strcmp(data->computation, "exact") == 0 || strcmp(data->computation, "diag_approx") == 0 ||
            strcmp(data->computation, "lr_approx") == 0)
            EXAGEOSTAT_dmle_Predict_Allocate(data, nZmiss, nZobs, ts, p_grid, q_grid, mse_flag);
    } else if (data->precision == 1) {
        if (strcmp(data->computation, "exact") == 0 || strcmp(data->computation, "diag_approx") == 0)
            EXAGEOSTAT_smle_Predict_Allocate(data, nZmiss, nZobs, ts, p_grid, q_grid, mse_flag);
    } else if (data->precision == 2) {
        if (strcmp(data->computation, "exact") == 0 || strcmp(data->computation, "diag_approx") == 0)
            EXAGEOSTAT_sdmle_Predict_Allocate(data, nZmiss, nZobs, ts, p_grid, q_grid, mse_flag);
    }
}


void mloe_mmom_init(MLE_data *data, int nZmiss, int nZobs, int ts, int p_grid, int q_grid)
//! initialize exageostat prediction allocation (allocate memory).
/* @param[in] data: MLE_data struct with different MLE inputs.
 * @param[in] nZmiss: number of missing measurements.
 * @param[in] nZobs: number of observed measurements.
 * @param[in] p_grid: p_grid in the case of distributed system.
 * @param[in] q_grid: q_grid in the case of distributed system.
 * @param[in] mse_flag: flag to enable or disable Mean Square Error (MSE) computing.
 */
{
    EXAGEOSTAT_dmle_mloe_mmom_Allocate(data, nZmiss, nZobs, ts, p_grid, q_grid);
}


void mloe_mmom_finalize(MLE_data *data)
//! Destory and free prediction memory
/*! allocations.
 * @param[in] data: MLE_data struct with different MLE inputs.
 */
{
    if (strcmp(data->computation, "exact") == 0 || strcmp(data->computation, "diag_approx") == 0) {
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->desck_t));
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->desck_a));
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descK_t));
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descK_a));
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descexpr1));
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descexpr2));
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descexpr3));
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descexpr4));
    }
}

void prediction_finalize(MLE_data *data)
//! Destory and free prediction memory
/*! allocations.
 * @param[in] data: MLE_data struct with different MLE inputs.
 */
{

    if (strcmp(data->kernel_fun, "univariate_matern_non_gaussian") == 0
        || strcmp(data->kernel_fun, "univariate_exp_non_gaussian") == 0) {
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descZobs));
        if ((CHAM_desc_t **) &(data->descZactual) != NULL)
            CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descZactual));
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descZmiss));
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descr));
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descrcpy));
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descC22));
    } else {
        if (strcmp(data->computation, "exact") == 0 || strcmp(data->computation, "diag_approx") == 0) {
            CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descZobs));
            CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descZactual));
            CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descZmiss));
            CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descC12));
            CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descC22));
        }
#if defined(EXAGEOSTAT_USE_HICMA)
        else if (strcmp(data->computation, "lr_approx") == 0) {
            CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descZobs));
            CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descZactual));
            CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descZmiss));
            HICMA_Desc_Destroy((HICMA_desc_t **) &(data->hicma_descC12D));
            HICMA_Desc_Destroy((HICMA_desc_t **) &(data->hicma_descC12UV));
            HICMA_Desc_Destroy((HICMA_desc_t **) &(data->hicma_descC12rk));
            HICMA_Desc_Destroy((HICMA_desc_t **) &(data->hicma_descC22D));
            HICMA_Desc_Destroy((HICMA_desc_t **) &(data->hicma_descC22UV));
            HICMA_Desc_Destroy((HICMA_desc_t **) &(data->hicma_descC22rk));
        }
#endif
    }
}

void MLE_Finalize(MLE_data *data)
//! Destory and free MLE memory
/*! allocations.
 * @param[in] data: MLE_data struct with different MLE inputs.
 */
{
    if (strcmp(data->computation, "exact") == 0 || strcmp(data->computation, "diag_approx") == 0) {
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descC));
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descproduct));
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descdet));
    }
#if defined(EXAGEOSTAT_USE_HICMA)
    else if (strcmp(data->computation, "lr_approx") == 0) {
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->hicma_descCD));
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->hicma_descCUV));
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->hicma_descCrk));
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->hicma_descproduct));
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->hicma_descdet));
    }
#endif
}

void MLOE_MMOM_Finalize(MLE_data *data) {

    CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->desck_t));
    CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->desck_a));
    CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->desck_atmp));
    CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descK_t));
    CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->desck_ttmp));
    CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descK_a));
    CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descexpr1));
    CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descexpr2));
    CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descexpr3));
    CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descexpr4));
    CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->desctruthalpha));
    CHAMELEON_Desc_Destroy((CHAM_desc_t **) &(data->descestimatedalpha));
}


void MLE_get_zobs(MLE_data *data, double* z, int n)
//! Convert measurements vector from Lapack
/*! represenation to Chameleon representation.
 * @param[in] data: MLE_data struct with different MLE inputs.
 * @param[in] z: Lapack measurements vector.
 * @param[in] n: number of meaurements.
 */
{
    if (strcmp(data->computation, "exact") == 0 || strcmp(data->computation, "diag_approx") == 0)
        CHAMELEON_Tile_to_Lapack(((MLE_data *) data)->descZcpy, z, n);
#if defined(EXAGEOSTAT_USE_HICMA)
    else if (strcmp(data->computation, "lr_approx") == 0)
        HICMA_Tile_to_Lapack(((MLE_data *) data)->hicma_descZcpy, z, n);
#endif
}
