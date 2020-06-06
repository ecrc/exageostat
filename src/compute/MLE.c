/**
 *
 * Copyright (c) 2017-2020, King Abdullah University of Science and Technology
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
 * @version 1.1.0
 *
 * @author Sameh Abdulah
 * @date 2020-06-06
 *
 **/
#include "../include/MLE.h"
//***************************************************************************************
/** ****************************************************************************
 *  Current software version.
 **/
const char *argp_program_version = "Version 1.1.0";


void MLE_sdregister(MLE_data *data)
{
    int precision           = data->precision;
    if(precision == 2)      //single/double precision.
        MORSE_MLE_sdregister_Tile(data);
    else
        printf("MLE_sdregister is only needed with mixed precision mode (no actioin required)\n");
}

void MLE_zvg(MLE_data *data,  double * Nrand, double * initial_theta, int n, int ts, int log, int p_grid, int q_grid)
    //! Generate Z observation vector using generated or given locations
    /*!  -- using only dense (single, double, or mixed) precision.
     * Returns Z observation vector
     * @param[in] data: MLE_data struct with different MLE inputs.
     * @param[in] Nrand: A uniform random vector with size n that is used to generate Z .
     * @param[in] initial_theta: Theta vector with three parameter (Variance, Range, Smoothness)
     *                           that is used to to generate the Covariance Matrix.
     * @param[in] n: Problem size (number spatial locations).
     * @param[in] ts: tile size (MB) is used only in the case of HiCMA not MORSE.
     * @param[in] log: equals one if the user needs to generate log files for his problem.
     * @param[in] p_grid: p_grid in the case of distributed system.
     * @param[in] q_grid: q_grid in the case of distributed system.
     * */
{
    data->iter_count	= 0;
    int async		= data->async;
    char *computation       = data->computation;
    int precision		= data->precision;

    //	printf("Precision: %d,\n", precision);

    if(precision == 0)      //double precision.
    {
        if ( strcmp (computation, "exact") == 0 && async == 0)
            MORSE_MLE_dzvg_Tile(data, Nrand, initial_theta, n, ts, log);
        else if ( strcmp (computation, "exact") == 0 && async == 1)
            MORSE_MLE_dzvg_Tile_Async(data, Nrand, initial_theta, n, ts, log);
        else if ( strcmp (computation, "diag_approx") == 0 && async == 0)
            MORSE_MLE_dzvg_Tile(data, Nrand, initial_theta, n, ts, log);
        else if ( strcmp (computation, "diag_approx") == 0 && async == 1)
            MORSE_MLE_dzvg_Tile_Async(data, Nrand, initial_theta, n, ts, log);    
#if defined( EXAGEOSTAT_USE_HICMA )
        else if ( strcmp (computation, "lr_approx") == 0 && async == 0)
            HICMA_MLE_dzvg_Tile(data, Nrand, initial_theta, n, ts, log, p_grid, q_grid);
        else if ( strcmp (computation, "lr_approx") == 0 && async == 1)
            HICMA_MLE_dzvg_Tile_Async(data, Nrand, initial_theta, n, ts, log, p_grid, q_grid);    
#endif
    }

    else if(precision == 1)    //single precision.
    {
        if ( strcmp (computation, "exact") == 0 && async == 0)
            MORSE_MLE_szvg_Tile(data, Nrand, initial_theta, n, ts, log);
        else if ( strcmp (computation, "exact") == 0 && async == 1)
            MORSE_MLE_szvg_Tile_Async(data, Nrand, initial_theta, n, ts, log);
        else if ( strcmp (computation, "diag_approx") == 0 && async == 0)
            MORSE_MLE_szvg_Tile(data, Nrand, initial_theta, n, ts, log);
        else if ( strcmp (computation, "diag_approx") == 0 && async == 1)
            MORSE_MLE_szvg_Tile_Async(data, Nrand, initial_theta, n, ts, log);
#if defined( EXAGEOSTAT_USE_HICMA )
        //else if ( strcmp (computation, "lr_approx") == 0 && async == 0)
        //        HICMA_MLE_szvg_Tile(data, Nrand, initial_theta, n, ts, log, p_grid, q_grid);
        //else if ( strcmp (computation, "lr_approx") == 0 && async == 1)
        //        HICMA_MLE_szvg_Tile_Async(data, Nrand, initial_theta, n, ts, log, p_grid, q_grid);
#endif
    }

    else if(precision == 2)     //Mixed-precision (single/double).
    {
        if ( strcmp (computation, "exact") == 0 && async == 0)
            MORSE_MLE_sdzvg_Tile(data, Nrand, initial_theta, n, ts, log);
        else if ( strcmp (computation, "exact") == 0 && async == 1)
            MORSE_MLE_sdzvg_Tile_Async(data, Nrand, initial_theta, n, ts, log);
        //else if ( strcmp (computation, "diag_approx") == 0 && async == 0)
        //MORSE_MLE_szvg_Tile(data, Nrand, initial_theta, n, ts, log);
        //else if ( strcmp (computation, "diag_approx") == 0 && async == 1)
        //MORSE_MLE_szvg_Tile_Async(data, Nrand, initial_theta, n, ts, log);
#if defined( EXAGEOSTAT_USE_HICMA )
        //else if ( strcmp (computation, "lr_approx") == 0 && async == 0)
        //        HICMA_MLE_szvg_Tile(data, Nrand, initial_theta, n, ts, log, p_grid, q_grid);
        //else if ( strcmp (computation, "lr_approx") == 0 && async == 1)
        //        HICMA_MLE_szvg_Tile_Async(data, Nrand, initial_theta, n, ts, log, p_grid, q_grid);
#endif
    }

    //else if(strcmp (precision, "mixed") == 0 && async == 0)
    //MORSE_MLE_mzvg_Tile(data, Nrand, initial_theta, n, ts, log);
    //else if(strcmp (precision, "mixed") == 0 && async == 1)
    //MORSE_MLE_mzvg_Tile_Async(data, Nrand, initial_theta, n, ts, log);

}


//void MLE_zvr(MLE_data *data, int n, char *format)
//! Read Z observation vector using a give file.
/*!  -- using dense or approximate computation
 * Returns Z observation vector
 * @param[in] data: MLE_data struct with different MLE inputs.
 * @param[in] n: Problem size (number spatial locations).
 * */
double MLE_alg(unsigned n, const double * theta, double * grad, void * data)
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
    char *computation	= ((MLE_data*)data)->computation;
    int async		= ((MLE_data*)data)->async;
    int precision         = ((MLE_data*)data)->precision;

    //printf("Precision: %d,\n", precision);
    // printf("%s - %s\n", computation, __func__);

    if (precision == 0 ){
        if(strcmp (computation, "exact") == 0)
        {
            if(async == 0)
                return MORSE_dmle_Tile(n, theta,  grad,  data);
            else
                return MORSE_dmle_Tile_Async(n, theta,  grad,  data);
        }
        else if (strcmp (computation, "diag_approx") == 0)
        {	
            if(async == 0)
                return MORSE_dmle_diag_Tile(n, theta,  grad,  data); 
            else
                return MORSE_dmle_diag_Tile_Async(n, theta,  grad,  data);
        }
#if defined(EXAGEOSTAT_USE_HICMA)
        else if (strcmp (computation, "lr_approx") == 0)				
        {
            if(async == 0)
                return HICMA_dmle_Tile(n, theta,  grad,  data);
            else
                return HICMA_dmle_Tile_Async(n, theta,  grad,  data);
        }
#endif
    }
    else if (precision == 1)    //single precision.
    {    
        if(strcmp (computation, "exact") == 0)
        {
            if(async == 0)
                return MORSE_smle_Tile(n, theta,  grad,  data);
            else
                return MORSE_smle_Tile_Async(n, theta,  grad,  data);
        }
    }
    //else if (strcmp (computation, "diag_approx") == 0)
    //if(async == 0)
    //return MORSE_smle_diag_Tile(n, theta,  grad,  data);
    //else
    //return MORSE_smle_diag_Tile_Async(n, theta,  grad,  data);
    //#if defined(EXAGEOSTAT_USE_HICMA)
    //else if (strcmp (computation, "lr_approx") == 0)
    //if(async == 0)
    //return HICMA_smle_Tile(n, theta,  grad,  data);
    //else
    //return HICMA_smle_Tile_Async(n, theta,  grad,  data);
    //#endif

    else if (precision == 2)   //mixed-precision (single/double).
    { 
        if(strcmp (computation, "exact") == 0 && async == 0)
            return MORSE_sdmle_Tile(n, theta,  grad,  data);
        else if(strcmp (computation, "exact") == 0 && async == 1)
            return MORSE_sdmle_Tile_Async(n, theta,  grad,  data);
    }
    // else if (strcmp (computation, "diag_approx") == 0)
    // if(async == 0)
    // return MORSE_smle_diag_Tile(n, theta,  grad,  data);
    // else
    // return MORSE_smle_diag_Tile_Async(n, theta,  grad,  data);
    //#if defined(EXAGEOSTAT_USE_HICMA)
    //else if (strcmp (computation, "lr_approx") == 0)
    //if(async == 0)
    //return HICMA_smle_Tile(n, theta,  grad,  data);
    //else
    //return HICMA_smle_Tile_Async(n, theta,  grad,  data);
    //#endif
    return 0;
}

void exageostat_init(int *ncores, int *gpus, int *dts, int *lts)
//! initialize exageostat (initiate underlying library)
    /* @param[in] ncores: number of used CPU cores.
     * @param[in] gpus: number of used GPU units.	
     * @param[in] ts: tile size.
     */
{
    MORSE_context_t *morse;
    morse = morse_context_self();
    if (morse != NULL) {
        printf("Another instance is already running...!");
    }
    else{

        MORSE_user_tag_size(31, 26);
        MORSE_Init(*ncores, *gpus);
    }
    // MORSE_Enable(MORSE_AUTOTUNING);
    //MORSE_Set(MORSE_TILE_SIZE, *dts);
    //MORSE_Set(HICMA_TILE_SIZE, *lts); //should be added to HiCMA
}

void exageostat_finalize()
//! finalize exageostat (initiate underlying library)
{
    MORSE_context_t *morse;
    morse = morse_context_self();
    if (morse == NULL) {
        printf("No active instance...please use exageostat_init() function to initiate a new instance!\n");
    }
    else
        MORSE_Finalize();
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

    if ( data->precision == 0)
    {        	if(strcmp(data->computation,"exact") == 0 || strcmp(data->computation,"diag_approx") == 0 || strcmp(data->computation,"lr_approx") == 0)
        MORSE_dmle_Predict_Allocate(data, nZmiss, nZobs, ts, p_grid, q_grid, mse_flag);
        //#if defined(EXAGEOSTAT_USE_HICMA)
        //	else if(strcmp(data->computation,"lr_approx") == 0)
        //		HICMA_dmle_Predict_Allocate(data, nZmiss, nZobs, ts, p_grid, q_grid, mse_flag);		
        //#endif
    }

    else if (data->precision == 1)
    {                if(strcmp(data->computation,"exact") == 0 || strcmp(data->computation,"diag_approx") == 0)
        MORSE_smle_Predict_Allocate(data, nZmiss, nZobs, ts, p_grid, q_grid, mse_flag);
    }


    else if (data->precision == 2)
    {                if(strcmp(data->computation,"exact") == 0 || strcmp(data->computation,"diag_approx") == 0)

        MORSE_sdmle_Predict_Allocate(data, nZmiss, nZobs, ts, p_grid, q_grid, mse_flag);

    }

    /*               #if defined(EXAGEOSTAT_USE_HICMA)
                     else if(strcmp(data->computation,"lr_approx") == 0)
                     HICMA_smle_Predict_Allocate(data, nZmiss, nZobs, ts, p_grid, q_grid, mse_flag);
#endif
*/
    /*
       else if (strcmp (data->precision, "mixed") == 0)
       if(strcmp(data->computation,"exact") == 0 || strcmp(data->computation,"diag_approx") == 0)
       MORSE_mmle_Predict_Allocate(data, nZmiss, nZobs, ts, p_grid, q_grid, mse_flag);
#if defined(EXAGEOSTAT_USE_HICMA)
else if(strcmp(data->computation,"lr_approx") == 0)
HICMA_mmle_Predict_Allocate(data, nZmiss, nZobs, ts, p_grid, q_grid, mse_flag);
#endif
*/
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

    //	if ( data->precision == 0)
    //	{
    //		if(strcmp(data->computation,"exact") == 0 || strcmp(data->computation,"diag_approx") == 0)
    MORSE_dmle_mloe_mmom_Allocate(data, nZmiss, nZobs, ts, p_grid, q_grid);
    //	}
    /*
#if defined(EXAGEOSTAT_USE_HICMA)
else if(strcmp(data->computation,"lr_approx") == 0)
HICMA_dmle_Predict_Allocate(data, nZmiss, nZobs, ts, p_grid, q_grid, mse_flag);
#endif
}

else if (data->precision == 1)
{                if(strcmp(data->computation,"exact") == 0 || strcmp(data->computation,"diag_approx") == 0)
MORSE_smle_Predict_Allocate(data, nZmiss, nZobs, ts, p_grid, q_grid, mse_flag);
}


else if (data->precision == 2)
{                if(strcmp(data->computation,"exact") == 0 || strcmp(data->computation,"diag_approx") == 0)

MORSE_sdmle_Predict_Allocate(data, nZmiss, nZobs, ts, p_grid, q_grid, mse_flag);

}
*/
    }


void mloe_mmom_finalize(MLE_data *data)
    //! Destory and free prediction memory
    /*! allocations.
     * @param[in] data: MLE_data struct with different MLE inputs.
     */
{
    if(strcmp(data->computation,"exact") == 0 || strcmp(data->computation,"diag_approx") == 0)
    {
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->desck_t) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->desck_a) );
        //        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descK_atmp) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descK_t) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descK_a) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descexpr1) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descexpr2) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descexpr3) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descexpr4) );
    }
    /*
#if defined(EXAGEOSTAT_USE_HICMA)
else if(strcmp(data->computation,"lr_approx") == 0)
{
MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descZobs) );
MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descZactual) );
MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descZmiss) );
MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->hicma_descC12D) );
MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->hicma_descC12UV) );
MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->hicma_descC12rk) );
MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->hicma_descC22D) );
MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->hicma_descC22UV) );
MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->hicma_descC22rk) );
}
#endif
*/
}

void prediction_finalize(MLE_data *data)
    //! Destory and free prediction memory
    /*! allocations.
     * @param[in] data: MLE_data struct with different MLE inputs.
     */
{
    if(strcmp(data->computation,"exact") == 0 || strcmp(data->computation,"diag_approx") == 0)
    {
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descZobs) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descZactual) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descZmiss) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descC12) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descC22) );


    }
#if defined(EXAGEOSTAT_USE_HICMA)
    else if(strcmp(data->computation,"lr_approx") == 0)
    {
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descZobs) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descZactual) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descZmiss) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->hicma_descC12D) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->hicma_descC12UV) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->hicma_descC12rk) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->hicma_descC22D) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->hicma_descC22UV) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->hicma_descC22rk) );
    }
#endif
}

void MLE_Finalize(MLE_data *data)
    //! Destory and free MLE memory
    /*! allocations.
     * @param[in] data: MLE_data struct with different MLE inputs.
     */
{
    if(strcmp(data->computation,"exact") == 0 || strcmp(data->computation,"diag_approx") == 0)
    {
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descC) );
        //MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descZ) );
        //MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descZcpy) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descproduct) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descdet) );
    }
#if defined(EXAGEOSTAT_USE_HICMA)
    else if(strcmp(data->computation,"lr_approx") == 0)
    {
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->hicma_descCD) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->hicma_descCUV) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->hicma_descCrk) );
        //MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->hicma_descZ) );
        //MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->hicma_descZcpy) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->hicma_descproduct) );
        MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->hicma_descdet) );
    }
#endif
    //	MORSE_Finalize();
}


void MLOE_MMOM_Finalize(MLE_data *data)
{


    MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->desck_t) );
    MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->desck_a) );
    MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descK_ttmp) );
    MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->desck_atmp) );
    MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descK_t) );
    MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->desck_ttmp) );
    MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descK_a) );
    MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descexpr1) );
    MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descexpr2) );
    MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descexpr3) );
    MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descexpr4) );
    MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->desctruthalpha) );
    MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->descestimatedalpha) );
    //      MORSE_Desc_Destroy( (MORSE_desc_t **)&(data->desc_mloe_mmom) );
}


void MLE_get_zobs(MLE_data *data, double *z, int n)
    //! Convert measurements vector from Lapack 
    /*! represenation to Chameleon representation.
     * @param[in] data: MLE_data struct with different MLE inputs.
     * @param[in] z: Lapack measurements vector.
     * @param[in] n: number of meaurements.
     */
{
    if(strcmp(data->computation,"exact") == 0 || strcmp(data->computation,"diag_approx") == 0)
        MORSE_Tile_to_Lapack( ((MLE_data*)data)->descZcpy, z, n);
#if defined(EXAGEOSTAT_USE_HICMA)
    else if(strcmp(data->computation,"lr_approx") == 0)
        MORSE_Tile_to_Lapack( ((MLE_data*)data)->hicma_descZcpy, z, n);
#endif
}
