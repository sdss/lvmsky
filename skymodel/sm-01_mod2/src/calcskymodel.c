/*
 *  This file is part of the Sky Background software package.
 *  Copyright (C) 2009-2018 European Southern Observatory
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

/*!
 * \ingroup sm_module2_run
 */

/**@{*/

/*!
 * \file calcskymodel.c
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  19 Oct 2009
 * \date   06 Oct 2015
 *
 * Stand-alone programme for ETC sky model
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <sm_general.h>
#include <sm_skyemcomp.h>
#include <sm_run.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

int main(void)
{
    /*!
     * Main routine for running the ETC sky model as stand-alone programmme.
     * It has to be started from the base directory of the sky model code by
     * typing "bin/calcskymodel". It requires the parameter files
     * "skymodel_etc.par" and "instrument_etc.par" in the "config/" folder and
     * returns a radiance spectrum and transmission curve in FITS format in
     * the folder "output/". For testing purposes, the run time of the code
     * and all errors and warnings are provided on stdout.
     */

    char configdir[FILENAME_MAX] = "config/";
    char modelparfile[FILENAME_MAX] = "skymodel_etc.par";
    char instruparfile[FILENAME_MAX] = "instrument_etc.par";
    char datadir[FILENAME_MAX] = "data/";
    char outputdir[FILENAME_MAX] = "output/";
    char radfitsfile[FILENAME_MAX] = "radspec.fits";
    char transfitsfile[FILENAME_MAX] = "transspec.fits";

    cpl_errorstate errstate = cpl_errorstate_get();
    cpl_parameterlist *params, *lamgrid;
    cpl_table *skytable;
    char modelpar[FILENAME_MAX], instrupar[FILENAME_MAX];
    char radfits[FILENAME_MAX], transfits[FILENAME_MAX];
    double cs, ce;

    /* Initialise CPL */
    cpl_init(CPL_INIT_DEFAULT);

    /* Check existence of config/ folder */
    if (sm_basic_access(configdir, F_OK) != CPL_ERROR_NONE) {
        cpl_errorstate_dump(errstate, CPL_FALSE, cpl_errorstate_dump_one);
        cpl_end();
        return EXIT_FAILURE;
    }

    /* Build full names of parameter files */
    sprintf(modelpar, "%s%s", configdir, modelparfile);
    sprintf(instrupar, "%s%s", configdir, instruparfile);

    /* Check existence of model parameter file */
    if (sm_basic_access(modelpar, F_OK) != CPL_ERROR_NONE) {
        cpl_errorstate_dump(errstate, CPL_FALSE, cpl_errorstate_dump_one);
        cpl_end();
        return EXIT_FAILURE;
    }

    /* Check existence of instrument parameter file */
    if (sm_basic_access(instrupar, F_OK) != CPL_ERROR_NONE) {
        cpl_errorstate_dump(errstate, CPL_FALSE, cpl_errorstate_dump_one);
        cpl_end();
        return EXIT_FAILURE;
    }

    /* Check existence of data/ folder */
    if (sm_basic_access(datadir, F_OK) != CPL_ERROR_NONE) {
        cpl_errorstate_dump(errstate, CPL_FALSE, cpl_errorstate_dump_one);
        cpl_end();
        return EXIT_FAILURE;
    }

    /* Check existence of output/ folder */
    if (sm_basic_access(outputdir, F_OK) != CPL_ERROR_NONE) {
        cpl_errorstate_dump(errstate, CPL_FALSE, cpl_errorstate_dump_one);
        cpl_end();
        return EXIT_FAILURE;
    }

    /* Build full names of output files */
    sprintf(radfits, "%s%s", outputdir, radfitsfile);
    sprintf(transfits, "%s%s", outputdir, transfitsfile);

    /* Create required CPL objects */
    params = cpl_parameterlist_new();
    lamgrid = cpl_parameterlist_new();
    skytable = cpl_table_new(0);

    /* Read sky model parameter file */
    sm_run_readmodelpar(params, modelpar);

    /* Read instrument parameter file */
    sm_run_readinstrupar(lamgrid, params, instrupar);

    /* Create CPL table for the input wavelength grid */
    sm_run_createtable(skytable, lamgrid);

    /* Check for errors */
    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        cpl_parameterlist_delete(params);
        cpl_parameterlist_delete(lamgrid);
        cpl_table_delete(skytable);
        cpl_errorstate_dump(errstate, CPL_FALSE, cpl_errorstate_dump_one);
        cpl_end();
        return EXIT_FAILURE;
    }

    /* Calculate sky brightness and transmission curve and perform
       convolution */
    cs = cpl_test_get_cputime();
    sm_run_calcmodel(skytable, params);
    ce = cpl_test_get_cputime();
    cpl_msg_info(cpl_func, "Run time: %g s", ce - cs);

    /* Write resulting spectra to FITS file */
    if (cpl_error_get_code() == CPL_ERROR_NONE) {
        sm_run_writefits(skytable, radfits, transfits);
    } else {
        cpl_msg_warning(cpl_func, "No output files due to error(s)!");
    }

    /* Free allocated memory */
    cpl_parameterlist_delete(params);
    cpl_parameterlist_delete(lamgrid);
    cpl_table_delete(skytable);

    /* List errors */
    if (cpl_errorstate_is_equal(errstate)) {
        cpl_msg_info(cpl_func, "No errors occurred");
    } else {
        cpl_errorstate_dump(errstate, CPL_TRUE, cpl_errorstate_dump_one);
    }

    /* Finish CPL */
    cpl_end();

    return EXIT_SUCCESS;
}

/**@}*/
