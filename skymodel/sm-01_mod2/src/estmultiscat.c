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
 * \file estmultiscat.c
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  21 Nov 2013
 * \date   06 Oct 2015
 *
 * Programme for estimating multiple scattering corrections for moonlight
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <sm_general.h>
#include <sm_scatmoonlight.h>
#include <sm_run.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

int main(void)
{
    /*!
     * Main routine for estimating multiple scattering corrections for
     * moonlight. It has to be started from the base directory of the sky
     * model code by typing "bin/estmultiscat". It requires the parameter
     * file "skymodel_etc.par" in the "config/" folder and returns the
     * multiple scattering correction ASCII files "sscatcor.dat" and
     * "dscatcor.dat" in the folder "output/". They are for correcting single
     * and double scattering calculations, respectively. For testing purposes,
     * the run time of the code and all errors and warnings are provided on
     * stdout.
     */

    char configdir[FILENAME_MAX] = "config/";
    char modelparname[FILENAME_MAX] = "skymodel_etc.par";
    char datadir[FILENAME_MAX] = "data/";
    char outputdir[FILENAME_MAX] = "output/";
    char sscatcorname[FILENAME_MAX] = "sscatcor.dat";
    char dscatcorname[FILENAME_MAX] = "dscatcor.dat";

    cpl_errorstate errstate = cpl_errorstate_get();
    cpl_parameterlist *params;
    smparmodel modelpar;
    char modelparfile[FILENAME_MAX], sscatcorfile[FILENAME_MAX];
    char dscatcorfile[FILENAME_MAX];
    double cs, ce;

    /* Initialise CPL */
    cpl_init(CPL_INIT_DEFAULT);

    /* Check existence of config/ folder */
    if (sm_basic_access(configdir, F_OK) != CPL_ERROR_NONE) {
        cpl_errorstate_dump(errstate, CPL_FALSE, cpl_errorstate_dump_one);
        cpl_end();
        return EXIT_FAILURE;
    }

    /* Build full name of sky model parameter file */
    sprintf(modelparfile, "%s%s", configdir, modelparname);

    /* Check existence of model parameter file */
    if (sm_basic_access(modelparfile, F_OK) != CPL_ERROR_NONE) {
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

    /* Create required CPL parameter list */
    params = cpl_parameterlist_new();

    /* Read sky model parameter file */
    if (sm_run_readmodelpar(params, modelparfile) != CPL_ERROR_NONE) {
        cpl_parameterlist_delete(params);
        cpl_errorstate_dump(errstate, CPL_FALSE, cpl_errorstate_dump_one);
        cpl_end();
        return EXIT_FAILURE;
    }

    /* Fill smparmodel structure with data of the input CPL parameter list
       and a special file containing crucial paths and file names */
    if (sm_etc_getparams(&modelpar, params) ||
        sm_etc_readfilenames(&modelpar, params)) {
        cpl_parameterlist_delete(params);
        cpl_errorstate_dump(errstate, CPL_FALSE, cpl_errorstate_dump_one);
        cpl_end();
        return EXIT_FAILURE;
    }

    /* Free allocated memory */
    cpl_parameterlist_delete(params);

    /* Build full names of output files */
    sprintf(sscatcorfile, "%s%s", outputdir, sscatcorname);
    sprintf(dscatcorfile, "%s%s", outputdir, dscatcorname);

    /* Calculate sky brightness and transmission curve and perform
       convolution */
    cs = cpl_test_get_cputime();
    sm_scat_estmultiscat(sscatcorfile, dscatcorfile, modelpar);
    ce = cpl_test_get_cputime();
    cpl_msg_info(cpl_func, "Run time: %g min", (ce - cs) / 60.);

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
