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
 * \file testskymodel.c
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  14 Dec 2012
 * \date   06 Oct 2015
 *
 * Test programme for ETC sky model
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
     * Main routine for testing the ETC sky model. It has to be started from
     * the base directory of the sky model code by typing "bin/testskymodel".
     * It requires the parameter files "test_skymodel_etc.par" and
     * "test_instrument_etc.par" in the "test/" folder and returns a radiance
     * spectrum and transmission curve in FITS format in the folder "output/".
     * The output data is compared to "test_radspec.fits" and
     * "test_transspec.fits" in the "test/" folder.
     */

    char testdir[FILENAME_MAX] = "test/";
    char tmodelparfile[FILENAME_MAX] = "test_skymodel_etc.par";
    char tinstruparfile[FILENAME_MAX] = "test_instrument_etc.par";
    char configdir[FILENAME_MAX] = "config/";
    char modelparfile[FILENAME_MAX] = "skymodel_etc.par";
    char instruparfile[FILENAME_MAX] = "instrument_etc.par";
    char tradfitsfile[FILENAME_MAX] = "test_radspec.fits";
    char ttransfitsfile[FILENAME_MAX] = "test_transspec.fits";
    char outputdir[FILENAME_MAX] = "output/";
    char radfitsfile[FILENAME_MAX] = "radspec.fits";
    char transfitsfile[FILENAME_MAX] = "transspec.fits";
    char skymodelbin[FILENAME_MAX] = "bin/calcskymodel";

    cpl_errorstate errstate = cpl_errorstate_get();
    cpl_table *tradspec, *ttransspec, *radspec, *transspec;
    char tmodelpar[FILENAME_MAX], tinstrupar[FILENAME_MAX];
    char modelpar[FILENAME_MAX], instrupar[FILENAME_MAX];
    char tradfits[FILENAME_MAX], ttransfits[FILENAME_MAX];
    char radfits[FILENAME_MAX], transfits[FILENAME_MAX];
    char sys[FILENAME_MAX];
    int d = 0, mismatch = 0, count = 0, i = 0;
    double rflux = 0., flux = 0., rtrans = 0., trans = 0.;

    /* Initialise CPL */
    cpl_init(CPL_INIT_DEFAULT);

    /* Check existence of test/ folder */
    if (sm_basic_access(testdir, F_OK) != CPL_ERROR_NONE) {
        cpl_errorstate_dump(errstate, CPL_FALSE, cpl_errorstate_dump_one);
        cpl_end();
        return EXIT_FAILURE;
    }

    /* Build full names of test parameter files */
    sprintf(tmodelpar, "%s%s", testdir, tmodelparfile);
    sprintf(tinstrupar, "%s%s", testdir, tinstruparfile);

    /* Check existence of test model parameter file */
    if (sm_basic_access(tmodelpar, F_OK) != CPL_ERROR_NONE) {
        cpl_errorstate_dump(errstate, CPL_FALSE, cpl_errorstate_dump_one);
        cpl_end();
        return EXIT_FAILURE;
    }

    /* Check existence of test instrument parameter file */
    if (sm_basic_access(tinstrupar, F_OK) != CPL_ERROR_NONE) {
        cpl_errorstate_dump(errstate, CPL_FALSE, cpl_errorstate_dump_one);
        cpl_end();
        return EXIT_FAILURE;
    }

    /* Check existence of config/ folder */
    if (sm_basic_access(configdir, F_OK) != CPL_ERROR_NONE) {
        cpl_errorstate_dump(errstate, CPL_FALSE, cpl_errorstate_dump_one);
        cpl_end();
        return EXIT_FAILURE;
    }

    /* Build full names of input parameter files */
    sprintf(modelpar, "%s%s", configdir, modelparfile);
    sprintf(instrupar, "%s%s", configdir, instruparfile);

    /* Copy and rename test parameter files */

    cpl_msg_info(cpl_func, "Copy parameter files:");

    cpl_msg_info(cpl_func, "%s -> %s", tmodelpar, modelpar);
    sprintf(sys, "cp %s %s", tmodelpar, modelpar);
    d = system(sys);
    if (d != 0) {
        cpl_msg_error(cpl_func, "%s failed!", sys);
        cpl_end();
        return EXIT_FAILURE;
    }

    cpl_msg_info(cpl_func, "%s -> %s", tinstrupar, instrupar);
    sprintf(sys, "cp %s %s", tinstrupar, instrupar);
    d = system(sys);
    if (d != 0) {
        cpl_msg_error(cpl_func, "%s failed!", sys);
        cpl_end();
        return EXIT_FAILURE;
    }

    /* Build full names of reference output files */
    sprintf(tradfits, "%s%s", testdir, tradfitsfile);
    sprintf(ttransfits, "%s%s", testdir, ttransfitsfile);

    /* Build full names of output files */
    sprintf(radfits, "%s%s", outputdir, radfitsfile);
    sprintf(transfits, "%s%s", outputdir, transfitsfile);

    /* Check existence of executable for sky model calculation */
    if (sm_basic_access(skymodelbin, F_OK) != CPL_ERROR_NONE) {
        cpl_errorstate_dump(errstate, CPL_FALSE, cpl_errorstate_dump_one);
        cpl_end();
        return EXIT_FAILURE;
    }

    /* Call sky model programme */
    cpl_msg_info(cpl_func, "Run sky model code:");
    d = system(skymodelbin);
    if (d != 0) {
        cpl_msg_error(cpl_func, "Sky model calculation failed!");
        cpl_end();
        return EXIT_FAILURE;
    }

    /* Read FITS files and write data into CPL tables */

    cpl_msg_info(cpl_func, "Load reference and output spectra:");

    tradspec = cpl_table_load(tradfits, 1, 0);
    if (tradspec == NULL) {
        cpl_msg_error(cpl_func, "-> Loading of %s failed!", tradfits);
        cpl_end();
        return EXIT_FAILURE;
    }

    ttransspec = cpl_table_load(ttransfits, 1, 0);
    if (ttransspec == NULL) {
        cpl_table_delete(tradspec);
        cpl_msg_error(cpl_func, "-> Loading of %s failed!", ttransfits);
        cpl_end();
        return EXIT_FAILURE;
    }

    radspec = cpl_table_load(radfits, 1, 0);
    if (radspec == NULL) {
        cpl_table_delete(tradspec);
        cpl_table_delete(ttransspec);
        cpl_msg_error(cpl_func, "-> Loading of %s failed!", radfits);
        cpl_end();
        return EXIT_FAILURE;
    }

    transspec = cpl_table_load(transfits, 1, 0);
    if (transspec == NULL) {
        cpl_table_delete(tradspec);
        cpl_table_delete(ttransspec);
        cpl_table_delete(radspec);
        cpl_msg_error(cpl_func, "-> Loading of %s failed!", transfits);
        cpl_end();
        return EXIT_FAILURE;
    }

    cpl_msg_info(cpl_func, "-> OK");

    /* Compare CPL table structures */

    cpl_msg_info(cpl_func, "Compare table structures:");

    if (cpl_table_compare_structure(tradspec, radspec) != 0) {
        cpl_msg_error(cpl_func, "-> Structures of output and reference "
                      "radiance spectra deviate!");
        mismatch = 1;
    }

    if (cpl_table_compare_structure(ttransspec, transspec) != 0) {
        cpl_msg_error(cpl_func, "-> Structures of output and reference "
                      "transmission curves deviate!");
        mismatch = 1;
    }

    if (mismatch == 1) {
        cpl_table_delete(tradspec);
        cpl_table_delete(ttransspec);
        cpl_table_delete(radspec);
        cpl_table_delete(transspec);
        cpl_end();
        return EXIT_FAILURE;
    }  else {
        cpl_msg_info(cpl_func, "-> OK");
    }

    /* Compare fluxes in output and reference radiance spectra */

    cpl_msg_info(cpl_func, "Compare fluxes in output and reference radiance "
                 "spectra:");

    for (count = 0, i = 0; i < cpl_table_get_nrow(tradspec); i++) {
        rflux = cpl_table_get(tradspec, "flux", i, NULL);
        flux = cpl_table_get(radspec, "flux", i, NULL);
        if (abs(flux - rflux) > DBL_EPSILON) {
            mismatch = 1;
            count++;
        }
    }

    if (count != 0) {
        cpl_msg_error(cpl_func, "-> Flux deviates for %d pixels!", count);
    } else {
        cpl_msg_info(cpl_func, "-> OK");
    }

    /* Compare output and reference transmission */

    cpl_msg_info(cpl_func, "Compare output and reference transmission:");

    for (count = 0, i = 0; i < cpl_table_get_nrow(ttransspec); i++) {
        rtrans = cpl_table_get(ttransspec, "trans", i, NULL);
        trans = cpl_table_get(transspec, "trans", i, NULL);
        if (abs(trans - rtrans) > DBL_EPSILON) {
            mismatch = 1;
            count++;
        }
    }

    if (count != 0) {
        cpl_msg_error(cpl_func, "-> Transmission deviates for %d pixels!",
                      count);
    } else {
        cpl_msg_info(cpl_func, "-> OK");
    }

    /* Result of comparison */
    if (mismatch == 1) {
        cpl_msg_error(cpl_func, "Comparison failed!");
    } else {
        cpl_msg_info(cpl_func, "Comparison successful!");
    }

    /* Free allocated memory */
    cpl_table_delete(tradspec);
    cpl_table_delete(ttransspec);
    cpl_table_delete(radspec);
    cpl_table_delete(transspec);

    /* Finish CPL */
    cpl_end();
    if (mismatch == 1) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}

/**@}*/
