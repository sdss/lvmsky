/*
 *  This file is part of the SKYCORR software package.
 *  Copyright (C) 2009-2013 European Southern Observatory
 *
 *  This programme is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This programme is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this programme. If not, see <http://www.gnu.org/licenses/>.
 */


/*!
 * \ingroup sky_correction
 */

/**@{*/

/*!
 * \file sc_skycorr.c
 *
 * Routines for handling CMPFIT
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  18 Feb 2011
 * \date   08 Sep 2013
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <sc_skycorr.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

cpl_error_code sc_skycorr(const char *parfile_in)
{
    /*!
     * \callgraph
     *
     * This is the top-level routine of the sky correction code that subtracts
     * an optimised 1D sky spectrum from an 1D science spectrum. The
     * optimisation process is a CMPFIT-based fitting procedure that adapts
     * airglow lines in the sky spectrum to those in the science spectrum.
     * Object-related lines in the science spectrum should not be affected by
     * this procedure. The sky continuum is taken from the sky spectrum
     * without optimisation, since a separation of object and sky continuum in
     * the science spectrum is not possible. All required information for the
     * sky correction procedure is provided by an input ASCII parameter file.
     *
     * \b INPUT:
     * \param parfile_in  name of parameter file
     *
     * \b ERRORS:
     * - Error in subroutine
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameterlist *parlist;
    cpl_table *scispec, *skyspec, *groups, *linetab;
    char parfile[SC_MAXLEN], errtxt[SC_MAXLEN];

    /* Read driver file */
    sc_basic_absfile(parfile, parfile_in);
    parlist = cpl_parameterlist_new();
    if ((status = sc_par_readfile(parlist, parfile)) != CPL_ERROR_NONE) {
        cpl_parameterlist_delete(parlist);
        return status;
    }

    /* Make output directory if required */
    sc_skycorr_makeoutputdir(parlist);

    /* Read spectral and header data from science and sky data files */
    scispec = cpl_table_new(0);
    skyspec = cpl_table_new(0);
    if ((status = sc_skycorr_readspec(scispec, skyspec, parlist)) !=
        CPL_ERROR_NONE) {
        cpl_parameterlist_delete(parlist);
        cpl_table_delete(scispec);
        cpl_table_delete(skyspec);
        return status;
    }

    /* Read line groups from ASCII file and adapt it to observing
       conditions */
    groups = cpl_table_new(0);
    if ((status = sc_lines(groups, parlist)) != CPL_ERROR_NONE) {
        cpl_parameterlist_delete(parlist);
        cpl_table_delete(scispec);
        cpl_table_delete(skyspec);
        cpl_table_delete(groups);
        return status;
    }

    /* Subtraction of continuum from science and sky spectrum and estimation
       of line FWHM in sky spectrum */

    /* Prepare science spectrum */
    cpl_msg_info(cpl_func, "Science spectrum:");
    linetab = cpl_table_new(0);
    sc_specdiss_init_linetab(linetab);
    sc_skycorr_subcont(scispec, linetab, parlist, groups);
    cpl_table_delete(linetab);
    if ((status = cpl_error_get_code()) != CPL_ERROR_NONE) {
        sprintf(errtxt, "%s: error while separating lines and continuum of "
                "science spectrum", SC_ERROR_EIS_TXT);
        cpl_error_set_message(cpl_func, SC_ERROR_EIS, "%s", errtxt);
    }

    /* Prepare sky spectrum */
    cpl_msg_info(cpl_func, "Sky spectrum:");
    linetab = cpl_table_new(0);
    sc_specdiss_init_linetab(linetab);
    sc_skycorr_subcont(skyspec, linetab, parlist, groups);
    cpl_table_delete(linetab);
    if ((status = cpl_error_get_code()) != CPL_ERROR_NONE) {
        sprintf(errtxt, "%s: error while separating lines and continuum of "
                "sky spectrum", SC_ERROR_EIS_TXT);
        cpl_error_set_message(cpl_func, SC_ERROR_EIS, "%s", errtxt);
    }

    /* CMPFIT-based fitting procedure to adapt the sky line spectrum to the
       science line spectrum */
    sc_skycorr_fit(scispec, skyspec, groups, parlist);

    /* Free allocated memory */
    cpl_table_delete(skyspec);
    cpl_table_delete(groups);

    /* Stop programme in the case of errors */
    if ((status = cpl_error_get_code()) != CPL_ERROR_NONE) {
        cpl_parameterlist_delete(parlist);
        cpl_table_delete(scispec);
        return status;
    }

    /* Perform sky subtraction */
    sc_skycorr_subsky(scispec);

    /* Write and plot sky correction results */
    sc_skycorr_writespec(scispec, parlist);

    /* Free allocated memory */
    cpl_parameterlist_delete(parlist);
    cpl_table_delete(scispec);

    /* Return error code of last error */
    return cpl_error_get_code();
}


cpl_error_code sc_skycorr_makeoutputdir(const cpl_parameterlist *parlist)
{
    /*!
     * Makes output directory as given by the general CPL parameter list if
     * required.
     *
     * \b INPUT:
     * \param parlist  general CPL parameter list
     *
     * \b OUTPUT:
     * -
     *
     * \b ERRORS:
     * - none
     */

    const cpl_parameter *p;
    char basedir[SC_MAXLEN] = "", outdir[SC_MAXLEN] = "";

    /* Get output path */
    p = cpl_parameterlist_find_const(parlist, "inst_dir");
    strncpy(basedir, cpl_parameter_get_string(p), SC_MAXLEN);
    p = cpl_parameterlist_find_const(parlist, "output_dir");
    sc_basic_abspath(outdir, cpl_parameter_get_string(p), basedir);

    /* Create output directory if not present */
    if (access(outdir, W_OK) != 0) {
        cpl_msg_info(cpl_func, "Make output directory %s", outdir);
        if (mkdir(outdir, 0777)) {};
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_skycorr_readspec(cpl_table *scispec, cpl_table *skyspec,
                                   cpl_parameterlist *parlist)
{
    /*!
     * Reads science and sky spectrum from input files and checks similarity
     * of the wavelength grids of both spectra.
     *
     * \b INPUT:
     * \param parlist  general CPL parameter list
     *
     * \b OUTPUT:
     * \param scispec  CPL table with science spectrum
     * \param skyspec  CPL table with sky spectrum
     * \param parlist  general CPL parameter list with FITS header values
     *                 if the observing instrument is given
     *
     * \b ERRORS:
     * - Invalid object structure
     * - Inconsistent data grids
     * - see ::sc_readspec
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameter *p;
    char errtxt[SC_MAXLEN];
    int nsci = 0, nsky = 0;
    double lsci = 0., dsci = 0., psci = 0., csci = 0.;
    double lsky = 0., dsky = 0., psky = 0., csky = 0.;
    double dpix = 0., dlam0 = 0.;

    /* Read spectral and header data from science data file */
    p = cpl_parameterlist_find(parlist, "spectype");
    cpl_parameter_set_string(p, "SCI");
    if ((status = sc_readspec(scispec, parlist)) != CPL_ERROR_NONE) {
        return status;
    }

    /* Read spectral and header data from sky data file
       (header data from science data file is overwritten) */
    cpl_parameter_set_string(p, "SKY");
    if ((status = sc_readspec(skyspec, parlist)) != CPL_ERROR_NONE) {
        return status;
    }

    /* Check agreement in optional columns */
    if (cpl_table_has_column(skyspec, "dflux") !=
        cpl_table_has_column(scispec, "dflux")) {
        sprintf(errtxt, "%s: cpl_table *skyspec != cpl table *scispec "
                "(flux error column)", SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }
    if (cpl_table_has_column(skyspec, "mask") !=
        cpl_table_has_column(scispec, "mask")) {
        sprintf(errtxt, "%s: cpl_table *skyspec != cpl table *scispec "
                "(mask column)", SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Characterise wavelength grid of science spectrum */
    nsci = cpl_table_get_nrow(scispec);
    lsci = cpl_table_get(scispec, "lambda", 0, NULL);
    dsci = cpl_table_get(scispec, "lambda", nsci-1, NULL) - lsci;
    psci = dsci / (nsci - 1.);
    csci = lsci + dsci / 2.;

    /* Characterise wavelength grid of sky spectrum */
    nsky = cpl_table_get_nrow(skyspec);
    lsky = cpl_table_get(skyspec, "lambda", 0, NULL);
    dsky = cpl_table_get(skyspec, "lambda", nsky-1, NULL) - lsky;
    psky = dsky / (nsky - 1.);
    csky = lsky + dsky / 2.;

    /* Compare wavelength grids of both spectra */
    dpix = fabs(psky / psci - 1);
    dlam0 = fabs(csky / csci - 1);
    if (dpix > SC_GRIDTOL || dlam0 > SC_GRIDTOL) {
        sprintf(errtxt, "%s: cpl_table *scispec and *skyspec too different",
                SC_ERROR_IDG_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IDG, "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_skycorr_subcont(cpl_table *spec, cpl_table *linetab,
                                  cpl_parameterlist *parlist,
                                  const cpl_table *groups)
{
    /*!
     * \callgraph
     *
     * Subtracts continuum from a spectrum and derives FWHM of spectral lines
     * in this spectrum.
     *
     * The FWHM is derived by an iterative procedure that identifies isolated
     * lines, subtracts the continuum from these lines, and derives a clipped
     * mean FWHM. Then, the resulting line width is used for another line
     * search and so on. The procedure stops either if the relative change of
     * the FWHM between two iterations is below the 'ltol' value from the
     * input parameter list or after the 10th iteration.
     *
     * The list of line pixels identified by means of the first derivative of
     * the line flux is supplemented by pixels belonging to lines taken from
     * the airglow line lists. The number of pixels per line depends on the
     * FWHM derived before. Only lines are considered that have peak fluxes
     * above a threshold. If desired (parameter 'fluxlim' = -1), the flux
     * limit will be searched by an iterative procedure which tries to find
     * the best compromise between sufficient continuum pixels and a high
     * number of considered lines. In the case of a positive value for
     * 'fluxlim', this value is multiplied by the median flux of the already
     * identified lines in order to obtain the actual threshold. Factors in
     * the order of 0.01 are expected to work for most kind of spectra.
     *
     * \b INPUT:
     * \param spec     CPL table with spectrum
     * \param linetab  CPL table for information on detected lines
     * \param parlist  general CPL parameter list
     * \param groups   CPL table with airglow line information
     *
     * \b OUTPUT:
     * \param spec     spectrum with separation of lines and continuum
     * \param linetab  CPL table with information on detected lines
     * \param parlist  parameter list with updated line FWHM and relative line
     *                 flux limit
     *
     * \b ERRORS:
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameter *p;
    //cpl_parameterlist *plottags;
    cpl_array *mask;
    int it = 0, i = 0, n = 10, nrow = 0, j = 0, niso = 0, search = 0;
    double fwhm0 = 0., fwhm = 1000., eps = 1e-2, tmp = 0., rms = 0.;
    double fluxlim = 0., minfluxlim = 0.005, maxfluxlim = 0.1;

    /* Print info message */
    cpl_msg_info(cpl_func, "Identify lines, estimate FWHM, "
                 "and subtract continuum");

    /* Get FWHM from parameter list */
    p = cpl_parameterlist_find(parlist, "fwhm");
    fwhm0 = cpl_parameter_get_double(p);

    /* Take FWHM from parameter list as start value if already fitted */
    p = cpl_parameterlist_find(parlist, "iteration");
    it = cpl_parameter_get_int(p);
    if (it > 0) {
        fwhm = fwhm0;
    }

    /* Get convergence criterion for line width */
    p = cpl_parameterlist_find(parlist, "ltol");
    eps = cpl_parameter_get_double(p);

    /* Iterative line identification and FWHM estimation */

    for (tmp = fwhm, i = 0; i < n; i++, tmp = fwhm) {

        /* Find emission lines in spectrum */
        status = sc_specdiss_find_emissionlines(spec, linetab, parlist,
                                                groups);
        if (status != CPL_ERROR_NONE) {
            return status;
        }

        /* Subtract continuum in spectrum */
        status = sc_contsub(spec);
        if (status != CPL_ERROR_NONE) {
            return status;
        }

        /* Count isolated lines */
        nrow = cpl_table_get_nrow(linetab);
        for (j = 0, niso = 0; j < nrow; j++) {
            if (cpl_table_get_int(linetab, "isol_flag", j, NULL) != 0) {
                niso++;
            }
        }

        /* No isolated lines -> break */
        if (niso == 0) {
            if (it == 0) {
                fwhm = fwhm0;
            }
            cpl_msg_warning(cpl_func, "No isolated lines found -> "
                            "Take initial FWHM");
            break;
        }

        /* Estimate FWHM of lines in spectrum */
        status = sc_fwhmest(&fwhm, &rms, spec, linetab, parlist);
        if (status != CPL_ERROR_NONE) {
            return status;
        }

        /* Update FWHM in parlist */
        p = cpl_parameterlist_find(parlist, "fwhm");
        cpl_parameter_set_double(p, fwhm);

        /* Check convergence */
        if (fabs((fwhm - tmp) / tmp) < eps) {
            break;
        }

    }

    /* Save pixel class column */
    mask = cpl_array_new(0, CPL_TYPE_INT);
    sc_basic_col2arr(mask, spec, "class");

    /* Get relative peak flux limit from parameter list */
    p = cpl_parameterlist_find(parlist, "fluxlim");
    fluxlim = cpl_parameter_get_double(p);

    /* Check whether automatic search is desired */
    if (fluxlim < 0) {
        fluxlim = minfluxlim;
        search = 1;
    }

    /* Get line and continuum pixels (iterative approach if desired) */

    do {

        /* Get initial pixel class column */
        sc_basic_arr2col(spec, "class", mask);

        /* Identify continuum windows by means of airglow line list */
        status = sc_contsub_identcont(spec, groups, linetab, fluxlim,
                                      parlist);
        if (status != CPL_ERROR_NONE) {
            cpl_array_delete(mask);
            return status;
        }

        /* Modify lower flux limit for line peaks */
        fluxlim *= 2;

    } while (sc_contsub_check(spec) != 0 && fluxlim <= maxfluxlim &&
             search == 1);

    /* Info message on fluxlim parameter */
    if (search == 1) {
        cpl_msg_info(cpl_func, "Search for line threshold: FLUXLIM = %g",
                     fluxlim / 2);
    }

    /* Avoid identification of lines in the thermal IR */
    sc_contsub_skipthermir(spec);

    /* Cleanup */
    cpl_array_delete(mask);

    /* Plotting spectrum with identified line peaks */
    //plottags=cpl_parameterlist_new();
    //sc_setplottags_single_spec(plottags, "lambda", "cflux", "CFLUX before",
    //                           "lambda", "cflux", parlist);
    //sc_plot_single_spec(spec, plottags);
    //cpl_parameterlist_delete(plottags);

    /* Subtract continuum in spectrum */
    status = sc_contsub(spec);
    if (status != CPL_ERROR_NONE) {
        return status;
    }

    /* Plotting spectrum and continuum flux */
    //plottags = cpl_parameterlist_new();
    //sc_setplottags_overplot_spec(plottags, "lambda", "flux", "cflux", "Flux",
    //                             "Continuum", "", "Lambda", "Flux",
    //                             parlist);
    //sc_overplot_spec(spec, plottags);
    //cpl_parameterlist_delete(plottags);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_skycorr_fit(cpl_table *scispec, cpl_table *skyspec,
                              cpl_table *groups, cpl_parameterlist *parlist)
{
    /*!
     * \callgraph
     *
     * CMPFIT-based fitting procedure to adapt the sky line spectrum to the
     * science line spectrum.
     *
     * \b INPUT:
     * \param scispec  CPL table with science spectrum
     * \param skyspec  CPL table with sky spectrum
     * \param groups   CPL table with airglow line information
     * \param parlist  general CPL parameter list
     *
     * \b OUTPUT:
     * \param scispec  science and best-fit modified sky spectrum
     * \param skyspec  sky spectrum with line group weights
     * \param groups   modified line table
     * \param parlist  parameter list with additional entries
     *
     * \b ERRORS:
     * - see ::sc_mpfit
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_table *fitpar;
    mp_result result;
    int nsci = 0;
    double ts = 0., te = 0., runtime = 0.;

    /* Create columns for modification of sky spectrum */
    cpl_table_new_column(skyspec, "mlambda", CPL_TYPE_DOUBLE);
    cpl_table_new_column(skyspec, "mlflux", CPL_TYPE_DOUBLE);
    if (cpl_table_has_column(skyspec, "dflux") == 1) {
        cpl_table_new_column(skyspec, "mdflux", CPL_TYPE_DOUBLE);
    }
    cpl_table_new_column(skyspec, "mweight", CPL_TYPE_DOUBLE);
    cpl_table_new_column(skyspec, "mpix", CPL_TYPE_INT);
    cpl_table_new_column(skyspec, "dpix", CPL_TYPE_DOUBLE);
    cpl_table_new_column(skyspec, "frat", CPL_TYPE_DOUBLE);

    /* Create columns for comparison of science and sky spectrum */
    cpl_table_new_column(scispec, "mcflux", CPL_TYPE_DOUBLE);
    cpl_table_new_column(scispec, "mlflux", CPL_TYPE_DOUBLE);
    cpl_table_new_column(scispec, "mflux", CPL_TYPE_DOUBLE);
    if (cpl_table_has_column(skyspec, "dflux") == 1) {
        cpl_table_new_column(scispec, "mdflux", CPL_TYPE_DOUBLE);
    }
    if (cpl_table_has_column(skyspec, "mask") == 1) {
        cpl_table_new_column(scispec, "mmask", CPL_TYPE_INT);
    }
    cpl_table_new_column(scispec, "mweight", CPL_TYPE_DOUBLE);
    cpl_table_new_column(scispec, "sigclip", CPL_TYPE_INT);
    cpl_table_new_column(scispec, "cweight", CPL_TYPE_DOUBLE);
    cpl_table_new_column(scispec, "dev", CPL_TYPE_DOUBLE);

    /* Initialise sigclip column with 1 (= ejected) */
    nsci = cpl_table_get_nrow(scispec);
    cpl_table_fill_column_window(scispec, "sigclip", 0, nsci, 1);

    /* Initialise fit parameter table */
    fitpar = cpl_table_new(0);

    /* Prepare the line group weights of each pixel of the input sky
       spectrum */
    if ((status = sc_weights(skyspec, fitpar, groups, parlist)) !=
        CPL_ERROR_NONE) {
        cpl_table_delete(fitpar);
        return status;
    }

    /* Start run time measurement */
    ts = cpl_test_get_walltime();

    /* Start CMPFIT and compare resulting deviations */
    status = sc_mpfit(&result, scispec, skyspec, fitpar, parlist);

    /* Get CMPFIT run time in s */
    te = cpl_test_get_walltime();
    runtime = te - ts;

    /* Print fit results */
    cpl_msg_info(cpl_func, "FIT RESULTS:");
    cpl_msg_info(cpl_func, "status: %d", result.status);
    if (status == CPL_ERROR_NONE) {
        cpl_msg_info(cpl_func, "npar: %d", result.npar);
        cpl_msg_info(cpl_func, "npix: %d", result.nfunc);
        cpl_msg_info(cpl_func, "niter: %d", result.niter);
        cpl_msg_info(cpl_func, "nfev: %d", result.nfev);
        cpl_msg_info(cpl_func, "fittime: %.2f s", runtime);
        cpl_msg_info(cpl_func, "orignorm: %.3e", result.orignorm);
        cpl_msg_info(cpl_func, "bestnorm: %.3e", result.bestnorm);
    }

    /* Free allocated memory */
    cpl_table_delete(fitpar);
    sc_mpfit_freememresult(&result);

    return status;
}


cpl_error_code sc_skycorr_subsky(cpl_table *scispec)
{
    /*!
     * Subtracts sky continuum and fit-optimised sky line flux from science
     * spectrum. The resulting sky-corrected flux is written into the new
     * column "scflux". If an error occurred before, the flux correction is
     * not carried out. If flux error and mask column exist in the input
     * spectra, the effect of sky subtraction on these columns is also
     * considered. The resulting columns are labelled "scdflux" and "scmask"
     * in this case.
     *
     * \b INPUT:
     * \param scispec  CPL table with science spectrum and best-fit modified
     *                 sky spectrum
     *
     * \b OUTPUT:
     * \param scispec  data table with column for sky-subtracted flux
     *
     * \b ERRORS:
     * - none
     */

    cpl_error_code status = CPL_ERROR_NONE;
    int nrow = 0, i = 0;
    int *mmask = NULL, *mask = NULL, *scmask = NULL;
    double *mweight = NULL, *dflux = NULL, *mdflux = NULL, *scdflux = NULL;

    /* Get number of rows in science spectrum table */
    nrow = cpl_table_get_nrow(scispec);

    /* Calculate combined sky line and continuum flux */
    cpl_table_fill_column_window(scispec, "mflux", 0, nrow, 0.);
    cpl_table_add_columns(scispec, "mflux", "mcflux");
    cpl_table_add_columns(scispec, "mflux", "mlflux");

    /* Derive mask for sky spectrum if present */
    if (cpl_table_has_column(scispec, "mask") == 1) {
        cpl_table_fill_column_window(scispec, "mmask", 0, nrow, 0);
        mweight = cpl_table_get_data_double(scispec, "mweight");
        mmask = cpl_table_get_data_int(scispec, "mmask");
        for (i = 0; i < nrow; i++) {
            if (mweight[i] == 0.) {
                mmask[i] = 0;
            } else {
                mmask[i] = 1;
            }
        }
    }

    /* Create column for sky-subtracted flux and initialise it with
       uncorrected flux */
    cpl_table_duplicate_column(scispec, "scflux", scispec, "flux");

    /* Create column for flux errors after sky subtraction and initialise it
       with input flux errors if present */
    if (cpl_table_has_column(scispec, "dflux") == 1) {
         cpl_table_duplicate_column(scispec, "scdflux", scispec, "dflux");
    }

    /* Create column for mask after sky subtraction and initialise it with
       input mask if present */
    if (cpl_table_has_column(scispec, "mask") == 1) {
         cpl_table_duplicate_column(scispec, "scmask", scispec, "mask");
    }

    /* Leave routine in the case of errors */
    if ((status = cpl_error_get_code()) != CPL_ERROR_NONE) {
        return status;
    }

    /* Subtract sky continuum and best-fit sky line flux */
    cpl_table_subtract_columns(scispec, "scflux", "mflux");

    /* Derive errors for sky-subtracted flux if errors are present in the
       input data */
    if (cpl_table_has_column(scispec, "dflux") == 1) {
        dflux = cpl_table_get_data_double(scispec, "dflux");
        mdflux = cpl_table_get_data_double(scispec, "mdflux");
        scdflux = cpl_table_get_data_double(scispec, "scdflux");
        for (i = 0; i < nrow; i++) {
            scdflux[i] = sqrt(dflux[i] * dflux[i] + mdflux[i] * mdflux[i]);
        }
    }

    /* Derive combined mask if masks are present in the input data */
    if (cpl_table_has_column(scispec, "mask") == 1) {
        mask = cpl_table_get_data_int(scispec, "mask");
        mmask = cpl_table_get_data_int(scispec, "mmask");
        scmask = cpl_table_get_data_int(scispec, "scmask");
        for (i = 0; i < nrow; i++) {
            if (mask[i] == 1 && mmask[i] == 1) {
                scmask[i] = 1;
            } else {
                scmask[i] = 0;
            }
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_skycorr_writespec(cpl_table *scispec,
                                    cpl_parameterlist *parlist)
{
    /*!
     * \callgraph
     *
     * Writes results of sky correction procedure to two different files.
     * The 'fit' file contains all fitting results as FITS table. The '_SC'
     * file equals the input science file except for a sky-subtracted flux
     * column (or extension in the case of FITS images) and some new header
     * keywords. The routine also plots a comparison of the input science and
     * optimised sky spectra.
     *
     * \b INPUT:
     * \param scispec  CPL table with sky correction results
     * \param parlist  general CPL parameter list
     *
     * \b OUTPUT:
     * - none
     *
     * \b ERRORS:
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_table *results = NULL;
    cpl_propertylist *pheader = NULL, *header = NULL;
    cpl_parameter *p;
    char basedir[SC_MAXLEN] = "", outdir[SC_MAXLEN] = "";
    char outname[SC_MAXLEN] = "", outfits[SC_MAXLEN] = "";

    /* Set type of spectrum to science */
    p = cpl_parameterlist_find(parlist, "spectype");
    cpl_parameter_set_string(p, "SCI");

    /* Get output directory */
    p = cpl_parameterlist_find(parlist, "inst_dir");
    strncpy(basedir, cpl_parameter_get_string(p), SC_MAXLEN);
    p = cpl_parameterlist_find(parlist, "output_dir");
    sc_basic_abspath(outdir, cpl_parameter_get_string(p), basedir);

    /* Write info message */
    cpl_msg_info(cpl_func, "Write fit results into output folder %s", outdir);

    /* Write full CPL table scispec to FITS file */
    p = cpl_parameterlist_find(parlist, "output_name");
    strncpy(outname, cpl_parameter_get_string(p), SC_MAXLEN);
    sprintf(outfits, "%s%s_fit.fits", outdir, outname);
    pheader = cpl_propertylist_new();
    header = cpl_propertylist_new();
    cpl_propertylist_append_int(pheader, "DUMMY", 0);
    cpl_propertylist_append_int(header, "DUMMY", 0);
    cpl_table_save(scispec, pheader, header, outfits, CPL_IO_CREATE);
    cpl_propertylist_delete(pheader);
    cpl_propertylist_delete(header);

    /* Plot results */
    sc_skycorr_plot(scispec, parlist);

    /* Create science data file that equals the input file except for a
       sky-subtracted flux column and some new header keywords */

    /* Read file with SKYCORR results and convert it into CPL table */
    results = cpl_table_new(0);
    if ((status = sc_conv_readresults(results, parlist)) != CPL_ERROR_NONE) {
        cpl_table_delete(results);
        return status;
    }

    /* Write CPL table and CPL property list to file with format of input data
       file */
    status = sc_conv_writefile(results, parlist);

    /* Free allocated memory */
    cpl_table_delete(results);

    return status;
}


cpl_error_code sc_skycorr_plot(cpl_table *scispec,
                               const cpl_parameterlist *parlist)
{
    /*!
     * Compares the input science and the optimised sky spectrum. GNUPLOT is
     * used for plotting.
     *
     * \b INPUT:
     * \param scispec  CPL table with sky correction results
     * \param parlist  general CPL parameter list
     *
     * \b ERRORS:
     * - Invalid object structure
     */

    FILE *specfile1, *specfile2, *specfile3, *gnufile;
    const cpl_parameter *p;
    char errtxt[SC_MAXLEN], tmpdir[SC_MAXLEN], tmpfile1[SC_MAXLEN];
    char tmpfile2[SC_MAXLEN], tmpfile3[SC_MAXLEN];
    char plottype[SC_MAXLEN], basedir[SC_MAXLEN], outdir[SC_MAXLEN];
    char outname[SC_MAXLEN], psfile[SC_MAXLEN], tmpfile4[SC_MAXLEN];
    char systemcall[SC_MAXLEN];
    int nrow = 0, i = 0, existdir = 0, dummy = 0, j = 0, plotopt = 2;
    double xmin = 0., xmax = 0., ymin = 0., ymax = 0., dy = 0.;
    double lam = 0., flux1 = 0., weight1 = 0., flux2 = 0., weight2 = 0.;

    /* Labels of required table columns */
    char collam[] = "lambda";
    char colflux1[] = "flux";
    char colweight1[] = "weight";
    char colflux2[] = "mflux";
    char colweight2[] = "mweight";

    /* Extension of y-axis in per cent */
    double del = 0.05;

    /* Temporary filenames */
    char filename1[] = "obsspec.dat";
    char filename2[] = "modspec.dat";
    char filename3[] = "diffspec.dat";
    char gnuname[] = "plot.gnu";

    /* Plot labels */
    char xlabel[] = "Wavelength [micron]";
    char ylabel[] = "Radiance";
    char title[] = "Quality of sky subtraction";

    /* Check existence of required columns */
    if (cpl_table_has_column(scispec, collam) != 1 ||
        cpl_table_has_column(scispec, colflux1) != 1 ||
        cpl_table_has_column(scispec, colweight1) != 1 ||
        cpl_table_has_column(scispec, colflux2) != 1 ||
        cpl_table_has_column(scispec, colweight2) != 1) {
        sprintf(errtxt, "%s: cpl_table *scispec (required columns not found)",
                SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Exit if plot range is zero */
    nrow = cpl_table_get_nrow(scispec);
    if (nrow == 0) {
        return CPL_ERROR_NONE;
    }

    /* Get wavelength range for x-axis */
    xmin = cpl_table_get(scispec, collam, 0, NULL);
    xmax = cpl_table_get(scispec, collam, nrow-1, NULL);

    /* Get flux range for y-axis */
    for (ymin = 1000., ymax = 0., i = 0; i < nrow; i++) {
        flux1 = cpl_table_get(scispec, colflux1, i, NULL);
        weight1 = cpl_table_get(scispec, colweight1, i, NULL);
        if (flux1 < ymin && weight1 > 0.) {
            ymin = flux1;
        }
        if (flux1 > ymax && weight1 > 0.) {
            ymax = flux1;
        }
        flux2 = cpl_table_get(scispec, colflux2, i, NULL);
        weight2 = cpl_table_get(scispec, colweight2, i, NULL);
        if (flux2 < ymin && weight2 > 0.) {
            ymin = flux2;
        }
        if (flux2 > ymax && weight2 > 0.) {
            ymax = flux2;
        }
    }

    /* All weights = 0 */
    if (ymax < ymin) {
         return CPL_ERROR_NONE;
    }

    /* Extend y-axis range */
    dy = ymax - ymin;
    ymin -= del * dy;
    ymax += del * dy;

    /* Create temporary directory */
    sprintf(tmpdir, "__tmpDIRtmp__");
    existdir = access(tmpdir, W_OK);
    if (existdir == 0) {
        cpl_msg_warning(cpl_func, "Directory %s already exists!", tmpdir);
    } else {
        if ((dummy = mkdir(tmpdir, 0777))) {};
    }

    /* Write ASCII files containing observed, fitted, and difference
       spectra */
    sprintf(tmpfile1, "%s/%s", tmpdir, filename1);
    specfile1 = fopen(tmpfile1, "w");
    sprintf(tmpfile2, "%s/%s", tmpdir, filename2);
    specfile2 = fopen(tmpfile2, "w");
    sprintf(tmpfile3, "%s/%s", tmpdir, filename3);
    specfile3 = fopen(tmpfile3, "w");
    for (i = 0; i < nrow; i++) {
        lam = cpl_table_get(scispec, collam, i, NULL);
        flux1 = cpl_table_get(scispec, colflux1, i, NULL);
        weight1 = cpl_table_get(scispec, colweight1, i, NULL);
        flux2 = cpl_table_get(scispec, colflux2, i, NULL);
        weight2 = cpl_table_get(scispec, colweight2, i, NULL);
        fprintf(specfile1, "%5.6g\t%5.6g\n", lam, flux1);
        fprintf(specfile2, "%5.6g\t%5.6g\n", lam, flux2);
        if (weight1 != 0. && weight2 != 0.) {
            fprintf(specfile3, "%5.6g\t%5.6g\n", lam, flux1 - flux2);
        }
    }
    fclose(specfile1);
    fclose(specfile2);
    fclose(specfile3);

    /* Check plot options */
    p = cpl_parameterlist_find_const(parlist, "plot_type");
    sprintf(plottype, "%s", cpl_parameter_get_string(p));

    /* Get path and name of POSTSCRIPT file */
    p = cpl_parameterlist_find_const(parlist, "inst_dir");
    strncpy(basedir, cpl_parameter_get_string(p), SC_MAXLEN);
    p = cpl_parameterlist_find_const(parlist, "output_dir");
    sc_basic_abspath(outdir, cpl_parameter_get_string(p), basedir);
    p = cpl_parameterlist_find_const(parlist, "output_name");
    strncpy(outname, cpl_parameter_get_string(p), SC_MAXLEN);
    sprintf(psfile, "%s%s_fit.ps", outdir, outname);

    /* Create GNUPLOT driver files and run them */

    for (j = 0; j < 3; j++) {

        /* Plot type requested? */
        if (j == 0 && strchr(plottype, 'W') != NULL) {
            plotopt = 1;
        } else if (j == 1 && strchr(plottype, 'X') != NULL) {
            plotopt = 2;
        } else if (j == 2) {
            plotopt = 3;
        } else {
            continue;
        }

        /* Open temporary GNUPLOT file */
        sc_basic_initstring(tmpfile4, SC_MAXLEN);
        sprintf(tmpfile4, "%s/%s", tmpdir, gnuname);
        gnufile = fopen(tmpfile4, "w");

        /* Write lines dependent on plot type */
        if (plotopt == 1) {
            /* Plot on WXT terminal */
            fprintf(gnufile, "set term wxt\n");
            fprintf(gnufile, "set termoption enhanced\n");
        } else if (plotopt == 2) {
            /* Plot on X11 terminal */
            fprintf(gnufile, "set term x11\n");
            fprintf(gnufile, "set termoption enhanced\n");
        } else if (plotopt == 3) {
            /* POSTSCRIPT file */
            fprintf(gnufile, "set term postscript enhanced color\n");
            fprintf(gnufile, "set output \"%s\"\n", psfile);
            fprintf(gnufile, "set termoption font \"Times,12\"\n");
        }

        /* Write lines independent of plot type */

        fprintf(gnufile, "# Plotting\n");
        fprintf(gnufile, "set key at screen 0.70, 0.55 autotitle column box "
                         "samplen 1 left\n");
        fprintf(gnufile, "set tmargin 0\n");
        fprintf(gnufile, "set bmargin 5\n");
        fprintf(gnufile, "set lmargin 12\n");
        fprintf(gnufile, "set rmargin 14\n");
        fprintf(gnufile, "set xrange [%g:%g]\n", xmin, xmax);
        fprintf(gnufile, "set yrange [%g:%g]\n", ymin, ymax);
        fprintf(gnufile, "unset title\n");
        fprintf(gnufile, "set multiplot layout 2,1 title \"%s\"\n", title);
        fprintf(gnufile, "set xlabel \"%s\"\n", xlabel);
        fprintf(gnufile, "set ylabel \"%s\" offset 1,0\n", ylabel);
        fprintf(gnufile, "set style data boxes\n");
        fprintf(gnufile, "plot '%s' using 1:2 title \"input\" with "
                         "lines lt -1, '%s' using 1:2 title "
                         "\"best-fit sky\" with lines lt 8\n",
                         tmpfile1, tmpfile2);
        fprintf(gnufile, "set key at screen 0.66, 0.07 autotitle column box "
                         "samplen 1 left\n");
        fprintf(gnufile, "set style line 1 lt 2 lc rgb \"red\" lw 3\n");
        fprintf(gnufile, "set style line 2 lt 1 lc rgb \"green\" lw 1\n");
        fprintf(gnufile, "set ytics nomirror\n");
        fprintf(gnufile, "set ylabel \"Residual (input-best-fit sky)\" "
                         "offset 1,0\n");
        fprintf(gnufile, "set y2tics nomirror textcolor lt 2\n");
        fprintf(gnufile, "set y2label \"Residual (input-best-fit sky)\" "
                         "textcolor lt 2\n");
        fprintf(gnufile, "plot '%s' using 1:2 title "
                         "\"original scaling\" with lines lw 2 lt -1 "
                         "axes x1y1, '%s' using 1:2 title "
                         "\"optimal scaling\" with lines ls 2 axes x2y2\n",
                         tmpfile3, tmpfile3);
        fprintf(gnufile, "set autoscale y2\n");
        fprintf(gnufile, "set title\n");
        fprintf(gnufile, "unset multiplot\n");

        /* Close temoprary GNUPLOT file */
        fclose(gnufile);

        /* Call GNUPLOT */
        sprintf(systemcall, "gnuplot -persist %s", tmpfile4);
        dummy = system(systemcall);

        /* Remove temporary GNUPLOT file */
        dummy = remove(tmpfile4);

    }

    /* Remove ASCII files with spectra and delete temporary directory */
    dummy = remove(tmpfile1);
    dummy = remove(tmpfile2);
    dummy = remove(tmpfile3);
    dummy = rmdir(tmpdir);

    return CPL_ERROR_NONE;
}

/**@}*/
