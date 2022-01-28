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
 * \file sc_mpfit.c
 *
 * Routines for handling CMPFIT
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  17 Feb 2011
 * \date   24 Feb 2014
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <sc_mpfit.h>


/*****************************************************************************
 *                                GLOBALS                                    *
 ****************************************************************************/

/* Definition of global variables */

/* Number of fitting function calls */
int nfev = 0;
/* Last modsky call? */
cpl_boolean lastcall = CPL_FALSE;


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

cpl_error_code sc_mpfit(mp_result *result, cpl_table *scispec,
                        cpl_table *skyspec, cpl_table *fitpar,
                        const cpl_parameterlist *parlist)
{
    /*!
     * \callgraph
     *
     * Handles the fitting routine CMPFIT (see mpfit.c for details on the
     * fitting algorithm).
     *
     * Needs the science line spectrum plus error-related weights, the sky
     * line spectrum plus error-related and line group weights, and
     * descriptions, fit flags, and initial values of the fit parameters
     * provided by a dedicated table and the general parameter list.
     *
     * The input data are used to compute a modified sky line spectrum which
     * is compared to the science line spectrum by deriving a vector of
     * weighted deviations. The sky spectrum is optimised by an initial
     * estimate of line group specific correction factors from the ratio of
     * science and sky line spectrum for selected, suitable line pixels.
     * As next step this estimate is repeated to improve the clipping of
     * unsuitable pixels. Finally, CMPFIT is called for a selection of line
     * group fit parameters which were flagged as crucial for a further fit
     * improvement by the previous flux correction estimate.
     *
     * The agreement of the wavelength grids of the science and sky spectrum
     * is improved by the fit of a Chebyshev polynomial by means of CMPFIT.
     * For obtaining an optimal polynomial the degree is increased
     * iteratively. For this reason, the fit plan described above is usually
     * repeated multiple times. After the first cycle the initial estimate
     * is substituted by the fit of the wavelength grid. In each run the
     * degree of the Chebyshev polynomial for the wavelength correction is
     * increased by 1 (exception: constant and linear term are considered
     * together) and the best-fit parameters of the previous run are used as
     * input. If \f${\chi^2}\f$ becomes worse or improves by less than the
     * relative value \e wtol, the loop is stopped and the best results of the
     * previous iterations are finally taken. Otherwise the procedure runs
     * until the maximum degree given by the parameter \e cheby_max is
     * reached. Moreover, the parameter \e cheby_min defines a minimum degree
     * for which a fit is performed. If \e cheby_min is greater than
     * \e cheby_max. A mode is applied which returns the results for the
     * degree \e cheby_max, even if the \f${\chi^2}\f$ is lower for an
     * iteration with a smaller degree.
     *
     * Information on the fit quality is written into a special results
     * structure. The CPL table with the science line spectrum is supplemented
     * by the best-fit sky line spectrum and the weighted deviations taken for
     * the \f${\chi^2}\f$ computation. Moreover, an ASCII file which
     * summarises the fit results (fit quality, run time, and best-fit
     * parameter values plus uncertainties) is written.
     *
     * \b INPUT:
     * \param scispec  CPL table with science spectrum
     * \param skyspec  CPL table with sky spectrum and line group weights
     * \param fitpar   CPL table of fit parameters
     * \param parlist  general CPL parameter list
     *
     * \b OUTPUT:
     * \param result   CMPFIT structure for fit results
     * \param scispec  CPL table with science and best-fit modified sky
     *                 spectrum
     * \param fitpar   fit parameter table with best-fit values
     *
     * \b ERRORS:
     * - Insufficient memory
     * - Error in subroutine
     */

    const cpl_parameter *pp;
    cpl_table *initfitpar = NULL, *tmpfitpar = NULL;
    cpl_vector *sinc = NULL;
    mp_config config = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    mp_result tmpresult;
    scpars fitpars;
    scvars v;
    char errtxt[SC_MAXLEN], type[SC_LENLINE+1];
    int m = 0, nmin = 0, nmax = 0, mode = 0, ncoef = 0, nloop = 0, pos = 0;
    int imin = 0, rebintype = 0, i = 0, status = 0, j = 0, calls = 0;
    int niter = 0, k = 0, nfit = 0, ibest = 0;
    double wtol = 0., norm = 0., dev = 0., orignorm = 0., ts = 0., te = 0.;
    double runtime = 0.;

    /* Set fit precision */
    pp = cpl_parameterlist_find_const(parlist, "ftol");
    config.ftol = cpl_parameter_get_double(pp);
    pp = cpl_parameterlist_find_const(parlist, "xtol");
    config.xtol = cpl_parameter_get_double(pp);

    /* Set maximum number of iterations */
    config.maxiter = 100;

    /* Get number of parameters */
    fitpars.n = cpl_table_get_nrow(fitpar);

    /* Allocate memory for scpars structure */
    if ((int) sc_mpfit_allocmempar(&fitpars, fitpars.n) ==
        (int) SC_ERROR_ISM) {
        result->status = -99;
        return SC_ERROR_ISM;
    }

    /* Get number of data points */
    m = cpl_table_get_nrow(scispec);

    /* Set return of fit residuals and parameter errors */
    if ((int) sc_mpfit_allocmemresult(result, m, fitpars.n) ==
        (int) SC_ERROR_ISM) {
        result->status = -99;
        sc_mpfit_freemempar(&fitpars);
        return SC_ERROR_ISM;
    }

    /* Initialise results structure */
    sc_mpfit_initresult(result, m, fitpars.n);

    /* Get minimum and maximum degree of polynomial and fitting mode for
       wavegrid correction */
    pp = cpl_parameterlist_find_const(parlist, "cheby_min");
    nmin = cpl_parameter_get_int(pp);
    pp = cpl_parameterlist_find_const(parlist, "cheby_max");
    nmax = cpl_parameter_get_int(pp);
    if (nmin > nmax) {
        mode = 1;
    }

    /* Allocate memory for temporary results structure if mode = 0 */
    if (mode == 0) {
        if ((int) sc_mpfit_allocmemresult(&tmpresult, m, fitpars.n) ==
            (int) SC_ERROR_ISM) {
            result->status = -99;
            sc_mpfit_freemempar(&fitpars);
            return SC_ERROR_ISM;
        }
    }

    /* Count coefficients for wavelength correction in fit parameter table and
       get maximum number of CMPFIT calls */
    if (nmax >= 0) {
        cpl_table_unselect_all(fitpar);
        cpl_table_or_selected_string(fitpar, "type", CPL_EQUAL_TO, "w");
        ncoef = cpl_table_count_selected(fitpar);
        nloop = ncoef;
        pos = fitpars.n - ncoef - 1;
        cpl_table_select_all(fitpar);
    } else {
        nloop = 1;
    }

    /* Get minimum iteration number */
    if (nmin < 1) {
        imin = 1;
    } else {
        imin = nmin;
    }

    /* Factor for chi^2 tolerance for estimate of polynomial degree */
    pp = cpl_parameterlist_find_const(parlist, "wtol");
    wtol = 1. + cpl_parameter_get_double(pp);

    /* Get method for rebinning from parameter list */
    pp = cpl_parameterlist_find_const(parlist, "rebintype");
    rebintype = cpl_parameter_get_int(pp);

    /* Pre-calculate damped sinc function for rebinning of data */
    sinc = cpl_vector_new(1);
    if (rebintype == 1) {
        sc_basic_calcsinc(sinc);
    }

    /* Info message */
    cpl_msg_info(cpl_func, "Fitting ...");

    /* Start time measurement */
    ts = cpl_test_get_walltime();

    /* Prepare science and sky spectrum table for first estimate */
    sc_modsky(scispec, skyspec, fitpar, sinc, parlist);

    /* Optimise initial fit parameters for line groups by considering the
       mean ratios of the science and sky line peak fluxes for each group */
    sc_mpfit_modinitpar(fitpar, scispec, skyspec, parlist);

    /* Update science and sky spectrum table */
    sc_modsky(scispec, skyspec, fitpar, sinc, parlist);

    for (i = 0; i < nloop; i++) {

        /* Loop info */
        if (i == 0) {
            cpl_msg_info(cpl_func, "Input wavelength grid (no fit)");
        } else {
            cpl_msg_info(cpl_func, "Chebyshev polynomial of degree %d", i);
        }

        /* Set errors to SC_DEFERRVAL */
        cpl_table_fill_column_window(fitpar, "err_fit", 0, fitpars.n,
                                     SC_DEFERRVAL);
        cpl_table_fill_column_window(fitpar, "err_est", 0, fitpars.n,
                                     SC_DEFERRVAL);

        /* Use CMPFIT to fit wavelength grid by a Chebyshev polynomial */

        if (ncoef == 0 || i == 0) {

            /* CMPFIT status for no fitting */
            status = 99;

            /* Calculate orignorm */
            for (norm = 0., j = 0; j < m; j++) {
                dev = cpl_table_get(scispec, "dev", j, NULL);
                norm += dev * dev;
            }
            result->orignorm = norm;

        } else {

            /* Set parameter vector and parameter constraints structure */
            sc_mpfit_setpar(&fitpars, fitpar, 'w');

            /* Pack all data and parameters into temporary structure in order
               to create a void pointer required by the CMPFIT user
               function */
            v.scispec = scispec;
            v.skyspec = skyspec;
            v.fitpar = fitpar;
            v.sinc = sinc;
            v.parlist = parlist;

            /* Call fitting function for m data points and n parameters */
            status = mpfit(sc_mpfit_calcdev, m, fitpars.n, fitpars.p,
                           fitpars.pars, &config, (void *) &v, result);

            if (status <= 0) {
                break;
            }

            /* Count number of CMPFIT calls and iterations */
            calls += 1;
            niter += result->niter;

            /* Write fit values to fit parameter table */
            cpl_table_copy_data_double(fitpar, "value", fitpars.p);

            /* Write fit errors to fit parameter table */
            for (k = 0; k < fitpars.n; k++) {
                sprintf(type, "%s", cpl_table_get_string(fitpar, "type", k));
                if (type[0] == 'w' &&
                    cpl_table_get_int(fitpar, "fit", k, NULL) == 1) {
                    cpl_table_set(fitpar, "err_fit", k, result->xerror[k]);
                }
            }

            /* Update science and sky spectrum table */
            sc_modsky(scispec, skyspec, fitpar, sinc, parlist);

        }

        /* Save orignorm */
        if (i == 0) {
            orignorm = result->orignorm;
        }

        /* Optimise fit parameters for line groups by considering the mean
           ratios of the science and sky line peak fluxes for each group */
        sc_mpfit_modinitpar(fitpar, scispec, skyspec, parlist);

        /* Update science and sky spectrum table */
        sc_modsky(scispec, skyspec, fitpar, sinc, parlist);

        /* Get number of free line group fit parameters */
        cpl_table_unselect_all(fitpar);
        cpl_table_or_selected_string(fitpar, "type", CPL_NOT_EQUAL_TO, "w");
        cpl_table_and_selected_int(fitpar, "fit", CPL_EQUAL_TO, 1);
        nfit = cpl_table_count_selected(fitpar);
        cpl_table_select_all(fitpar);

        /* Fit uncertain line group parameters if present */

        if (nfit == 0) {

            /* CMPFIT status for no fitting */
            status = 99;

            /* Calculate bestnorm if no CMPFIT call */
            for (norm = 0., j = 0; j < m; j++) {
                dev = cpl_table_get(scispec, "dev", j, NULL);
                norm += dev * dev;
            }
            result->bestnorm = norm;

        } else {

            /* Set parameter vector and parameter constraints structure */
            sc_mpfit_setpar(&fitpars, fitpar, 'l');

            /* Pack all data and parameters into temporary structure in order
               to create a void pointer required by the CMPFIT user
               function */
            v.scispec = scispec;
            v.skyspec = skyspec;
            v.fitpar = fitpar;
            v.sinc = sinc;
            v.parlist = parlist;

            /* Save initial fit parameter table */
            initfitpar = cpl_table_duplicate(fitpar);

            /* Call fitting function for m data points and n parameters */
            status = mpfit(sc_mpfit_calcdev, m, fitpars.n, fitpars.p,
                           fitpars.pars, &config, (void *) &v, result);

            /* Leave loop in the case of bad mpfit status */
            if (status <= 0) {
                cpl_table_delete(initfitpar);
                break;
            }

            /* Count number of CMPFIT calls and iterations */
            calls += 1;
            niter += result->niter;

            /* Write fit values to fit parameter table */
            cpl_table_copy_data_double(fitpar, "value", fitpars.p);

            /* Write fit errors to fit parameter table */
            for (k = 0; k < fitpars.n; k++) {
                sprintf(type, "%s", cpl_table_get_string(fitpar, "type", k));
                if (type[0] != 'w' &&
                    cpl_table_get_int(fitpar, "fit", k, NULL) == 1) {
                    cpl_table_set(fitpar, "err_fit", k, result->xerror[k]);
                }
            }

            /* Substitute uncertain fit parameters by initial estimates */
            sc_mpfit_substbadfitpar(fitpar, initfitpar, parlist);
            cpl_table_delete(initfitpar);

            /* Update science and sky spectrum table */
            sc_modsky(scispec, skyspec, fitpar, sinc, parlist);

        }

        /* Message on chi^2 */
        cpl_msg_info(cpl_func, "bestnorm: %.3e", result->bestnorm);

        /* Condition for break of loop */

        if (mode == 0 && ncoef > 0) {

            /* If at least minimum degree is reached and results become worse,
               stop loop */

            if (i >= imin && result->bestnorm * wtol >= tmpresult.bestnorm) {
                /* Take best results */
                if (result->bestnorm >= tmpresult.bestnorm) {
                    sc_mpfit_copyresult(result, &tmpresult);
                    sc_basic_copytable_content(fitpar, tmpfitpar);
                } else {
                    ibest = i;
                }
                break;
            } else if (i == nloop - 1 &&
                       result->bestnorm < tmpresult.bestnorm) {
                /* Take results of final degree */
                ibest = i;
            }

        } else {

            /* Case for run until maximum degree (no break) */

            if (i == nloop - 1) {
                ibest = i;
            }

        }

        /* Store fitting data in temporary structures if mode = 0 */
        if (mode == 0 && i < nloop - 1) {
            /* Store only better results than obtained by all previous
               iterations */
            if (i == 0 || result->bestnorm < tmpresult.bestnorm) {
                ibest = i;
                sc_mpfit_copyresult(&tmpresult, result);
                cpl_table_delete(tmpfitpar);
                tmpfitpar = cpl_table_duplicate(fitpar);
            }
        }

        /* Set relevance and fit flag = 1 for the next 'w' coefficient in the
           parameter table */
        if (i < nloop - 1) {
            if (i == 0) {
                /* Set constant term together with linear term */
                pos++;
                cpl_table_set(fitpar, "relevance", pos, 1);
                cpl_table_set(fitpar, "fit", pos, 1);
            }
            pos++;
            cpl_table_set(fitpar, "relevance", pos, 1);
            if (nmax != 0) {
                /* No fit of linear term if degree of poynomial is 0 */
                cpl_table_set(fitpar, "fit", pos, 1);
            }
            v.fitpar = fitpar;
        }

    }

    /* Print resulting degree of polynomial */
    if (ibest == 0) {
        cpl_msg_info(cpl_func, "STOP -> No wavegrid correction");
    } else {
        cpl_msg_info(cpl_func, "STOP -> Take results of degree %d", ibest);
    }

    /* Recover chi^2 before first fit and set final number of iterations */
    if (status > 0) {
        result->orignorm = orignorm;
        result->niter = niter;
    }

    /* Get CMPFIT run time in s */
    te = cpl_test_get_walltime();
    runtime = te - ts;

    /* Fill CPL table scispec with best-fit modified sky line spectrum and
       deviations */
    lastcall = CPL_TRUE;
    sc_modsky(scispec, skyspec, fitpar, sinc, parlist);

    /* Write a summary of the CMPFIT results into an ASCII file */
    if (status > 0) {
        sc_mpfit_writeresults(result, scispec, fitpar, parlist, calls,
                              runtime);
    }

    /* Free memory */
    if (status > 0 && mode == 0) {
        sc_mpfit_freememresult(&tmpresult);
        cpl_table_delete(tmpfitpar);
    }
    sc_mpfit_freemempar(&fitpars);
    cpl_vector_delete(sinc);

    /* Error message if CMPFIT fails */
    if (status <= 0) {
        result->status = status;
        sprintf(errtxt, "%s: mpfit (internal error %d)", SC_ERROR_EIS_TXT,
                status);
        return cpl_error_set_message(cpl_func, SC_ERROR_EIS, "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


int sc_mpfit_calcdev(int m, int n, double *p, double *dy, double **dvec,
                     void *vars)
{
    /*!
     * \callgraph
     *
     * User function for CMPFIT. Provides spectral data and fit parameters
     * to ::sc_modsky. Gets weighted deviations between the modified
     * sky line spectrum and the science line spectrum from ::sc_modsky.
     * The syntax of the function is predefined by CMPFIT.
     *
     * \b INPUT:
     * \param m     number of data points
     * \param n     number of parameters
     * \param p     array of fit parameters
     * \param dvec  derivatives (not used)
     * \param vars  private data -> observed spectrum and driver file
     *                              parameters
     *
     * \b OUTPUT:
     * \param dy    array of residuals ([model - obs. spectrum] * weight)
     *
     * \b ERRORS:
     * - none
     */

    scvars *v = (scvars *) vars;
    cpl_table *scispec, *skyspec, *fitpar;
    cpl_vector *sinc;
    const cpl_parameterlist *parlist;

    if (n) {};
    if (dvec) {};

    /* Update number of fitting function call */
    nfev++;

    /* Unpack observed spectral data and input parameters */
    scispec = v->scispec;
    skyspec = v->skyspec;
    fitpar = v->fitpar;
    sinc = v->sinc;
    parlist = v->parlist;

    /* Put fit parameters in CPL table */
    cpl_table_copy_data_double(fitpar, "value", p);

    /* Modification of sky spectrum */
    sc_modsky(scispec, skyspec, fitpar, sinc, parlist);

    /* Fill array of residuals */
    assert(m == cpl_table_get_nrow(scispec));
    cpl_table_fill_invalid_double(scispec, "dev", 0.);
    memcpy(dy, cpl_table_get_data_double_const(scispec, "dev"),
           cpl_table_get_nrow(scispec) * sizeof(double));

    return CPL_ERROR_NONE;
}


cpl_error_code sc_mpfit_modinitpar(cpl_table *fitpar,
                                   cpl_table *scispec,
                                   cpl_table *skyspec,
                                   const cpl_parameterlist *parlist)
{
    /*!
     * Optimises initial fit parameters for line groups by considering the
     * mean ratio of group-specific line fluxes in the science and sky
     * spectrum. The pixels identified as line peaks in the science spectrum
     * (see ::sc_specdiss_find_emissionlines) and those at a distance of not
     * more than FWHM / 2 (computed by ::sc_fwhmest) from these peaks are
     * taken for the ratio calculation only. Moreover, outliers (e.g. object
     * lines) indicating a deviation of the ratio of more than \e siglim
     * \f${\sigma}\f$ (see input parameter list) from the mean ratio of the
     * line peaks (as computed by ::sc_basic_clipmean) are excluded. Then,
     * only those lines are used that dominate a spectrum pixel by a weight of
     * at least the value of the input parameter \e weightlim. This criterion
     * is separately applied to A and B line groups. A \f${\sigma}\f$-clipping
     * procedure using ::sc_basic_clipmean is also performed for the pixels of
     * each group. If the resulting clipped pixels of a line outnumber the
     * valid ones by a factor of 2 and more, the entire line is clipped. If it
     * turns out that a line group does not have suitable lines for estimating
     * the line correction factor, the mean ratio of all lines is taken for A
     * groups and a value of 1 is taken for B groups.
     *
     * The routine prepares the use of CMPFIT. Apart from the start values for
     * the line flux correction factors, the selected line peak pixels are
     * also set for the fitting procedure. No other pixels are considered for
     * the \f${\chi^2}\f$ calculation. Moreover, it is decided which
     * parameters are flagged as free and which parameters are hold. Line
     * groups for which the initial estimate failed are also not considered
     * for the fitting procedure. On the other hand, the user can influence
     * the selection by the parameter \e fitlim which only selects line groups
     * that indicate a ratio of RMS to mean correction factor of at least the
     * value given. Consequently, the fitting can be carried out for the
     * uncertain but fittable groups only.
     *
     * \b INPUT:
     * \param fitpar   CPL table of fit parameters
     * \param scispec  CPL table with science spectrum
     * \param skyspec  CPL table with sky spectrum
     * \param parlist  input CPL parameter list
     *
     * \b OUTPUT:
     * \param fitpar   fit parameter table with adapted start values for line
     *                 groups
     *
     * \b ERRORS:
     * - none
     */

    /* Parameter list */
    const cpl_parameter *p;
    double fwhm = 0., siglim = 0., fitlim = 0.;

    /* Science spectrum */
    cpl_array *ratio, *tmpratio;
    cpl_boolean ispeak = CPL_FALSE;
    int nsci = 0, pix = 0, npeak = 0;
    double lflux = 0., mlflux = 0., rat = 0., initrat = 0., rms = 0.;

    /* Sky spectrum */
    char dcolname[SC_LENLINE+1];
    int nsky = 0, i = 0, io = 0, im = 0;

    /* Fit parameters */
    cpl_table *tmpfitpar;
    int row = 0;

    /* Line group type */
    char grouptype[2] = "AB", ngroupid[SC_LENLINE+1];
    int h = 0;

    /* Line group systems */
    cpl_array *corrsyst, *corrfac, *meanratio, *groupsum, *maxrms[2];
    int nsyst = 0, syst = 0, gsum = 0;
    double ratsum = 0., fac = 0., meanrat = 0.;

    /* Line groups */
    cpl_array *gratio, *sgratio;
    int ngroup[2] = {0, 0}, ming = 0, g = 0;
    double mrat = 0., msig = 0., omrat = 0.;
    double nmrat = 0.;

    /* Lines */
    cpl_array *selpeak;
    int maxpix = 0, d = 0, nlin = 0, n = 0, k = 0;

    /* Pixels */
    cpl_array *selpix, *linid;
    int sum = 0, j = 0, nlinpix = 0, nclippix = 0;

    /* Get number of rows in science and sky spectrum */
    nsci = cpl_table_get_nrow(scispec);
    nsky = cpl_table_get_nrow(skyspec);

    /* Ratio of science and sky line flux */

    ratio = cpl_array_new(nsci, CPL_TYPE_DOUBLE);

    for (pix = 0, npeak = 0; pix < nsci; pix++) {
        if (cpl_table_get(scispec, "class", pix, NULL) >= 1 &&
            cpl_table_get(scispec, "weight", pix, NULL) > 0. &&
            cpl_table_get(scispec, "mweight", pix, NULL) > 0.) {
            /* Line pixel with non-zero weight in science spectrum */
            lflux = cpl_table_get(scispec, "lflux", pix, NULL);
            mlflux = cpl_table_get(scispec, "mlflux", pix, NULL);
            if (mlflux == 0) {
                rat = 1e6;
            } else {
                rat = lflux / mlflux;
            }
            cpl_array_set_double(ratio, pix, rat);
            /* Count line peaks */
            if (cpl_table_get(scispec, "class", pix, NULL) >= 2) {
                npeak += 1;
            }
        } else {
            /* Set array element invalid if pixel is not suitable */
            cpl_array_set_invalid(ratio, pix);
        }
    }

    /* Do not modify parameter values in the case of no line peaks */
    if (npeak == 0) {
        cpl_array_delete(ratio);
        return CPL_ERROR_NONE;
    }

    /* FWHM from parameter list */
    p = cpl_parameterlist_find_const(parlist, "fwhm");
    fwhm = cpl_parameter_get_double(p);

    /* Number of pixels within FWHM / 2 */
    maxpix = (int) floor(fwhm / 2. + 0.5);

    /* Use only line pixels within a distance of FWHM / 2 from peak */
    for (pix = 0; pix < nsci; pix++) {
        if (cpl_array_is_valid(ratio, pix) == 1 &&
            cpl_table_get(scispec, "class", pix, NULL) == 1) {
            ispeak = CPL_FALSE;
            for (d = -maxpix; d <= maxpix; d++) {
                if (pix + d >= 0 && pix + d < nsci &&
                    cpl_table_get(scispec, "class", pix + d, NULL) >= 2 &&
                    cpl_array_is_valid(ratio, pix + d) == 1) {
                    ispeak = CPL_TRUE;
                    break;
                }
            }
            if (ispeak == CPL_FALSE) {
                cpl_array_set_invalid(ratio, pix);
            }
        }
    }

    /* Save input fit parameter table */
    tmpfitpar = cpl_table_duplicate(fitpar);

    /* Temporary ratio array for sigma-clipping */
    tmpratio = cpl_array_duplicate(ratio);

    /* Use line peaks for mean ratio calculation only */
    for (pix = 0; pix < nsci; pix++) {
        if (cpl_array_is_valid(tmpratio, pix) == 1 &&
            cpl_table_get(scispec, "class", pix, NULL) == 1) {
            cpl_array_set_invalid(tmpratio, pix);
        }
    }

    /* Perform sigma-clipping for entire spectrum */
    sc_basic_clipmean(&initrat, &rms, tmpratio, CPL_FALSE);

    /* Get parameter for sigma limit */
    p = cpl_parameterlist_find_const(parlist, "siglim");
    siglim = cpl_parameter_get_double(p);

    /* Set elements of ratio array invalid that show deviations greater than
       siglim sigma from mean */
    for (pix = 0; pix < nsci; pix++) {
        rat = cpl_array_get(ratio, pix, NULL);
        if (rat < initrat - siglim * rms ||
            rat > initrat + siglim * rms) {
            cpl_array_set_invalid(ratio, pix);
        }
    }

    /* Initialise flux ratio column with invalid elements */
    cpl_table_set_column_invalid(skyspec, "frat", 0, nsky);

    /* Get ratio of line flux in science and sky spectrum for valid pixels */
    for (i = 0; i < nsky; i++) {
        pix = cpl_table_get(skyspec, "mpix", i, NULL);
        if (pix >= 0 && pix < nsci && cpl_array_is_valid(ratio, pix) == 1) {
            rat = cpl_array_get(ratio, pix, NULL);
            cpl_table_set(skyspec, "frat", i, rat);
        }
    }

    /* Create array for system B group correction flag  */
    p = cpl_parameterlist_find_const(parlist, "n_system");
    nsyst = cpl_parameter_get_int(p);
    corrsyst = cpl_array_new(nsyst, CPL_TYPE_DOUBLE);
    cpl_array_fill_window(corrsyst, 0, nsyst, 1);

    /* Create array for correction factors of "0" groups */
    corrfac = cpl_array_new(nsyst, CPL_TYPE_DOUBLE);

    /* Create system-related array for mean factor and number of reliable
       A groups */
    meanratio = cpl_array_new(nsyst + 1, CPL_TYPE_DOUBLE);
    cpl_array_fill_window(meanratio, 0, nsyst + 1, 0.);
    groupsum = cpl_array_new(nsyst + 1, CPL_TYPE_INT);
    cpl_array_fill_window(groupsum, 0, nsyst + 1, 0);

    /* Get parameter for relative RMS limit for setup of fit parameters */
    p = cpl_parameterlist_find_const(parlist, "fitlim");
    fitlim = cpl_parameter_get_double(p);

    /* Get numbers of dominating A and B groups for each pixel in sky
       spectrum */

    for (h = 0; h < 2; h++) {

        /* Get number of A or B groups */
        sprintf(ngroupid, "n_group%c", grouptype[h]);
        p = cpl_parameterlist_find_const(parlist, ngroupid);
        ngroup[h] = cpl_parameter_get_int(p);

        /* No data -> skip group type */
        if (ngroup[h] == 0) {
            continue;
        }

        /* Build name for column of dominating group IDs */
        sprintf(dcolname, "dg%c", grouptype[h]);

        /* Get minimum group number */
        if (h == 1) {
            /* n_system negative B group numbers */
            p = cpl_parameterlist_find_const(parlist, "n_system");
            ming = -cpl_parameter_get_int(p);
        } else {
            ming = 0;
        }

        /* Create array for maximum RMS of A or B group systems */
        maxrms[h] = cpl_array_new(nsyst, CPL_TYPE_DOUBLE);
        cpl_array_fill_window(maxrms[h], 0, nsyst, 0.);

        /* Get mean flux ratio for each group */

        for (g = ming; g <= ngroup[h]; g++) {
            int * pig;
            /* Get pixel numbers of group pixels and line peaks */
            cpl_table_unselect_all(skyspec);
            cpl_table_or_selected_int(skyspec, dcolname, CPL_EQUAL_TO, g);
            cpl_table_and_selected_double(skyspec, "frat", CPL_NOT_EQUAL_TO,
                                          HUGE_VAL);
            selpix = cpl_table_where_selected(skyspec);
            cpl_table_and_selected_int(skyspec, "class",
                                       CPL_NOT_LESS_THAN, 2);
            selpeak = cpl_table_where_selected(skyspec);
            cpl_table_select_all(skyspec);

            /* Get number of group pixels and lines */
            sum = cpl_array_get_size(selpix);
            nlin = cpl_array_get_size(selpeak);

            /* Get mean flux ratio */

            if (nlin == 0) {

                /* No correction for groups without suitable lines */
                mrat = 1.;
                msig = SC_DEFERRVAL;
                n = 0;

            } else {

                /* Assign selected pixels to lines (IDs: increasing numbers
                   from 1 to N_lin) */
                linid = cpl_array_new(sum, CPL_TYPE_INT);
                io = cpl_array_get(selpeak, 0, NULL);
                for (j = 0, k = 1; k <= nlin; k++) {
                    if (k == nlin) {
                        im = nsky;
                    } else {
                        i = cpl_array_get(selpeak, k, NULL);
                        im = (int) ceil(0.5 * (io + i));
                    }
                    while (j < sum) {
                        if (cpl_array_get(selpix, j, NULL) < im) {
                            cpl_array_set(linid, j, k);
                        } else {
                            break;
                        }
                        j++;
                      }
                    io = i;
                }

                /* Create flux ratio array for line group */
                gratio = cpl_array_new(sum, CPL_TYPE_DOUBLE);
                for (j = 0; j < sum; j++) {
                    i = cpl_array_get(selpix, j, NULL);
                    cpl_array_set(gratio, j,
                                  cpl_table_get(skyspec, "frat", i, NULL));
                }

                /* Save flux ratio array */
                sgratio = cpl_array_duplicate(gratio);

                /* Sigma-clipping procedure for N_pixel >= 3 */
                if (sum >= 3) {
                    sc_basic_clipmean(&mrat, &msig, gratio, CPL_FALSE);
                }

                /* Correct clipping of pixels for each line */

                for (n = 0, k = 0; k < nlin; k++) {

                    /* Count valid and clipped line pixels */
                    for (nlinpix = 0, nclippix = 0, j = 0; j < sum; j++) {
                        if (cpl_array_get(linid, j, NULL) == k + 1) {
                            if (cpl_array_is_valid(gratio, j) == 1) {
                                nlinpix++;
                            } else {
                                nclippix++;
                            }
                        }
                    }

                    /* Avoid complete clipping of one line if only two lines
                       are available */
                    if (nlin == 2 && nlinpix == 0) {
                        for (j = 0; j < sum; j++) {
                            if (cpl_array_get(linid, j, NULL) == k + 1 &&
                                cpl_array_is_valid(gratio, j) != 1) {
                                cpl_array_set(gratio, j,
                                              cpl_array_get(sgratio, j,
                                                            NULL));
                                nlinpix++;
                                nclippix--;
                            }
                        }
                    }

                    /* Clip complete line if clipped pixels of a line
                       outnumber valid pixels by a factor of 2 and more */
                    if (nclippix > 2 * nlinpix) {
                        for (nlinpix = 0, j = 0; j < sum; j++) {
                            if (cpl_array_get(linid, j, NULL) == k + 1 &&
                                cpl_array_is_valid(gratio, j) == 1) {
                                cpl_array_set_invalid(gratio, j);
                                nlinpix--;
                                nclippix++;
                            }
                        }
                    }

                    /* Count number of unclipped lines of a group */
                    if (nlinpix > 0) {
                        n++;
                    }

                }

                /* Get group-specific mean and RMS */
                if (sum == 1) {
                    mrat = cpl_array_get(gratio, 0, NULL);
                    msig = 0;
                } else {
                    mrat = cpl_array_get_mean(gratio);
                    msig = cpl_array_get_stdev(gratio);
                }

                /* Consider unclipped line pixels for RMS calculation */
                for (j = 0; j < sum; j++) {
                    if (cpl_array_is_valid(gratio, j) == 1) {
                        i = cpl_array_get(selpix, j, NULL);
                        pix = cpl_table_get(skyspec, "mpix", i, NULL);
                        if (h == 0 && g > 0 && pix >= 0 && pix < nsci) {
                            cpl_table_set(scispec, "sigclip", pix, 0);
                        }
                    }
                }

                /* Delete temporary array */
                cpl_array_delete(linid);
                cpl_array_delete(sgratio);
                cpl_array_delete(gratio);

                /* Make sure that ratio is not beyond limits */
                if (mrat < SC_CORRFAC_MIN || mrat > SC_CORRFAC_MAX) {
                    mrat = 1.;
                }

            }

            /* Delete temporary array */
            cpl_array_delete(selpix);
            cpl_array_delete(selpeak);

            /* Divide flux ratios by derived group-specific mean ratios */
            pig = cpl_table_get_data_int(skyspec, dcolname);
            for (i = 0; i < nsky; i++) {
                if (pig[i] == g && cpl_table_is_valid(skyspec, dcolname, i) &&
                    cpl_table_is_valid(skyspec, "frat", i)) {
                    rat = cpl_table_get(skyspec, "frat", i, NULL);
                    cpl_table_set(skyspec, "frat", i, rat / mrat);
                }
            }

            /* Row in fit parameter list */
            if (h == 0) {
                row = g - 1;
            } else if (h == 1) {
                row = ngroup[0] + g - 1;
            }

            /* Set system correction flag to 0 (no correction) if a B group
               of a system cannot be fitted */
            if (h == 1 && g > 0 && nlin == 0) {
                syst = cpl_table_get(fitpar, "system", row, NULL);
                if (syst > 0) {
                    cpl_array_set(corrsyst, syst-1, 0);
                }
            }

            /* Get maximum RMS for all A or B groups of a system */
            if (g > 0 && n >= 2) {
                syst = cpl_table_get(fitpar, "system", row, NULL);
                if (syst > 0) {
                    if (rms > cpl_array_get(maxrms[h], syst-1, NULL)) {
                        cpl_array_set(maxrms[h], syst-1, rms);
                    }
                }
            }

            if (g > 0) {
                /* Set initial correction factor for group fluxes */
                cpl_table_set(fitpar, "value", row, mrat);
                /* Put RMS in estimate error column */
                cpl_table_set(fitpar, "err_est", row, msig);
                /* Put number of valid lines in N_lin column */
                cpl_table_set(fitpar, "N_lin", row, n);
            } else if (g < 0) {
                /* Save correction factors for "0" groups */
                cpl_array_set(corrfac, -g-1, mrat);
            }

        }

    }

    /* Set RMS of groups with only one valid line to maximum of the system */
    for (row = 0; row < ngroup[0] + ngroup[1]; row++) {
        syst = cpl_table_get(fitpar, "system", row, NULL);
        n = cpl_table_get(fitpar, "N_lin", row, NULL);
        if (row >= ngroup[0]) {
            h = 1;
        } else {
            h = 0;
        }
        if (syst > 0 && n == 1) {
            cpl_table_set(fitpar, "err_est", row,
                          cpl_array_get(maxrms[h], syst-1, NULL));
        }
    }

    /* No correction of B groups of systems with correction flag = 0 */
    for (row = ngroup[0]; row < ngroup[0] + ngroup[1]; row++) {
        syst = cpl_table_get(fitpar, "system", row, NULL);
        if (syst > 0 && cpl_array_get(corrsyst, syst-1, NULL) == 0) {
            cpl_table_set(fitpar, "value", row, 1.);
            cpl_table_set(fitpar, "err_est", row, SC_DEFERRVAL);
        }
    }

    /* Modify flux correction factors in order to get a B group factor of 1
       for the "0" group lines that do not belong to a B group */
    for (row = 0; row < ngroup[0] + ngroup[1]; row++) {
        syst = cpl_table_get(fitpar, "system", row, NULL);
        if (syst == 0) {
            continue;
        }
        mrat = cpl_table_get(fitpar, "value", row, NULL);
        msig = cpl_table_get(fitpar, "err_est", row, NULL);
        fac = cpl_array_get(corrfac, syst-1, NULL);
        if (fac != 0 && msig != SC_DEFERRVAL &&
            cpl_array_get(corrsyst, syst-1, NULL) != 0) {
            if (row < ngroup[0]) {
                cpl_table_set(fitpar, "value", row, mrat * fac);
                cpl_table_set(fitpar, "err_est", row, msig * fac);
            } else {
                cpl_table_set(fitpar, "value", row, mrat / fac);
                cpl_table_set(fitpar, "err_est", row, msig / fac);
            }
        }
    }

    /* Convert relative correction factors into absolute ones */
    for (row = 0; row < ngroup[0] + ngroup[1]; row++) {
        n = cpl_table_get(fitpar, "N_lin", row, NULL);
        if (n > 0) {
            omrat = cpl_table_get(tmpfitpar, "value", row, NULL);
            mrat = cpl_table_get(fitpar, "value", row, NULL);
            nmrat = omrat * mrat;
            if (nmrat < SC_CORRFAC_MIN) {
        nmrat = SC_CORRFAC_MIN;
            } else if (nmrat > SC_CORRFAC_MAX) {
                nmrat = SC_CORRFAC_MAX;
            }
            cpl_table_set(fitpar, "value", row, nmrat);
            msig = cpl_table_get(fitpar, "err_est", row, NULL);
            if (msig != SC_DEFERRVAL) {
                cpl_table_set(fitpar, "err_est", row, omrat * msig);
            }
        }
    }

    /* Calculate mean value of A groups of a system (only consider groups with
       at least SC_MINNLIN lines) */
    for (row = 0; row < ngroup[0]; row++) {
        n = cpl_table_get(fitpar, "N_lin", row, NULL);
        if (n >= SC_MINNLIN) {
            syst = cpl_table_get(fitpar, "system", row, NULL);
            mrat = cpl_table_get(fitpar, "value", row, NULL);
            if (syst > 0) {
                ratsum = cpl_array_get(meanratio, syst, NULL) + mrat;
                cpl_array_set(meanratio, syst, ratsum);
                gsum = cpl_array_get(groupsum, syst, NULL) + 1;
                cpl_array_set(groupsum, syst, gsum);
            }
            /* Take mean of all suitable group factors for system = 0 */
            ratsum = cpl_array_get(meanratio, 0, NULL) + mrat;
            cpl_array_set(meanratio, 0, ratsum);
            gsum = cpl_array_get(groupsum, 0, NULL) + 1;
            cpl_array_set(groupsum, 0, gsum);
        }
    }

    /* Get mean ratio for A groups of a system */
    cpl_array_divide(meanratio, groupsum);

    /* For A groups without a valid line take mean factor of the system  */
    for (row = 0; row < ngroup[0]; row++) {
        n = cpl_table_get(fitpar, "N_lin", row, NULL);
        if (n == 0) {
            syst = cpl_table_get(fitpar, "system", row, NULL);
            if (cpl_array_is_valid(meanratio, syst) != 1) {
                /* Take global ratio if no system-specific mean */
                syst = 0;
            }
            if (cpl_array_is_valid(meanratio, syst) != 1) {
                /* Take initial ratio if no global mean */
                cpl_table_set(fitpar, "value", row, 1.);
            } else {
                cpl_table_set(fitpar, "value", row,
                              cpl_array_get(meanratio, syst, NULL));
            }
        }
    }

    /* No correction of a system A group with less than SC_MINNLIN lines if
       deviation of scaling factor from system mean is greater than allowed by
       SC_MAXRELFAC */
    for (row = 0; row < ngroup[0]; row++) {
        n = cpl_table_get(fitpar, "N_lin", row, NULL);
        syst = cpl_table_get(fitpar, "system", row, NULL);
        if (n > 0 && n < SC_MINNLIN && syst > 0 &&
            cpl_array_is_valid(meanratio, syst) == 1) {
            mrat = cpl_table_get(fitpar, "value", row, NULL);
            meanrat = cpl_array_get(meanratio, syst, NULL);
            if (meanrat == 0 || mrat / meanrat < 1. / SC_MAXRELFAC ||
                mrat / meanrat > SC_MAXRELFAC) {
                cpl_table_set(fitpar, "value", row, meanrat);
                cpl_table_set(fitpar, "err_est", row, SC_DEFERRVAL);
            }
        }
    }

    /* Set group fit flag to 0 if rel. RMS < fitlim and/or no line peaks can
       be fitted */
    for (row = 0; row < ngroup[0] + ngroup[1]; row++) {
        mrat = cpl_table_get(fitpar, "value", row, NULL);
        msig = cpl_table_get(fitpar, "err_est", row, NULL);
        if (mrat == 0 || msig == SC_DEFERRVAL || msig / mrat < fitlim) {
            cpl_table_set(fitpar, "fit", row, 0);
        } else {
            cpl_table_set(fitpar, "fit", row, 1);
        }
    }

    /* Free memory */
    cpl_array_delete(tmpratio);
    cpl_array_delete(ratio);
    cpl_array_delete(corrsyst);
    cpl_array_delete(corrfac);
    cpl_array_delete(meanratio);
    cpl_array_delete(groupsum);
    cpl_array_delete(maxrms[0]);
    cpl_array_delete(maxrms[1]);
    cpl_table_delete(tmpfitpar);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_mpfit_setpar(scpars *fitpars, const cpl_table *fitpar,
                               const char fittype)
{
    /*!
     * Provides a vector of parameters which are variables of the CMPFIT
     * fitting process. Moreover, constraints for these parameters are
     * delivered by the CMPFIT structure mp_par. Both objects are returned by
     * a container structure. The parameter information is taken from the
     * fit parameter table fitpar. The input parameter \e fittype decides
     * whether line groups ('l') or the wavelength grid ('w') are fitted.
     *
     * \b INPUT:
     * \param fitpar   CPL table of fit parameters
     * \param fittype  'l' for line groups or 'w' for wavelength grid
     *
     * \b OUTPUT:
     * \param fitpars  structure containing the fit parameters
     *
     * \b ERRORS:
     * - Inconsistent data grids
     */

    char errtxt[SC_MAXLEN], name[SC_LENLINE+1], type[SC_LENLINE+1];
    int npar = 0, i = 0, fit = 0;

    /* Get and check number of parameters */
    npar = cpl_table_get_nrow(fitpar);
    if (fitpars->n != npar) {
        sprintf(errtxt, "%s: # of parameters: "
                "scpars *fitpars != cpl_table *fitpar", SC_ERROR_IDG_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IDG, "%s", errtxt);
    }

    /* Set parameter values and constraints for fitting */

    for (i = 0; i < npar; i++) {
        /* Type of parameter */
        sprintf(type, "%s", cpl_table_get_string(fitpar, "type", i));
        /* Parameter value */
        fitpars->p[i] = cpl_table_get(fitpar, "value", i, NULL);
        /* Fit flag */
        fit = cpl_table_get(fitpar, "fit", i, NULL);
        /* Parameter limits (0 -> no limits) and relative steps */
        if (type[0] == 'w') {
            /* Coefficients for wavelength correction */
            if (fittype == 'w') {
                fitpars->pars[i].fixed = 1 - fit;
            } else {
                fitpars->pars[i].fixed = 1;
            }
            fitpars->pars[i].limited[0] = 0;
            fitpars->pars[i].limited[1] = 0;
            fitpars->pars[i].limits[0] = 0.;
            fitpars->pars[i].limits[1] = 0.;
            fitpars->pars[i].relstep = 0.01;
        } else {
            /* Line group flux correction factors */
            if (fittype == 'l') {
                fitpars->pars[i].fixed = 1 - fit;
            } else {
                fitpars->pars[i].fixed = 1;
            }
            fitpars->pars[i].limited[0] = 1;
            fitpars->pars[i].limited[1] = 1;
            fitpars->pars[i].limits[0] = SC_CORRFAC_MIN - SC_TOL;
            fitpars->pars[i].limits[1] = SC_CORRFAC_MAX + SC_TOL;
            fitpars->pars[i].relstep = 0.01;
        }
        /* Parameter label */
        sprintf(name, "%c%d", type[0],
                cpl_table_get_int(fitpar, "N", i, NULL));
        strcpy(fitpars->pars[i].parname, name);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_mpfit_allocmempar(scpars *fitpars, const int npar)
{
    /*!
     * Allocates memory for an ::scpars structure.
     *
     * \b INPUT:
     * \param npar     number of fit parameters
     *
     * \b OUTPUT:
     * \param fitpars  ::scpars structure containing npar fit parameters
     *
     * \b ERRORS:
     * - Insufficient memory
     */

    cpl_boolean fl_mem = CPL_TRUE;
    char errtxt[SC_MAXLEN];
    int it = 0, i = 0, nchar = SC_LENLINE+1;

    /* Number of fit parameters */
    fitpars->n = npar;

    /* Allocate memory for parameter vector */
    fitpars->p = (double *) calloc(fitpars->n, sizeof(double));
    if (fitpars->p == NULL) {
        fitpars->n = 0;
        sprintf(errtxt, "%s: scpars *fitpars", SC_ERROR_ISM_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_ISM, "%s", errtxt);
    }

    /* Allocate memory for parameter constraints structure */
    fitpars->pars = (mp_par *) calloc(fitpars->n, sizeof(mp_par));
    if (fitpars->pars == NULL) {
        fitpars->n = 0;
        free(fitpars->p);
        fitpars->p = NULL;
        sprintf(errtxt, "%s: scpars *fitpars", SC_ERROR_ISM_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_ISM, "%s", errtxt);
    }

    /* Allocate memory for parameter names in parameter constraints
       structure */
    for (it = 0; it < 2; it++) {
        for (i = 0; i < fitpars->n; i++) {
            fitpars->pars[i].parname = (char *) calloc(nchar, sizeof(char));
            if (it == 0 && fitpars->pars[i].parname == NULL) {
                nchar = 0;
                fl_mem = CPL_FALSE;
                continue;
            }
        }
        if (it == 0 && fl_mem == CPL_TRUE) {
            break;
        } else if (it == 1 && fl_mem == CPL_FALSE) {
            fitpars->n = 0;
            free(fitpars->p);
            fitpars->p = NULL;
            free(fitpars->pars);
            fitpars->pars = NULL;
            sprintf(errtxt, "%s: scpars *fitpars", SC_ERROR_ISM_TXT);
            return cpl_error_set_message(cpl_func, SC_ERROR_ISM, "%s",
                                         errtxt);
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_mpfit_freemempar(scpars *fitpars)
{
    /*!
     * Frees memory occupied by an ::scpars structure.
     *
     * \b INPUT:
     * \param fitpars  ::scpars structure with fit parameters
     *
     * \b OUTPUT:
     * \param fitpars  ::scpars structure without allocated memory
     *
     * \b ERRORS:
     * - none
     */

    int i;

    /* Free memory occupied by parameter vector */
    if (fitpars->p != NULL) {
        free(fitpars->p);
        fitpars->p = NULL;
    }

    /* Free memory occupied by parameter constraints structure */
    if (fitpars->pars != NULL) {
        for (i = 0; i < fitpars->n; i++) {
            free(fitpars->pars[i].parname);
            fitpars->pars[i].parname = NULL;
        }
        free(fitpars->pars);
        fitpars->pars = NULL;
    }

    /* Set number of parameters to 0 */
    fitpars->n = 0;

    return CPL_ERROR_NONE;
}


cpl_error_code sc_mpfit_allocmemresult(mp_result *result, const int m,
                                       const int n)
{
    /*!
     * Allocates memory for fit residuals and parameter errors in the CMPFIT
     * results structure.
     *
     * \b INPUT:
     * \param m       number of data points
     * \param n       number of parameters
     *
     * \b OUTPUT:
     * \param result  CMPFIT structure for fit results with allocated memory
     *                for residuals and parameter errors
     *
     * \b ERRORS:
     * - Insufficient memory
     */

    char errtxt[SC_MAXLEN];

    /* No consideration of the covariance matrix */
    result->covar = NULL;

    /* Memory allocation for fit residuals */
    result->resid = (double *) calloc(m, sizeof(double));
    if (result->resid == NULL) {
        result->status = -99;
        result->xerror = NULL;
        sprintf(errtxt, "%s: mp_result *result", SC_ERROR_ISM_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_ISM, "%s", errtxt);
    }

    /* Memory allocation for parameter errors */
    result->xerror = (double *) calloc(n, sizeof(double));
    if (result->xerror == NULL) {
        result->status = -99;
        free(result->resid);
        result->resid = NULL;
        sprintf(errtxt, "%s: mp_result *result", SC_ERROR_ISM_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_ISM, "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_mpfit_initresult(mp_result *result, const int m,
                                   const int n)
{
    /*!
     * Initialises CMPFIT structure for fit results.
     *
     * \b INPUT:
     * \param result  CMPFIT structure for fit results
     * \param m       number of data points
     * \param n       number of parameters
     *
     * \b OUTPUT:
     * \param result  initialised CMPFIT structure for fit results
     *
     * \b ERRORS:
     * - none
     */

    int i = 0, j = 0;

    /* Set default values */
    result->bestnorm = HUGE_VAL;
    result->orignorm = HUGE_VAL;
    result->niter = 0;
    result->nfev = 0;
    result->status = 99;
    result->npar = n;
    result->nfree = 0;
    result->npegged = 0;
    result->nfunc = m;
    for (i = 0; i < result->nfunc; i++) {
        result->resid[i] = 0;
    }
    for (j = 0; j < result->npar; j++) {
        result->xerror[j] = 0;
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_mpfit_copyresult(mp_result *outresult,
                                   const mp_result *inresult)
{
    /*!
     * Copies content of CMPFIT structure for fit results.
     *
     * \b INPUT:
     * \param inresult   CMPFIT structure for fit results
     *
     * \b OUTPUT:
     * \param outresult  copied CMPFIT structure for fit results
     *
     * \b ERRORS:
     * - none
     */

    int i = 0, j = 0;

    /* Copy data */
    outresult->bestnorm = inresult->bestnorm;
    outresult->orignorm = inresult->orignorm;
    outresult->niter = inresult->niter;
    outresult->nfev = inresult->nfev;
    outresult->status = inresult->status;
    outresult->npar = inresult->npar;
    outresult->nfree = inresult->nfree;
    outresult->npegged = inresult->npegged;
    outresult->nfunc = inresult->nfunc;
    for (i = 0; i < outresult->nfunc; i++) {
        outresult->resid[i] = inresult->resid[i];
    }
    for (j = 0; j < outresult->npar; j++) {
        outresult->xerror[j] = inresult->xerror[j];
    }
    outresult->covar = inresult->covar;
    strcpy(outresult->version, inresult->version);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_mpfit_freememresult(mp_result *result)
{
    /*!
     * Frees memory occupied by a CMPFIT structure for fit residuals and
     * parameter errors. Memory allocated for carrying the covariance matrix
     * is not freed.
     *
     * \b INPUT:
     * \param result  CMPFIT structure for fit results with allocated memory
     *                for residuals and parameter errors
     *
     * \b OUTPUT:
     * \param result  CMPFIT structure for fit results without memory for
     *                residuals and parameter errors
     *
     * \b ERRORS:
     * - none
     */

    /* Free memory for residuals */
    if (result->resid != NULL) {
        free(result->resid);
        result->resid = NULL;
    }

    /* Free memory for parameter errors */
    if (result->xerror != NULL) {
        free(result->xerror);
        result->xerror = NULL;
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_mpfit_substbadfitpar(cpl_table *fitpar,
                                       const cpl_table *initfitpar,
                                       const cpl_parameterlist *parlist)
{
    /*!
     * Substitutes bad fits of flux correction factors by value of initial
     * estimate. Bad fits are recognised by a relative uncertainty greater
     * than ::SC_MAXPARERR. If the A group fits of a system are identified as
     * bad, then all B group fits of this system are also rejected. On the
     * other hand, B groups are not touched if the A group fits of the same
     * system are sufficiently good. The fit uncertainties of the substituted
     * parameters are set to ::SC_DEFERRVAL.
     *
     * \b INPUT:
     * \param fitpar      CPL table with best fit parameters
     * \param initfitpar  CPL table with initial estimates of fit parameters
     * \param parlist     general CPL parameter list
     *
     * \b OUTPUT:
     * \param fitpar      fit parameter table with initial estimates and no
     *                    fit errors for substituted parameters
     *
     * \b ERRORS:
     * - none
     */

    const cpl_parameter *p;
    cpl_array *minrat;
    int ngroupa = 0, ngroupb = 0, ngroup = 0, i = 0, syst = 0;
    double val = 0., err = 0., rat = 0.;

    /* Get number of line groups */
    p = cpl_parameterlist_find_const(parlist, "n_groupA");
    ngroupa = cpl_parameter_get_int(p);
    p = cpl_parameterlist_find_const(parlist, "n_groupB");
    ngroupb = cpl_parameter_get_int(p);
    ngroup = ngroupa + ngroupb;

    /* Create array for minimum errors of systems and initialise it */
    p = cpl_parameterlist_find_const(parlist, "n_system");
    minrat = cpl_array_new(cpl_parameter_get_int(p), CPL_TYPE_DOUBLE);
    cpl_array_fill_window(minrat, 0, cpl_parameter_get_int(p), HUGE_VAL);

    /* Find bad A group parameter fits and substitute them by initial
       estimate */

    for (i = 0; i < ngroupa; i++) {

        if (cpl_table_get(fitpar, "fit", i, NULL) == 1) {

            /* Derive ratio of error and flux */
            val = cpl_table_get(fitpar, "value", i, NULL);
            err = cpl_table_get(fitpar, "err_fit", i, NULL);
            if (val <= 0) {
                rat = HUGE_VAL;
            } else {
                rat = err / val;
            }

            /* Identify bad parameters and modify parameter values and
               errors */
            if (rat > SC_MAXPARERR && err != SC_DEFERRVAL) {
                cpl_table_set(fitpar, "value", i,
                              cpl_table_get(initfitpar, "value", i, NULL));
                cpl_table_set(fitpar, "err_fit", i, SC_DEFERRVAL);
            }

            /* Find minimum ratio for each group system */
            syst = cpl_table_get(fitpar, "system", i, NULL);
            if (syst > 0 && rat < cpl_array_get(minrat, syst - 1, NULL)) {
                cpl_array_set(minrat, syst - 1, rat);
            }

        }

    }

    /* Substitute B group parameter fits if the minimum ratio of a system is
       higher than SC_MAXPARERR */

    for (i = ngroupa; i < ngroup; i++) {

        if (cpl_table_get(fitpar, "fit", i, NULL) == 1) {

            /* Identify parameters of bad systems and modify parameter values
               and errors */
            syst = cpl_table_get(fitpar, "system", i, NULL);
            if (syst > 0 &&
                cpl_array_get(minrat, syst - 1, NULL) > SC_MAXPARERR) {
                cpl_table_set(fitpar, "value", i,
                              cpl_table_get(initfitpar, "value", i, NULL));
                cpl_table_set(fitpar, "err_fit", i, SC_DEFERRVAL);
            }

        }

    }

    /* Free memory */
    cpl_array_delete(minrat);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_mpfit_writeresults(const mp_result *result,
                                     const cpl_table *scispec,
                                     const cpl_table *fitpar,
                                     const cpl_parameterlist *parlist,
                                     const int calls, const double runtime)
{
    /*!
     * Writes a summary of the CMPFIT results into an ASCII file.
     *
     * \b INPUT:
     * \param result   CMPFIT structure for fit results
     * \param scispec  CPL table with science and best-fit sky spectrum
     * \param fitpar   CPL table with best fit parameters
     * \param parlist  general CPL parameter list
     * \param calls    number of CMPFIT calls
     * \param runtime  fit run time
     *
     * \b OUTPUT:
     * - none
     *
     * \b ERRORS:
     * - File opening failed
     * - Insufficient data points
     * - Invalid object value(s)
     */

    FILE *stream;
    cpl_error_code err = CPL_ERROR_NONE;
    const cpl_parameter *p;
    cpl_table *tmpspec;
    cpl_array *ratio;
    char basedir[SC_MAXLEN], output_dir[SC_MAXLEN];
    char output_name[SC_MAXLEN], outfile[SC_MAXLEN];
    char errtxt[SC_MAXLEN], type[SC_LENLINE+1];
    int npix = 0, i = 0, nw = 0, npar = 0, j = 0, nrelpar = 0, dof = 0;
    int sum = 0;
    double chi2red = 0., wf2sum = 0., devf2sum = 0., wsum = 0., wfsum = 0.;
    double w2sum = 0., dev2sum = 0., wl2sum = 0., devl2sum = 0., arms = 0.;
    double wmean = 0., frms = 0., lrms = 0., rms = 0., lflux = 0.;
    double mlflux = 0., rat = 0., mrat = 0., sig = 0., fwhm = 0.;
    const double *weight, *flux, *dev;

    /* Get output file path and name */
    p = cpl_parameterlist_find_const(parlist, "inst_dir");
    strcpy(basedir, cpl_parameter_get_string(p));
    p = cpl_parameterlist_find_const(parlist, "output_dir");
    sc_basic_abspath(output_dir, cpl_parameter_get_string(p), basedir);
    p = cpl_parameterlist_find_const(parlist, "output_name");
    strcpy(output_name, cpl_parameter_get_string(p));
    sprintf(outfile, "%s%s_fit.res", output_dir, output_name);

    /* Open output file */
    if ((stream = fopen(outfile, "w+")) == NULL) {
        sprintf(errtxt, "%s: %s", SC_ERROR_FOF_TXT, outfile);
        return cpl_error_set_message(cpl_func, SC_ERROR_FOF, "%s", errtxt);
    }

    /* Write file names of science and sky spectrum */
    fprintf(stream, "INPUT DATA FILES:\n");
    p = cpl_parameterlist_find_const(parlist, "scispec");
    fprintf(stream, "Science: %s\n", cpl_parameter_get_string(p));
    p = cpl_parameterlist_find_const(parlist, "skyspec");
    fprintf(stream, "Sky:     %s\n\n", cpl_parameter_get_string(p));

    /* Write CMPFIT status */
    fprintf(stream, "MPFIT RESULTS:\n");
    fprintf(stream, "Status:                 %d\n", result->status);

    /* Get number of data points with non-zero weight (= DOF+1) */

    npix = cpl_table_get_nrow(scispec);
    for (i = 0; i < npix; i++) {
        if (cpl_table_get(scispec, "cweight", i, NULL) > 0) {
            nw++;
        }
    }

    if (nw == 1) {
        sprintf(errtxt, "%s: cpl_table *scispec "
                "(only 1 data point with weight > 0)", SC_ERROR_ISD_TXT);
        //err = cpl_error_set_message(cpl_func, SC_ERROR_ISD, "%s", errtxt);
        cpl_msg_warning(cpl_func, "cpl_table *scispec (only 1 data point "
                        "with weight > 0)");
    }

    /* Get number of relevant fit parameters */
    npar = cpl_table_get_nrow(fitpar);
    for (j = 0; j < npar; j++) {
        if (cpl_table_get(fitpar, "relevance", j, NULL) == 1) {
            nrelpar++;
        }
    }

    /* Write CMPFIT results */

    fprintf(stream, "Fit parameters:         %d\n", nrelpar);
    fprintf(stream, "Data points:            %d\n", npix);
    fprintf(stream, "Weight > 0:             %d\n", nw);
    fprintf(stream, "MPFIT calls:            %d\n", calls);
    fprintf(stream, "Iterations:             %d\n", result->niter);
    fprintf(stream, "Function evaluations:   %d\n", nfev);
    fprintf(stream, "Fit run time in s:      %.2f\n", runtime);
    fprintf(stream, "Initial chi2:           %.3e\n", result->orignorm);
    fprintf(stream, "Best chi2:              %.3e\n", result->bestnorm);
    if (nw - nrelpar <= 1) {
        fprintf(stream, "Reduced chi2:           UNDEF\n");
        fprintf(stream, "RMS rel. to error:      UNDEF\n");
    } else {
        dof = nw - nrelpar - 1;
        chi2red = result->bestnorm / (double) dof;
        fprintf(stream, "Reduced chi2:           %.3e\n", chi2red);
        fprintf(stream, "RMS rel. to error:      %.3e\n", sqrt(chi2red));
    }

    /* Save results and set sigclip to 0 to calculate deviations for all
       pixels */
    tmpspec = cpl_table_duplicate(scispec);
    cpl_table_fill_column_window(tmpspec, "sigclip", 0, npix, 0);
    sc_modsky_calcdev(tmpspec);

    /* Get pointers to table columns */
    flux = cpl_table_get_data_double_const(tmpspec, "lflux");
    weight = cpl_table_get_data_double_const(tmpspec, "cweight");
    dev = cpl_table_get_data_double_const(tmpspec, "dev");

    /* Compute weighted deviations for line peaks, lines (peak +- HWHM),
       and the full spectrum and weighted line peak flux */
    for (i = 0; i < npix; i++) {
        wf2sum += weight[i] * weight[i];
        devf2sum += dev[i] * dev[i];
        if (cpl_table_get(scispec, "sigclip", i, NULL) == 0) {
            if (cpl_table_get(scispec, "class", i, NULL) >= 2) {
                wsum += weight[i];
                wfsum += weight[i] * flux[i];
                w2sum += weight[i] * weight[i];
                dev2sum += dev[i] * dev[i];
            }
            wl2sum += weight[i] * weight[i];
            devl2sum += dev[i] * dev[i];
        }
    }

    /* Delete temporary table */
    cpl_table_delete(tmpspec);

    /* Compute RMS of full spectrum and write it to file */
    if (wf2sum == 0) {
        fprintf(stream, "Full RMS:               UNDEF\n");
    } else {
        arms = sqrt(devf2sum / wf2sum);
        fprintf(stream, "Full RMS:               %.3e\n", arms);
    }

    /* Compute three different RMS relative to weighted mean of line peaks and
       write them to file */
    if (wsum == 0 || wfsum == 0) {
        fprintf(stream, "Full RMS rel. to peaks: UNDEF\n");
        fprintf(stream, "Line RMS rel. to peaks: UNDEF\n");
        fprintf(stream, "Peak RMS rel. to peaks: UNDEF\n");
        if (wsum == 0) {
            sprintf(errtxt, "%s: cpl_table *scispec (all weights = 0)",
                    SC_ERROR_IOV_TXT);
            //err = cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s",
            //errtxt);
            cpl_msg_warning(cpl_func, "cpl_table *scispec (all weights = 0)");
        } else {
            sprintf(errtxt, "%s: cpl_table *scispec (all fluxes = 0)",
                    SC_ERROR_IOV_TXT);
            err = cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
        }
    } else {
        wmean = wfsum / wsum;
        frms = arms / wmean;
        fprintf(stream, "Full RMS rel. to peaks: %.3e\n", frms);
        lrms = sqrt(devl2sum / wl2sum) / wmean;
        fprintf(stream, "Line RMS rel. to peaks: %.3e\n", lrms);
        rms = sqrt(dev2sum / w2sum) / wmean;
        fprintf(stream, "Peak RMS rel. to peaks: %.3e\n", rms);
    }

    /* Compute sigma-clipped mean of ratio between residual and input line
       flux for valid line peaks and write it to file */

    ratio = cpl_array_new(npix, CPL_TYPE_DOUBLE);

    for (sum = 0, i = 0; i < npix; i++) {
        if (cpl_table_get(scispec, "sigclip", i, NULL) == 0 &&
            cpl_table_get(scispec, "class", i, NULL) >= 2) {
            /* Valid (i.e. unclipped) line peak found in science spectrum */
            lflux = cpl_table_get(scispec, "lflux", i, NULL);
            mlflux = cpl_table_get(scispec, "mlflux", i, NULL);
            if (lflux == 0) {
                rat = 1e6;
            } else {
                rat = fabs((mlflux - lflux) / lflux);
            }
            cpl_array_set_double(ratio, i, rat);
            sum++;
        } else {
            /* Set array element invalid if pixel is not valid line peak */
            cpl_array_set_invalid(ratio, i);
        }
    }

    /* Perform sigma-clipping for entire spectrum */
    if (sum > 0) {
        sc_basic_clipmean(&mrat, &sig, ratio, CPL_FALSE);
    }

    cpl_array_delete(ratio);

    if (sum == 0) {
        fprintf(stream, "Mean rel. residual:     UNDEF\n\n");
        sprintf(errtxt, "%s: cpl_table *scispec (no valid line peaks)",
                SC_ERROR_IOV_TXT);
        //err = cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
        cpl_msg_warning(cpl_func, "cpl_table *scispec (no valid line peaks)");
    } else {
        fprintf(stream, "Mean rel. residual:     %.3e\n\n", mrat);
    }

    /* FWHM from parameter list */
    p = cpl_parameterlist_find_const(parlist, "fwhm");
    fwhm = cpl_parameter_get_double(p);

    /* Spectral resolution from line width measurements in pixels */
    fprintf(stream, "ESTIMATED SPECTRAL RESOLUTION:\n");
    fprintf(stream, "FWHM in pixels: %.3f\n\n", fwhm);

    /* Write best-fit parameters */
    fprintf(stream, "BEST-FIT PARAMETERS:\n\n");

    /* Write parameter type, number, and best-fit value with uncertainty */
    fprintf(stream,
            "Type N  Fit     Value                RMS        N_lin\n");

    /* Loop over all parameters but write relvant parameters only */
    for (j = 0; j < npar; j++) {
        if (cpl_table_get_int(fitpar, "relevance", j, NULL) == 1) {
            sprintf(type, "%s", cpl_table_get_string(fitpar, "type", j));
            fprintf(stream, "%c    %.2d  %d  %10.4g +- %-9.4g  %-9.4g  %d\n",
                    type[0], cpl_table_get_int(fitpar, "N", j, NULL),
                    cpl_table_get_int(fitpar, "fit", j, NULL),
                    cpl_table_get(fitpar, "value", j, NULL),
                    cpl_table_get(fitpar, "err_fit", j, NULL),
                    cpl_table_get(fitpar, "err_est", j, NULL),
                    cpl_table_get_int(fitpar, "N_lin", j, NULL));
        }
    }

    /* Explanations for printed table */
    fprintf(stream, "\nREMARKS:\n");
    fprintf(stream, "Type: A/B = line group A/B, w = wavelength fit coef.\n");
    fprintf(stream, "Fit: 1 = free MPFIT par., 0 = only initial estimate\n");
    fprintf(stream, "RMS: uncertainty of initial estimate\n");
    fprintf(stream, "%g: no error available\n", SC_DEFERRVAL);
    fprintf(stream, "N_lin: number of lines for fitting \n");

    fclose(stream);

    return err;
}

/**@}*/
