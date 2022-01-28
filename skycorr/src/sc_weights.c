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
 * \file sc_weights.c
 *
 * Routines for preparing the line group weights of each pixel of the input
 * sky spectrum
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  14 Feb 2011
 * \date   19 Sep 2013
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <sc_weights.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

cpl_error_code sc_weights(cpl_table *skyspec, cpl_table *fitpar,
                          cpl_table *groups, cpl_parameterlist *parlist)
{
    /*!
     * \callgraph
     *
     * Derives line group weights for each pixel of the input sky spectrum.
     *
     * \b INPUT:
     * \param skyspec  CPL table with sky spectrum
     * \param groups   CPL table with airglow line information
     * \param parlist  input CPL parameter list
     *
     * \b OUTPUT:
     * \param skyspec  sky spectrum with line group weights
     * \param fitpar   CPL table of fit parameters
     * \param groups   modified line table
     * \param parlist  parameter list complemented by weight-specific
     *                 parameters
     *
     * \b ERRORS:
     * - none
     */

    /* Create table of fit parameters */
    sc_weights_initfitpar(fitpar, parlist, groups);

    /* Write info message */
    cpl_msg_info(cpl_func, "Derive line group weights for each pixel "
                 "of the sky spectrum");

    /* Calculate spectrum for each line group and get group contributions to
       each pixel of the input sky spectrum */
    sc_weights_getpixcontrib(skyspec, fitpar, groups, parlist);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_weights_initfitpar(cpl_table *fitpar,
                                     cpl_parameterlist *parlist,
                                     const cpl_table *groups)
{
    /*!
     * Creates table of fit parameters. There are table columns for parameter
     * type, ID number, system flag, relevance flag, fit flag, value, fit
     * error, estimate error, and number of lines. The parameter types are
     * either the line groups 'A' and 'B' or 'w' for the coefficients of the
     * Chebyshev polynomial for the wavelength grid correction. The number N
     * is a counter that starts with 1 in the case of the line groups and with
     * 0 in the case of the polynomial and the \f${\sigma}\f$-clipping limit.
     * The system flag indicates fit parameters that are connected (e.g. A and
     * B groups of OH lines). By default the relevance and fit flags are set
     * to 1 (valid fit parameter) for all parameters. The initial parameter
     * value is 1 for the line groups. For the 'w' type the coefficient with
     * N = 1 is set to 1 and the remaining ones become 0. This setup ensures
     * that there is no wavelength correction by default. The number of
     * wavelength coefficients is given by the driver parameter \e cheby_max.
     * An exception is a selected degree of 0 for which also the linear
     * coefficient has to be provided to reproduce the input wavelength grid.
     * For the coefficients of the polynomial the initial relevance and fit
     * flags are set to 0. The fit and estimate error columns as well as the
     * number of lines column are initialised with 0 for all parameters.
     * Finally, the general CPL parameter list is supplemented by the maximum
     * N of the parameter types 'A' and 'B' and the number of parameter
     * systems.
     *
     * \b INPUT:
     * \param parlist  input CPL parameter list
     * \param groups   CPL table with airglow line information
     *
     * \b OUTPUT:
     * \param fitpar   CPL table of fit parameters
     * \param parlist  parameter list complemented by line group counts
     *
     * \b ERRORS:
     * - none
     */

    cpl_parameter *p;
    int nlin = 0, i = 0, groupa = 0, groupb = 0, maxa = 0, maxb = 0, minb = 0;
    int deg = 0, ncoef = 0, npar = 0, j = 0, group = 0, syst = 0;

    /* Create CPL table columns */
    cpl_table_new_column(fitpar, "type", CPL_TYPE_STRING);
    cpl_table_new_column(fitpar, "N", CPL_TYPE_INT);
    cpl_table_new_column(fitpar, "system", CPL_TYPE_INT);
    cpl_table_new_column(fitpar, "relevance", CPL_TYPE_INT);
    cpl_table_new_column(fitpar, "fit", CPL_TYPE_INT);
    cpl_table_new_column(fitpar, "value", CPL_TYPE_DOUBLE);
    cpl_table_new_column(fitpar, "err_fit", CPL_TYPE_DOUBLE);
    cpl_table_new_column(fitpar, "err_est", CPL_TYPE_DOUBLE);
    cpl_table_new_column(fitpar, "N_lin", CPL_TYPE_INT);

    /* Get number of airglow lines */
    nlin = cpl_table_get_nrow(groups);

    /* Find maximum group numbers for both group types and minimum number for
       B groups (= number of systems) */
    for (i = 0; i < nlin; i++) {
        groupa = cpl_table_get(groups, "groupA", i, NULL);
        groupb = cpl_table_get(groups, "groupB", i, NULL);
        if (groupa > maxa) {
            maxa = groupa;
        }
        if (groupb > maxb) {
            maxb = groupb;
        }
        if (groupb < minb) {
            minb = groupb;
        }
    }

    /* Write group numbers into general parameter list */
    p = cpl_parameter_new_value("n_groupA", CPL_TYPE_INT, "", "", maxa);
    cpl_parameterlist_append(parlist, p);
    p = cpl_parameter_new_value("n_groupB", CPL_TYPE_INT, "", "", maxb);
    cpl_parameterlist_append(parlist, p);
    p = cpl_parameter_new_value("n_system", CPL_TYPE_INT, "", "", -minb);
    cpl_parameterlist_append(parlist, p);

    /* Get maximum number of coefficients of the Chebyshev polynomial for the
       wavelength grid correction */
    p = cpl_parameterlist_find(parlist, "cheby_max");
    deg = cpl_parameter_get_int(p);
    if (deg < 0) {
        ncoef = 0;
    } else if (deg == 0) {
        /* Constant term is not sufficient to reproduce input wavelength
           grid */
        ncoef = 2;
    } else {
        ncoef = deg + 1;
    }

    /* Resize fit parameter table */
    npar = maxa + maxb + ncoef;
    cpl_table_set_size(fitpar, npar);

    /* Put parameter type, number, and initial value into table */

    for (j = 0; j < npar; j++) {

        if (j >= maxa + maxb) {
            cpl_table_set_string(fitpar, "type", j, "w");
            if (j == maxa + maxb) {
                p = cpl_parameterlist_find(parlist, "cheby_const");
                cpl_table_set(fitpar, "value", j,
                              cpl_parameter_get_double(p));
            } else if (j == maxa + maxb + 1) {
                cpl_table_set(fitpar, "value", j, 1.);
            } else {
                cpl_table_set(fitpar, "value", j, 0.);
            }
        } else {
            if (j < maxa) {
                cpl_table_set_string(fitpar, "type", j, "A");
            } else {
                cpl_table_set_string(fitpar, "type", j, "B");
            }
            cpl_table_set(fitpar, "value", j, 1.);
        }

        if (j == maxa + maxb) {
            group = 0;
        } else if (j == maxa) {
            group = 1;
        } else {
            group++;
        }
        cpl_table_set_int(fitpar, "N", j, group);

    }

    /* Set all relevance and fit flags to 1 for 'A' and 'B' parameters and to
       0 for 'w' parameters */
    cpl_table_fill_column_window(fitpar, "relevance", 0, maxa + maxb, 1);
    cpl_table_fill_column_window(fitpar, "fit", 0, maxa + maxb, 1);
    cpl_table_fill_column_window(fitpar, "relevance", maxa + maxb, ncoef, 0);
    cpl_table_fill_column_window(fitpar, "fit", maxa + maxb, ncoef, 0);

    /* Get system flags from line list */
    cpl_table_fill_column_window(fitpar, "system", 0, npar, 0);
    for (i = 0; i < nlin; i++) {
        syst = cpl_table_get(groups, "system", i, NULL);
        groupa = cpl_table_get(groups, "groupA", i, NULL);
        groupb = cpl_table_get(groups, "groupB", i, NULL);
        cpl_table_set_int(fitpar, "system", groupa - 1, syst);
        if (groupb >= 0) {
            cpl_table_set_int(fitpar, "system", maxa + groupb - 1, syst);
        }
    }

    /* Set the two error columns to 0 */
    cpl_table_fill_column_window(fitpar, "err_fit", 0, npar, 0.);
    cpl_table_fill_column_window(fitpar, "err_est", 0, npar, 0.);

    /* Set number of lines for each group to 0 */
    cpl_table_fill_column_window(fitpar, "N_lin", 0, npar, 0);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_weights_getpixcontrib(cpl_table *skyspec, cpl_table *fitpar,
                                        cpl_table *groups,
                                        const cpl_parameterlist *parlist)
{
    /*!
     * \callgraph
     *
     * Calculates spectrum for each line group and gets group contributions to
     * each pixel of the input sky spectrum. Depending on the parameter
     * \e varfwhm, the line widths are constant or linearly-increasing with
     * wavelength.
     *
     * \b INPUT:
     * \param skyspec    CPL table with sky spectrum
     * \param fitpar     CPL table of fit parameters
     * \param groups     CPL table with airglow line information
     * \param parlist    input CPL parameter list
     *
     * \b OUTPUT:
     * \param skyspec    sky spectrum with line group weights
     * \param fitpar     table of fit parameters with adapted fit flags
     * \param groups     line list for wavelength range of sky spectrum only
     *
     * \b ERRORS:
     * - none
     */

    const cpl_parameter *p;
    cpl_array *array;
    cpl_table *groupspec;
    char grouptype[2] = "AB", ncolname[SC_LENLINE+1], wcolname[SC_LENLINE+1];
    char dcolname[SC_LENLINE+1], grouppar[SC_LENLINE+1];
    int nrow = 0, varfwhm = 0, i = 0, ngroup = 0, mingroup = 0, group = 0;
    double lcut = 0., ucut = 0., fwhm = 0., speedpar = 0., meanlam = 0.;

    /* Create temporary columns for fluxes of selected group and all groups */
    nrow = cpl_table_get_nrow(skyspec);
    cpl_table_new_column(skyspec, "f_group", CPL_TYPE_DOUBLE);
    cpl_table_new_column(skyspec, "f_all", CPL_TYPE_DOUBLE);

    /* Create wavelength grid for temporary line group spectra */
    groupspec = cpl_table_new(0);
    sc_weights_createwavegrid(groupspec, skyspec);

    /* Erase lines of the line table at wavelengths outside the range of the
       line group spectrum */
    cpl_table_unselect_all(groups);
    lcut = cpl_table_get(groupspec, "minlam", 0, NULL);
    cpl_table_or_selected_double(groups, "lambda", CPL_LESS_THAN, lcut);
    cpl_table_erase_selected(groups);
    cpl_table_unselect_all(groups);
    ucut = cpl_table_get(groupspec, "minlam", cpl_table_get_nrow(groupspec)-1,
                         NULL);
    cpl_table_or_selected_double(groups, "lambda", CPL_GREATER_THAN, ucut);
    cpl_table_erase_selected(groups);
    cpl_table_select_all(groups);

    /* Get FWHM for mean wavelength from parameter list in pixels
       (consider oversampling in group spectrum) */
    p = cpl_parameterlist_find_const(parlist, "fwhm");
    fwhm = cpl_parameter_get_double(p) * SC_SAMPFAC;

    /* Gaussian kernel: linearly increasing with wavelength or constant? */
    p = cpl_parameterlist_find_const(parlist, "varfwhm");
    varfwhm = cpl_parameter_get_int(p);
    if (varfwhm == 1) {
        /* variable kernel */
        speedpar = SC_LIMRELLAMVAR;
    } else {
        /* constant kernel */
        speedpar = HUGE_VAL;
    }

    /* Get mean wavelength */
    p = cpl_parameterlist_find_const(parlist, "meanlam");
    meanlam = cpl_parameter_get_double(p);

    /* Calculate pixel-specific group contributions for each group type */

    for (i = 0; i < 2; i++) {

        /* Initialise temporary flux columns in sky spectrum with 0 */
        cpl_table_fill_column_window(skyspec, "f_group", 0, nrow, 0.);
        cpl_table_fill_column_window(skyspec, "f_all", 0, nrow, 0.);

        /* Create new table columns for line group IDs, weights, and IDs of
           dominating line groups */
        sprintf(ncolname, "ng%c", grouptype[i]);
        cpl_table_new_column_array(skyspec, ncolname, CPL_TYPE_DOUBLE,
                                   SC_COLDEPTH);
        sprintf(wcolname, "wg%c", grouptype[i]);
        cpl_table_new_column_array(skyspec, wcolname, CPL_TYPE_DOUBLE,
                                   SC_COLDEPTH);
        sprintf(dcolname, "dg%c", grouptype[i]);
        cpl_table_new_column(skyspec, dcolname, CPL_TYPE_INT);

        /* Initialise new array table columns with arrays of SC_COLDEPTH
           invalid elements */
        array = cpl_array_new(SC_COLDEPTH, CPL_TYPE_DOUBLE);
        cpl_array_fill_window_invalid(array, 0, SC_COLDEPTH);
        cpl_table_fill_column_window_array(skyspec, ncolname, 0, nrow, array);
        cpl_table_fill_column_window_array(skyspec, wcolname, 0, nrow, array);
        cpl_array_delete(array);

        /* Get number of line groups */
        sprintf(grouppar, "n_group%c", grouptype[i]);
        p = cpl_parameterlist_find_const(parlist, grouppar);
        ngroup = cpl_parameter_get_int(p);

        /* Skip group type if number of groups is zero */
        if (ngroup == 0) {
            continue;
        }

        /* Get minimum group number */
        if (i == 1) {
            /* n_system negative B group numbers */
            p = cpl_parameterlist_find_const(parlist, "n_system");
            mingroup = -cpl_parameter_get_int(p);
        } else {
            mingroup = 0;
        }

        /* Get pixel contributions for each group */

        for (group = mingroup; group <= ngroup; group++) {

            /* Write group lines into temporary spectrum and update group fit
               flags */
            sc_weights_getgrouplines(groupspec, fitpar, groups, parlist,
                                     grouptype[i], group);

            /* Convolve temporary spectrum with Gaussian kernel */
            sc_weights_convolve(groupspec, fwhm, meanlam, speedpar);

            /* Rebin temporary line group spectrum to wavelength grid of input
               sky spectrum */
            sc_weights_rebinspec(skyspec, groupspec);

            /* Write group contributions into pixel-specific arrays */
            sc_weights_fillgrouparrays(skyspec, grouptype[i], group);

        }

        /* Divide group flux by flux sum for each pixel to get weights */
        sc_weights_normpixcontrib(skyspec, grouptype[i]);

        /* Get dominating group for each pixel */
        sc_weights_getdomgroups(skyspec, parlist, grouptype[i]);

    }

    /* Free allocated memory */
    cpl_table_delete(groupspec);

    /* Delete temporary columns */
    cpl_table_erase_column(skyspec, "f_group");
    cpl_table_erase_column(skyspec, "f_all");

    return CPL_ERROR_NONE;
}


cpl_error_code sc_weights_createwavegrid(cpl_table *groupspec,
                                         const cpl_table *skyspec)
{
    /*!
     * Creates a wavelength grid for line group spectra. The original grid of
     * the input sky spectrum is extended by ::SC_EXTRACOVER pixels on each
     * margin and split up in ::SC_SAMPFAC smaller pixels for a better
     * sampling. The output spectrum consists of three columns that provide
     * the central wavelength, the lower wavelength limit, and the flux for
     * each pixel. The latter is always zero.
     *
     * \b INPUT:
     * \param skyspec    CPL table with sky spectrum
     *
     * \b OUTPUT:
     * \param groupspec  CPL table with wavelength grid for line group spectra
     *
     * \b ERRORS:
     * - none
     */

    int nrow = 0, npix = 0, i = 0, j = 0, k = 0;
    double llam = 0., clam = 0., ulam = 0., dlam = 0., llim = 0., lam = 0.;

    /* Create columns of CPL table for line group spectra */
    cpl_table_new_column(groupspec, "lam", CPL_TYPE_DOUBLE);
    cpl_table_new_column(groupspec, "minlam", CPL_TYPE_DOUBLE);
    cpl_table_new_column(groupspec, "flux", CPL_TYPE_DOUBLE);

    /* Get number of pixels by considering extra pixels at the margins and a
       sampling factor that causes a split-up of each intial pixel */
    nrow = cpl_table_get_nrow(skyspec);
    npix = (nrow + 2 * SC_EXTRACOVER) * SC_SAMPFAC + 1;
    cpl_table_set_size(groupspec, npix);

    /* Set values of flux column = 0 */
    cpl_table_fill_column_window(groupspec, "flux", 0, npix, 0.);

    /* Derive wavelength grid from input sky spectrum */

    for (i = 1; i < nrow - 1; i++) {

        /* Get wavelengths of sky spectrum pixels i-1, i, and i+1 */

        if (i == 1) {
            llam = cpl_table_get(skyspec, "lambda", i - 1, NULL);
            clam = cpl_table_get(skyspec, "lambda", i, NULL);
        } else {
            llam = clam;
            clam = ulam;
        }

        ulam = cpl_table_get(skyspec, "lambda", i + 1, NULL);

        /* Handle pixels of lower margin */

        if (i == 1) {

            j = -1;
            dlam = (clam - llam) / SC_SAMPFAC;
            llim = (llam + clam) / 2
                   - dlam * (SC_EXTRACOVER + 1) * SC_SAMPFAC;
            lam = llim + dlam / 2;

            for (k = 0; k < (SC_EXTRACOVER + 1) * SC_SAMPFAC; k++) {
                j++;
                cpl_table_set(groupspec, "lam", j, lam);
                cpl_table_set(groupspec, "minlam", j, llim);
                lam += dlam;
                llim += dlam;
            }

        }

        /* Handle regular pixels */

        dlam = (ulam - llam) / (2 * SC_SAMPFAC);
        llim = (llam + clam) / 2;
        lam = llim + dlam / 2;

        for (k = 0; k < SC_SAMPFAC; k++) {
            j++;
            cpl_table_set(groupspec, "lam", j, lam);
            cpl_table_set(groupspec, "minlam", j, llim);
            lam += dlam;
            llim += dlam;
        }

        /* Handle pixels of upper margin */

        if (i == nrow - 2) {

            dlam = (ulam - clam) / SC_SAMPFAC;
            llim = (clam + ulam) / 2;
            lam = llim + dlam / 2;

            for (k = 0; k < (SC_EXTRACOVER + 1) * SC_SAMPFAC + 1; k++) {
                j++;
                cpl_table_set(groupspec, "lam", j, lam);
                cpl_table_set(groupspec, "minlam", j, llim);
                lam += dlam;
                llim += dlam;
            }

        }

    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_weights_getgrouplines(cpl_table *groupspec,
                                        cpl_table *fitpar, cpl_table *groups,
                                        const cpl_parameterlist *parlist,
                                        const char grouptype, const int group)
{
    /*!
     * Writes line fluxes of a selected group into line group spectrum if the
     * group can be found in the wavelength range covered by the output
     * spectrum. Moreover, for missing groups the relevance and fit flags in
     * the fit parameter table are set to 0.
     *
     * \note The input line table remains in a selected state after the end
     *       of the routine.
     *
     * \b INPUT:
     * \param groupspec  CPL table with wavelength grid for line group spectra
     * \param fitpar     CPL table of fit parameters
     * \param groups     CPL table with lines and group classification
     * \param parlist    input CPL parameter list
     * \param grouptype  line group type (either 'A' or 'B')
     * \param group      number of line group
     *
     * \b OUTPUT:
     * \param groupspec  spectrum with the fluxes of a single line group
     * \param fitpar     table of fit parameters with updated fit flags
     *
     * \b ERRORS:
     * - none
     */

    const cpl_parameter *p;
    char groupcol[SC_LENLINE+1];
    cpl_boolean found = CPL_FALSE;
    int nsel = 0, npix = 0, row = 0, i = 0, j = 0;
    double llim = 0., lam = 0., lflux = 0., ulim = 0., dlam = 0., pflux = 0.;

    /* Select rows that belong to given line group */
    cpl_table_unselect_all(groups);
    sprintf(groupcol, "group%c", grouptype);
    nsel = cpl_table_or_selected_int(groups, groupcol, CPL_EQUAL_TO, group);

    /* Set flux in group spectrum to zero */
    npix = cpl_table_get_nrow(groupspec);
    cpl_table_fill_column_window(groupspec, "flux", 0, npix, 0);

    /* Change fit flag to 0 and return if desired group could not be found */
    if (nsel <= 0) {
        if (grouptype == 'B') {
            p = cpl_parameterlist_find_const(parlist, "n_groupA");
            row = cpl_parameter_get_int(p) + group - 1;
        } else {
            row = group - 1;
        }
        if (group > 0) {
            cpl_table_set(fitpar, "relevance", row, 0);
            cpl_table_set(fitpar, "fit", row, 0);
        }
        return CPL_ERROR_NONE;
    }

    /* Find lines of selected group and enter their fluxes in rows of the
       output table that cover the wavelengths of these lines */
    llim = cpl_table_get(groupspec, "lam", 0, NULL);
    for (i = 0; i < cpl_table_get_nrow(groups); i++) {
        if (cpl_table_is_selected(groups, i) == 1) {
            lam = cpl_table_get(groups, "lambda", i, NULL);
            lflux = cpl_table_get(groups, "flux", i, NULL);
            found = CPL_FALSE;
            do {
                if (j >= npix - 1) break; // strange behaviour of 'while'
                ulim = cpl_table_get(groupspec, "lam", j+1, NULL);
                if (lam >= llim && lam < ulim) {
                    found = CPL_TRUE;
                    dlam = ulim - llim;
                    pflux = cpl_table_get(groupspec, "flux", j, NULL);
                    pflux += lflux / dlam;
                    cpl_table_set(groupspec, "flux", j, pflux);
                } else {
                    llim = ulim;
                    j++;
                }
            } while (found == CPL_FALSE && j < npix - 1);
        }
    }

    /* Cancel row selection */
    cpl_table_select_all(groups);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_weights_convolve(cpl_table *groupspec, const double fwhm,
                                   const double reflam, const double speedpar)
{
    /*!
     * Convolves line group spectrum with a wavelength-dependent Gaussian
     * kernel. The scaling of the kernel width is forced by the assumption of
     * constant resolution. The \e speedpar parameter rules the number of
     * kernel calculations. It provides the relative wavelength change that
     * causes a recalculation. If a constant kernel is desired, \e speedpar
     * has to be set to HUGE_VAL.
     *
     * \b INPUT:
     * \param groupspec  CPL table with input line group spectrum
     * \param fwhm       width of Gaussian in pixels
     * \param reflam     reference wavelength for \e fwhm
     * \param speedpar   criterion for recalculation of kernel
     *                   (relative change of kernel width)
     *
     * \b OUTPUT:
     * \param groupspec  convolved spectrum
     *
     * \b ERRORS:
     * - No data
     * - Invalid input parameter(s)
     * - Invalid object structure
     */

    cpl_array *flux = NULL, *convflux = NULL, *kernel = NULL;
    char errtxt[SC_MAXLEN];
    int m = 0, range[2] = {0, 0}, i = 0;
    double llam = 0., ulam = 0., clam = 0., cfwhm = 0.;
    double *vec, *lam;

    /* Check number of data points in spectrum */
    m = cpl_table_get_nrow(groupspec);
    if (m <= 0) {
        sprintf(errtxt, "%s: cpl_table *groupspec", SC_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_NDA, "%s", errtxt);
    }

    /* Check for invalid input parameters */
    if (fwhm < 0.) {
        sprintf(errtxt, "%s: fwhm < 0", SC_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IIP, "%s", errtxt);
    }
    if (reflam <= 0.) {
        sprintf(errtxt, "%s: reflam <= 0", SC_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IIP, "%s", errtxt);
    }
    if (speedpar < 0.) {
        sprintf(errtxt, "%s: speedpar < 0", SC_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IIP, "%s", errtxt);
    }

    /* Check existence of 'lambda' and 'flux' columns */
    if (cpl_table_has_column(groupspec, "lam") != 1) {
        sprintf(errtxt, "%s: cpl_table *groupspec (no 'lam' column)",
                SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }
    if (cpl_table_has_column(groupspec, "flux") != 1) {
        sprintf(errtxt, "%s: cpl_table *groupspec (no 'flux' column)",
                SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Copy flux column into CPL array */
    vec = cpl_table_get_data_double(groupspec, "flux");
    flux = cpl_array_new(m, CPL_TYPE_DOUBLE);
    cpl_array_copy_data_double(flux, vec);

    /* Create CPL array for convolved flux */
    convflux = cpl_array_new(m, CPL_TYPE_DOUBLE);
    cpl_array_fill_window_double(convflux, 0, m, 0.);

    /* Get pointer to wavelength */
    lam = cpl_table_get_data_double(groupspec, "lam");

    /* CPL array for kernel */
    kernel = cpl_array_new(0, CPL_TYPE_DOUBLE);

    /* Convolve spectral ranges with kernel */

    do {

        /* Get lower and upper wavelengths of range */
        llam = lam[range[0]];
        ulam = llam * (1 + speedpar);

        /* Find upper range pixel */
        i = range[0];
        while (i < m && lam[i] < ulam) {
            i++;
        }
        if (i == m) {
            range[1] = m-1;
        } else {
            range[1] = i;
        }

        /* Get central wavelength of range */
        ulam = lam[range[1]];
        clam = (ulam + llam) / 2;

        /* Scale FWHM of kernel if variable kernel is requested */
        if (speedpar == HUGE_VAL) {
            cfwhm = fwhm;
        } else {
            cfwhm = fwhm * clam / reflam;
        }

        /* Calculate Gaussian kernel */
        sc_weights_calckernel(kernel, cfwhm);

        /* Convolve spectrum with Gaussian kernel */
        sc_basic_convolvewindow(convflux, flux, range, kernel);

        /* Set next lower range pixel */
        range[0] = range[1] + 1;

    } while (range[0] < m);

    /* Copy resulting flux array into "flux" column of input CPL table */
    vec = cpl_array_get_data_double(convflux);
    cpl_table_copy_data_double(groupspec, "flux", vec);

    /* Free memory */
    cpl_array_delete(flux);
    cpl_array_delete(convflux);
    cpl_array_delete(kernel);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_weights_calckernel(cpl_array *kernel, const double fwhm)
{
    /*!
     * Calculates Gaussian kernel depending on the FWHM provided by the
     * general parameter list. The number of kernel pixels is the upper odd
     * integer number of the input FWHM times ::SC_KERNFAC. The sum of the
     * kernel values is normalised to 1.
     *
     * \b INPUT:
     * \param fwhm    FWHM of Gaussian in pixels
     *
     * \b OUTPUT:
     * \param kernel  CPL array with convolution kernel elements
     *
     * \b ERRORS:
     * - Invalid object value(s)
     */

    char errtxt[SC_MAXLEN];
    int nkpix = 0, k = 0, nx = 0, refpix = 0, i = 0;
    double *kern;
    double sigma = 0, xmax = 0, xmin = 0, dx = 0, x = 0, sum = 0;

    /* Check FWHM */
    if (fwhm < 0.) {
        kernel = NULL;
        sprintf(errtxt, "%s: fwhm < 0", SC_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IIP, "%s", errtxt);
    }

    /* sigma of Gaussian */
    sigma = fwhm / CPL_MATH_FWHM_SIG;

    /* Number of kernel pixels */
    nkpix = 2 * ceil(fwhm * SC_KERNFAC / 2 - 0.5) + 1;
    cpl_array_set_size(kernel, nkpix);

    /* Kernel with one pixel only */
    if (nkpix == 1) {
        cpl_array_set(kernel, 0, 1.);
        return CPL_ERROR_NONE;
    }

    /* Get pointer to CPL array */
    cpl_array_fill_window_double(kernel, 0, nkpix, 0.);
    kern = cpl_array_get_data_double(kernel);

    /* Integration limits for upper bin */
    xmax = 0.5 * nkpix;
    xmin = xmax - 1;

    /* Number of wavelengths for integration of Gaussian */
    nx = ceil(SC_BINS_PER_FWHM / fwhm);

    /* Wavelength step for integration of Gaussian */
    dx = 1. / nx;

    /* Reference pixel for mirroring of kernel values */
    refpix = floor(nkpix / 2.);

    /* Calculate kernel up to reference pixel */

    for (k = nkpix - 1; k >= refpix; k--) {

        if (xmax <= 0.) {

            /* Skip integration */
            kern[k] = 1.;

            } else {

                /* First relative wavelength for integration */
                x = xmin + dx / 2;

                /* Perform integration */

                kern[k] = 0.;

                for (i = 0; i < nx; i++) {
                    kern[k] += exp(-0.5 * pow(x / sigma, 2));
                    x += dx;
                }

                kern[k] /= (double) nx;

            }

            /* Shift integration limits for next bin */
            xmax = xmin;
            xmin = xmax - 1;

        }

        /* Mirror right wing of kernel */
        for (k = refpix - 1; k >= 0; k--) {
            kern[k] = kern[nkpix - k - 1];
        }


    /* Add all kernel values */
    for (k = 0; k < nkpix; k++) {
        if (kern[k] < 0.) {
            kern[k] = 0.;
        }
        sum += kern[k];
    }

    /* Normalise kernel values */
    for (k = 0; k < nkpix; k++) {
        kern[k] /= sum;
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_weights_rebinspec(cpl_table *skyspec,
                                    const cpl_table *groupspec)
{
    /*!
     * Rebins line group spectrum to wavelength grid of input sky spectrum.
     * The resulting fluxes are written in two temporary columns in the sky
     * spectrum: "f_group" for the selected group and "f_all" for the sum
     * flux of all groups. This routine reverses the operations performed in
     * ::sc_weights_createwavegrid.
     *
     * \b INPUT:
     * \param skyspec    CPL table with sky spectrum
     * \param groupspec  CPL table with line group spectrum
     *
     * \b OUTPUT:
     * \param skyspec   CPL table with sky and rebinned line group spectrum
     *
     * \b ERRORS:
     * - Invalid object structure
     */

    char errtxt[SC_MAXLEN];
    int nrow = 0, npix = 0, jmin = 0, jmax = 0, jrow = 0, j = 0, k = 0, i = 0;
    double flux = 0.;
    const double * g_flux = NULL;

    /* Get number of pixels in input and output spectrum */
    nrow = cpl_table_get_nrow(skyspec);
    npix = cpl_table_get_nrow(groupspec);

    /* Get limiting pixels and check number of pixels */
    jmin = SC_EXTRACOVER * SC_SAMPFAC;
    jmax = npix - SC_EXTRACOVER * SC_SAMPFAC - 2;
    jrow = (jmax - jmin + 1) / SC_SAMPFAC;
    if (jrow != nrow) {
        sprintf(errtxt, "%s: cpl_table *groupspec (unexpected pixel number)",
                SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Sum up fluxes in line group spectrum to get pixel-specific fluxes in
       sky spectrum */
    if (cpl_table_count_invalid(groupspec, "flux") == 0)
        g_flux = cpl_table_get_data_double_const(groupspec, "flux");
    for (j = jmin; j <= jmax; j++) {
        if (k < SC_SAMPFAC) {
             flux += g_flux ?
                 g_flux[j] : cpl_table_get(groupspec, "flux", j, NULL);
             k++;
        }
        if (k == SC_SAMPFAC) {
             cpl_table_set(skyspec, "f_group", i, flux);
             flux = 0.;
             k = 0;
             i++;
        }
    }

    /* Add line group fluxes to total airglow flux */
    cpl_table_add_columns(skyspec, "f_all", "f_group");

    return CPL_ERROR_NONE;
}


cpl_error_code sc_weights_fillgrouparrays(cpl_table *skyspec,
                                          const char grouptype,
                                          const int group)
{
    /*!
     * Writes ID and flux contribution of the provided group into
     * pixel-specific arrays of the sky spectrum columns belonging to the
     * given group type. The array size depends on the pixel with the highest
     * number of contributing line groups. Unused array elements remain
     * invalid.
     *
     * \b INPUT:
     * \param skyspec    CPL table with sky and line group spectrum
     * \param grouptype  line group type (either 'A' or 'B')
     * \param group      number of line group
     *
     * \b OUTPUT:
     * \param skyspec    sky spectrum with updated pixel-specific group flux
     *                   arrays
     *
     * \b ERRORS:
     * - none
     */

    cpl_array *narray = NULL, *warray = NULL;
    char ncolname[SC_LENLINE+1], wcolname[SC_LENLINE+1];
    int nrow = 0, i = 0, size = 0, nfree = 0, ngroup = 0;
    double flux = 0.;

    /* Build table column names */
    sprintf(ncolname, "ng%c", grouptype);
    sprintf(wcolname, "wg%c", grouptype);

    /* Get number of table rows */
    nrow = cpl_table_get_nrow(skyspec);

    /* Fill the array columns with group IDs and line fluxes */

    for (i = 0; i < nrow; i++) {

        flux = cpl_table_get(skyspec, "f_group", i, NULL);

        if (flux > 0) {

            /* Get arrays */
            narray = cpl_array_duplicate(cpl_table_get_array(skyspec,
                                                             ncolname, i));
            warray = cpl_array_duplicate(cpl_table_get_array(skyspec,
                                                             wcolname, i));

            /* Number of unfilled array elements */
            size = cpl_array_get_size(narray);
            if (size == 0) {
                nfree = 0;
            } else {
                nfree = cpl_array_count_invalid(narray);
            }

            /* Resize array columns if more groups contribute to a pixel than
               allowed by the set depth of the array table columns in the sky
               spectrum */
            if (nfree == 0) {
                cpl_table_set_column_depth(skyspec, ncolname, size+1);
                cpl_table_set_column_depth(skyspec, wcolname, size+1);
                cpl_array_set_size(narray, size+1);
                cpl_array_set_size(warray, size+1);
            }

            /* New number of groups in the array */
            ngroup = size - nfree + 1;

            /* Add ID and flux of given group */
            cpl_array_set(narray, ngroup-1, group);
            cpl_array_set(warray, ngroup-1, flux);

            /* Fill table columns with modified arrays */
            cpl_table_set_array(skyspec, ncolname, i, narray);
            cpl_table_set_array(skyspec, wcolname, i, warray);

            /* Delete arrays */
            cpl_array_delete(narray);
            cpl_array_delete(warray);

        }

    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_weights_normpixcontrib(cpl_table *skyspec,
                                         const char grouptype)
{
    /*!
     * Divides group flux by flux sum for each pixel to get weights. The
     * operation is performed for the sky spectrum group-weight column
     * belonging to the given group type.
     *
     * \b INPUT:
     * \param skyspec    CPL table of sky spectrum with pixel-specific group
     *                   flux arrays
     * \param grouptype  line group type (either 'A' or 'B')
     *
     * \b OUTPUT:
     * \param skyspec    sky spectrum with pixel-specific group weight arrays
     *
     * \b ERRORS:
     * - none
     */

    cpl_array *warray = NULL;
    char wcolname[SC_LENLINE+1];
    int nrow = 0, i = 0;
    double sumflux = 0.;

    /* Build weight column name */
    sprintf(wcolname, "wg%c", grouptype);

    /* Get number of table rows */
    nrow = cpl_table_get_nrow(skyspec);

    /* Normalise fluxes in arrays by total flux */

    for (i = 0; i < nrow; i++) {

        sumflux = cpl_table_get(skyspec, "f_all", i, NULL);

        if (sumflux > 0) {

            /* Get flux array for pixel i */
            warray = cpl_array_duplicate(cpl_table_get_array(skyspec,
                                                             wcolname, i));
            /* Divide fluxes by total flux */
            cpl_array_divide_scalar(warray, sumflux);

            //printf("type %c: pixel %d:\n", grouptype, i+1);
            //cpl_array_dump(warray, 0, cpl_array_get_size(warray), NULL);

            /* Fill weight column with corrected array */
            cpl_table_set_array(skyspec, wcolname, i, warray);

            /* Delete array */
            cpl_array_delete(warray);

        }

    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_weights_getdomgroups(cpl_table *skyspec,
                                       const cpl_parameterlist *parlist,
                                       const char grouptype)
{
    /*!
     * Gets the IDs of the dominating line group of each pixel and writes them
     * into the sky spectrum table. A line group is only considered if it has
     * a minimum weight of \e weightlim which is provided by the general
     * parameter list. If such a group does not exist for a pixel, -99 is
     * written.
     *
     * \b INPUT:
     * \param skyspec    CPL table of sky spectrum with pixel-specific group
     *                   weight arrays
     * \param parlist    input CPL parameter list
     * \param grouptype  line group type (either 'A' or 'B')
     *
     * \b OUTPUT:
     * \param skyspec    sky spectrum table with column of dominating groups
     *
     * \b ERRORS:
     * - none
     */

    const cpl_parameter *p;
    cpl_array *narray = NULL, *warray = NULL;
    char ncolname[SC_LENLINE+1], wcolname[SC_LENLINE+1];
    char dcolname[SC_LENLINE+1];
    int nrow = 0, initid = -99, depth = 0, i = 0, n = 0, j = 0, group = 0;
    int mgroup = 0;
    double weightlim = 0., mwgroup = 0., wgroup = 0.;

    /* Build column names */
    sprintf(ncolname, "ng%c", grouptype);
    sprintf(wcolname, "wg%c", grouptype);
    sprintf(dcolname, "dg%c", grouptype);

    /* Get number of table rows */
    nrow = cpl_table_get_nrow(skyspec);

    /* Initialise column for dominating groups with initid */
    cpl_table_fill_column_window(skyspec, dcolname, 0, nrow, initid);

    /* Get depth of line group column */
    depth = cpl_table_get_column_depth(skyspec, ncolname);

    /* No groups -> return */
    if (depth == 0) {
        return CPL_ERROR_NONE;
    }

    /* Get weight limit from parameter list */
    p = cpl_parameterlist_find_const(parlist, "weightlim");
    weightlim = cpl_parameter_get_double(p);

    /* Get numbers of pixel-dominating groups */

    for (i = 0; i < nrow; i++) {

        /* Get arrays */
        narray = cpl_array_duplicate(cpl_table_get_array(skyspec, ncolname,
                                                         i));
        warray = cpl_array_duplicate(cpl_table_get_array(skyspec, wcolname,
                                                         i));

        /* Get number of contributing line groups */
        n = depth - cpl_array_count_invalid(narray);

        /* Skip pixel if no line group contributes */
        if (n == 0) {
            cpl_array_delete(narray);
            cpl_array_delete(warray);
            continue;
        }

        /* Find group with highest weight for given pixel */
        for (mwgroup = 0., j = 0; j < n; j++) {
            group = cpl_array_get(narray, j, NULL);
            wgroup = cpl_array_get(warray, j, NULL);
            if (wgroup > mwgroup) {
                mgroup = group;
                mwgroup = wgroup;
            }
        }

        /* Write number of strongest group in dgA or dgB column
           if weight >= weightlim from general parameter list */
        if (mwgroup >= weightlim) {
            cpl_table_set(skyspec, dcolname, i, mgroup);
        }

        /* Delete arrays */
        cpl_array_delete(narray);
        cpl_array_delete(warray);

    }

    return CPL_ERROR_NONE;
}

/**@}*/
