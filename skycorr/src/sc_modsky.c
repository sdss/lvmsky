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
 * \file sc_modsky.c
 *
 * Routines related to the modification of the sky spectrum
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  20 Feb 2011
 * \date   12 Jun 2014
 *
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <sc_mpfit.h>
#include <sc_modsky.h>

#include <math.h>
#include <assert.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

cpl_error_code sc_modsky(cpl_table *scispec, cpl_table *skyspec,
                         cpl_table *fitpar, const cpl_vector *sinc,
                         const cpl_parameterlist *parlist)
{
    /*!
     * \callgraph
     *
     * Modifies sky spectrum by scaling the different line groups and by
     * changing the wavelength grid by means of a Chebyshev polynomial. The
     * resulting spectrum is rebinned to match the wavelength grid of the
     * science spectrum. Finally, weighted deviations between modified sky
     * and science spectrum are computed, which is required for the
     * \f${\chi^2}\f$ derivation. Spectral lines with unreasonable ratios
     * (probably due to an object emission line) are excluded from this
     * computation.
     *
     * \b INPUT:
     * \param scispec  CPL table with science spectrum
     * \param skyspec  CPL table with sky spectrum and line group weights
     * \param fitpar   CPl table of fit parameters
     * \param sinc     CPL vector with damped sinc kernel
     * \param parlist  general CPL parameter list
     *
     * \b OUTPUT:
     * \param scispec  table with science spectrum and adapted sky spectrum
     * \param skyspec  table with original and modified sky spectrum
     *
     * \b ERRORS:
     * - none
     */

    int nrow = 0;

    /* Prepare columns in sky spectrum table */
    nrow = cpl_table_get_nrow(skyspec);
    cpl_table_fill_column_window(skyspec, "mlflux", 0, nrow, 0.);
    cpl_table_add_columns(skyspec, "mlflux", "lflux");
    if (lastcall == CPL_TRUE &&
        cpl_table_has_column(skyspec, "mdflux") == 1) {
        cpl_table_fill_column_window(skyspec, "mdflux", 0, nrow, 0.);
        cpl_table_add_columns(skyspec, "mdflux", "dflux");
    }
    cpl_table_fill_column_window(skyspec, "mlflux", 0, nrow, 0.);
    cpl_table_add_columns(skyspec, "mlflux", "lflux");
    cpl_table_fill_column_window(skyspec, "mweight", 0, nrow, 0.);
    cpl_table_add_columns(skyspec, "mweight", "weight");

    /* Modify line group fluxes in sky spectrum */
    sc_modsky_modlines(skyspec, fitpar, 'A'); /* A groups */
    sc_modsky_modlines(skyspec, fitpar, 'B'); /* B groups */

    /* Modify flux errors if provided by sky spectrum table */
    sc_modsky_moderrors(skyspec);

    /* Modify wavelength grid of sky spectrum by means of a Chebyshev
       polynomial */
    sc_modsky_modwavegrid(skyspec, fitpar);

    /* Rebin modified sky spectrum to wavelength grid of science spectrum */
    sc_modsky_rebin(scispec, skyspec, sinc, parlist);

    /* Calculate weighted deviations between modified sky and science
       spectrum */
    sc_modsky_calcdev(scispec);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_modsky_modlines(cpl_table *skyspec, cpl_table *fitpar,
                                  const char grouptype)
{
    /*!
     * Modifies sky spectrum by scaling the different line groups of a given
     * type. The group-specific factors are taken from the fit parameter table
     * fitpar. The pixel-specific weights of the different line groups have
     * been added to the sky spectrum by ::sc_weights.
     *
     * \b INPUT:
     * \param skyspec    CPL table with sky spectrum and line group weights
     * \param fitpar     CPl table of fit parameters
     * \param grouptype  line group type (either 'A' or 'B')
     *
     * \b OUTPUT:
     * \param skyspec    sky spectrum table with modified line fluxes
     *
     * \b ERRORS:
     * - none
     */

    cpl_table *grouppar;
    char type[2], ncolname[SC_LENLINE+1], wcolname[SC_LENLINE+1];
    int *rel = NULL;
    int ngroup = 0, depth = 0, nrow = 0, i = 0, n = 0, j = 0, group = 0;
    double *val = NULL;
    const cpl_array ** pa_narray, ** pa_warray;
    double * pd_flux;

    /* Create subtable from fit parameter table that only contains the
       coefficients for the wavelength correction */
    sprintf(type, "%c", grouptype);
    cpl_table_unselect_all(fitpar);
    cpl_table_or_selected_string(fitpar, "type", CPL_EQUAL_TO, type);
    grouppar = cpl_table_extract_selected(fitpar);
    cpl_table_select_all(fitpar);

    /* Get number of line groups */
    ngroup = cpl_table_get_nrow(grouppar);

    /* Build group number and weight column name */
    sprintf(ncolname, "ng%c", grouptype);
    sprintf(wcolname, "wg%c", grouptype);

    /* Get depth of line group column */
    depth = cpl_table_get_column_depth(skyspec, ncolname);

    /* No data -> print warning message and return without change */
    if (ngroup == 0 || depth == 0) {
        if (nfev == 1) {
            /* Print message only once */
            cpl_msg_warning(cpl_func, "No %c line groups", grouptype);
        }
        cpl_table_delete(grouppar);
        return CPL_ERROR_NONE;
    }

    /* Get pointer to fit parameter values and relevance flags */
    val = cpl_table_get_data_double(grouppar, "value");
    rel = cpl_table_get_data_int(grouppar, "relevance");

    /* Get number of rows in sky spectrum table */
    nrow = cpl_table_get_nrow(skyspec);

    /* Modify flux of sky line spectrum depending on pixel-specific line
       group weights */

    assert(cpl_table_count_invalid(skyspec, wcolname) == 0);
    assert(cpl_table_count_invalid(skyspec, ncolname) == 0);
    pa_narray = cpl_table_get_data_array_const(skyspec, ncolname);
    pa_warray = cpl_table_get_data_array_const(skyspec, wcolname);
    pd_flux = cpl_table_get_data_double(skyspec, "mlflux");
    cpl_table_fill_invalid_double(skyspec, "mlflux", 0.);
    for (i = 0; i < nrow; i++) {
        double flux = pd_flux[i];

        /* Skip pixel if flux is zero */
        if (flux > 0) {
            double wsum, wgroup;

            /* Get arrays */
            const cpl_array * narray = pa_narray[i];
            const cpl_array * warray = pa_warray[i];

            /* Get number of contributing line groups */
            n = depth - cpl_array_count_invalid(narray);

            /* Skip pixel if no line group contributes */
            if (n == 0) {
                continue;
            }

            /* Get weight of each line group and multiply it by the
               corresponding fit parameter value */
            wsum = 0;
            for (j = 0; j < n; j++) {
                group = cpl_array_get(narray, j, NULL);
                wgroup = cpl_array_get(warray, j, NULL);
                if (group <= 0 || rel[group-1] == 0) {
                    /* Do not change flux for group with relevance flag = 0
                       and group 0 */
                  wsum += wgroup;
                } else {
                  wsum += wgroup * val[group-1];
                }
            }

            /* Multiply flux by weight sum to get modified pixel flux */
            flux *= wsum;

            /* must be valid or flux > 0 check would have been false */
            pd_flux[i] = flux;

        }

    }

    /* Free memory */
    cpl_table_delete(grouppar);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_modsky_moderrors(cpl_table *skyspec)
{
    /*!
     * Modifies flux errors if a corresponding column is provided by the
     * input sky spectrum table. The modification is related to the scaling
     * of the line group fluxes, which was calculated by ::sc_modsky_modlines.
     *
     * \b INPUT:
     * \param skyspec    CPL table with sky spectrum
     *
     * \b OUTPUT:
     * \param skyspec    sky spectrum table with modified flux errors
     *
     * \b ERRORS:
     * - none
     */

    int nrow = 0, i = 0;
    double *cflux = NULL, *lflux = NULL, *mlflux = NULL, *dflux = NULL;
    double *mdflux = NULL;

    /* Return without change if the modsky call is not the last one or an
       error column is lacking */
    if (lastcall == CPL_FALSE ||
        cpl_table_has_column(skyspec, "mdflux") != 1) {
        return CPL_ERROR_NONE;
    }

    /* Get pointers to table columns */
    cflux = cpl_table_get_data_double(skyspec, "cflux");
    lflux = cpl_table_get_data_double(skyspec, "lflux");
    mlflux = cpl_table_get_data_double(skyspec, "mlflux");
    dflux = cpl_table_get_data_double(skyspec, "dflux");
    mdflux = cpl_table_get_data_double(skyspec, "mdflux");

    /* Get number of rows in sky spectrum table */
    nrow = cpl_table_get_nrow(skyspec);

    /* Modify flux errors */
    for (i = 0; i < nrow; i++) {
        mdflux[i] = dflux[i];
        if (cflux[i] + lflux[i] > 0. && cflux[i] + mlflux[i] > 0.) {
            mdflux[i] *= (cflux[i] + mlflux[i]) / (cflux[i] + lflux[i]);
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_modsky_modwavegrid(cpl_table *skyspec, cpl_table *fitpar)
{
    /*!
     * Modifies sky spectrum by changing the wavelength grid by means of a
     * Chebyshev polynomial. The degree of this polynomial is provided by the
     * general parameter \e cheby_max and should be at least 1. If -1 is
     * given, the wavelength correction is skipped. The coefficients are
     * provided by the fit parameter table fitpar. Only coefficients with fit
     * flag = 1 are considered. The polynomial is applied to a temporary,
     * normalised wavelength grid that ranges from -1 to 1.
     *
     * \b INPUT:
     * \param skyspec  CPL sky spectrum table with modified line fluxes
     * \param fitpar   CPl table of fit parameters
     *
     * \b OUTPUT:
     * \param skyspec  sky spectrum table with modified line fluxes and
     *                 modified wavelength grid
     *
     * \b ERRORS:
     * - none
     */

    cpl_array *cheby;
    cpl_table *wpar;
    int nlam = 0, ncoef = 0, i = 0, maxpix = -1, j = 0;
    int *class;
    double limlam[2] = {0., 0.}, dlam = 0., wmean = 0., dellam = 0.;
    const double *par;
    double *slam, *mlam, *t;

    /* Default output wavelength grid = input wavelength grid */
    nlam = cpl_table_get_nrow(skyspec);
    cpl_table_fill_column_window(skyspec, "mlambda", 0, nlam, 0.);
    cpl_table_add_columns(skyspec, "mlambda", "lambda");

    /* Create subtable from fit parameter table that only contains the
       relevant coefficients for the wavelength correction */
    cpl_table_unselect_all(fitpar);
    cpl_table_or_selected_string(fitpar, "type", CPL_EQUAL_TO, "w");
    cpl_table_and_selected_int(fitpar, "relevance", CPL_EQUAL_TO, 1);
    wpar = cpl_table_extract_selected(fitpar);
    cpl_table_select_all(fitpar);

    /* Get number of coefficients for wavelength correction */
    ncoef = cpl_table_get_nrow(wpar);

    /* Fit flag = 0 for coefficient < maximum degree in selected fit parameter
       table -> print warning message and return without change */
    if (ncoef > 0 && ncoef - 1 < cpl_table_get(wpar, "N", ncoef - 1, NULL)) {
        if (nfev == 1) {
            /* Print message only once */
            cpl_msg_warning(cpl_func, "Invalid set of coefficients "
                            "for wavelength correction -> skip");
        }
        cpl_table_delete(wpar);
        return CPL_ERROR_NONE;
    }

    /* N of coef. < 2 -> print warning message if N = 1 and return without
       change */
    if (ncoef < 2) {
        if (nfev == 1 && ncoef == 1) {
            /* Print message only once */
            cpl_msg_warning(cpl_func, "Only one coefficient "
                            "for wavelength correction -> skip");
        }
        cpl_table_delete(wpar);
        return CPL_ERROR_NONE;
    }

    /* Get pointer to fit parameters */
    par = cpl_table_get_data_double(wpar, "value");

    /* Scale wavelengths in the way that the interval [-1,1] is covered */
    limlam[0] = cpl_table_get(skyspec, "lambda", 0, NULL);
    limlam[1] = cpl_table_get(skyspec, "lambda", nlam-1, NULL);
    dlam = limlam[1] - limlam[0];
    wmean = (limlam[0] + limlam[1]) / 2;
    cpl_table_duplicate_column(skyspec, "slambda", skyspec, "lambda");
    cpl_table_subtract_scalar(skyspec, "slambda", wmean);
    cpl_table_divide_scalar(skyspec, "slambda", dlam / 2.);

    /* Get pointer to normalised wavelength grid */
    slam = cpl_table_get_data_double(skyspec, "slambda");

    /* Set wavelengths in output column "mlambda" to 0 and get pointer */
    cpl_table_fill_column_window(skyspec, "mlambda", 0, nlam, 0.);
    mlam = cpl_table_get_data_double(skyspec, "mlambda");

    /* Find last line pixel for wavelength shift in the thermal IR */
    class = cpl_table_get_data_int(skyspec, "class");
    maxpix = nlam - 1;
    if (limlam[1] >= SC_THERMIRLIM) {
        for (i = nlam-1; i >= 0; i--) {
            if (class[i] > 0) {
                maxpix = i;
                break;
            }
        }
    }

    /* Create array for Chebyshev polynomials */
    cheby = cpl_array_new(ncoef, CPL_TYPE_DOUBLE);
    t = cpl_array_get_data_double(cheby);
    t[0] = 1;

    /* Compute Chebyshev polynomials */
    for (i = 0; i < nlam; i++) {
        if (i <= maxpix) {
            for (j = 0; j < ncoef; j++) {
                if (j == 1) {
                    t[j] = slam[i];
                } else if (j > 1) {
                    t[j] = 2 * slam[i] * t[j-1] - t[j-2];
                }
                mlam[i] += par[j] * t[j];
            }
            if (i == maxpix) {
                dellam = mlam[i] - slam[i];
            }
        } else {
            /* constant shift for thermal infrared */
            mlam[i] = slam[i] + dellam;
        }
    }

    /* Rescale wavelengths */
    cpl_table_multiply_scalar(skyspec, "mlambda", dlam / 2.);
    cpl_table_add_scalar(skyspec, "mlambda", wmean);

    /* Free memory */
    cpl_array_delete(cheby);
    cpl_table_delete(wpar);
    cpl_table_erase_column(skyspec, "slambda");

    return CPL_ERROR_NONE;
}


cpl_error_code sc_modsky_rebin(cpl_table *scispec, cpl_table *skyspec,
                               const cpl_vector *sinc,
                               const cpl_parameterlist *parlist)
{
    /*!
     * Rebins the modified sky spectrum (continuum, line flux, and possible
     * flux errors) to the wavelength grid of the science spectrum by means of
     * either a summation of fractional input pixel contributions to the
     * output pixels or a convolution with an asymmetric, exponentially-damped
     * sinc kernel (default), which avoids undesired line broadening. The
     * rebinning approach is provided by the integer parameter REBINTYPE of
     * the general parameter list. The results are written into the science
     * spectrum table. The weight is rebinned in a traditional way to avoid
     * negative values. Another difference is the handling of weight = 0,
     * which is substituted in the sky spectrum table by a minimum weight
     * depending on the factor ::SC_RELMAXERR. In the science spectrum table,
     * the resulting weight is then converted back to weight = 0 if it is
     * below a threshold. This approach allows one to consider rebinned pixels
     * with only a small contribution from an initial sky pixel with
     * weight = 0.
     *
     * \b INPUT:
     * \param scispec  CPL table with science spectrum
     * \param skyspec  CPL table with modified sky spectrum
     * \param sinc     CPL vector with damped sinc kernel
     * \param parlist  general CPL parameter list
     *
     * \b OUTPUT:
     * \param scispec  science spectrum table with rebinned modified sky
     *                 spectrum
     *
     * \b ERRORS:
     * - none
     */

    const cpl_parameter *p;
    int rebintype = 0, nsky = 0, i = 0, pmin = 0, pmax = 0, nsci = 0;
    double wmin = 0., maxerr = 0.;
    double * pd_mweight , * pd_tdflux, * pdsky_weight;

    /* Derive pixel shifts */
    sc_modsky_getpixelshifts(skyspec, scispec);

    /* Get method for rebinning from parameter list */
    p = cpl_parameterlist_find_const(parlist, "rebintype");
    rebintype = cpl_parameter_get_int(p);

    /* Rebin sky continuum spectrum by using selected type of rebinning
       (only in the case of last modsky call) */
    if (lastcall == CPL_TRUE) {
        if (rebintype == 1) {
            sc_modsky_sincrebin(scispec, "mcflux", skyspec, "cflux", sinc);
        } else {
            sc_basic_rebin(scispec, "lambda", "mcflux",
                           skyspec, "mlambda", "cflux");
        }
    }

    /* Rebin modified sky line spectrum by using selected type of rebinning */
    if (rebintype == 1) {
        sc_modsky_sincrebin(scispec, "mlflux", skyspec, "mlflux", sinc);
    } else {
        sc_basic_rebin(scispec, "lambda", "mlflux",
                       skyspec, "mlambda", "mlflux");
    }

    /* Rebin flux errors (if present) by using selected type of rebinning
       (only in the case of last modsky call) */
    if (lastcall == CPL_TRUE &&
        cpl_table_has_column(skyspec, "mdflux") == 1) {
        if (rebintype == 1) {
            sc_modsky_sincrebin(scispec, "mdflux", skyspec, "mdflux", sinc);
        } else {
            sc_basic_rebin(scispec, "lambda", "mdflux",
                           skyspec, "mlambda", "mdflux");
        }
    }

    /* Create temporary error columns in sky and science spectrum */
    cpl_table_new_column(scispec, "tdflux", CPL_TYPE_DOUBLE);
    cpl_table_new_column(skyspec, "tdflux", CPL_TYPE_DOUBLE);

    /* Get minimum non-zero weight of sky spectrum */
    nsky = cpl_table_get_nrow(skyspec);
    cpl_table_fill_invalid_double(skyspec, "weight", 0.);
    pdsky_weight = cpl_table_get_data_double(skyspec, "weight");
    for (i = 0; i < nsky; i++) {
        double weight = pdsky_weight[i];
        if (weight != 0 && weight < wmin) {
            wmin = weight;
        }
    }

    /* Derive maximum error for weight = 0 substitution */
    wmin /= (1 + SC_TOL);
    maxerr = SC_RELMAXERR / wmin;

    /* Fill error column of sky spectrum */
    pd_tdflux = cpl_table_get_data_double(skyspec, "tdflux");
    /* make sure column is valid */
    cpl_table_fill_column_window_double(skyspec, "tdflux", 0, nsky, 0.);
    for (i = 0; i < nsky; i++) {
        double weight = pdsky_weight[i];
        /* Handle weight = 0 */
        if (weight == 0) {
            pd_tdflux[i] = maxerr;
        } else {
            pd_tdflux[i] = 1. / weight;
        }
    }

    /* Rebin errors of sky spectrum
       (no sinc function to avoid negative errors) */
    sc_basic_rebin(scispec, "lambda", "tdflux", skyspec, "mlambda", "tdflux");

    /* Get minimum and maximum pixel in column "mpix" to identify
       non-overlapping wavelength ranges */
    pmin = cpl_table_get(skyspec, "mpix", 0, NULL);
    pmax = cpl_table_get(skyspec, "mpix", nsky-1, NULL);

    /* Convert errors back to weights */
    nsci = cpl_table_get_nrow(scispec);
    pd_mweight = cpl_table_get_data_double(scispec, "mweight");
    assert(cpl_table_count_invalid(scispec, "tdflux") == 0);
    pd_tdflux = cpl_table_get_data_double(scispec, "tdflux");
    /* make sure column is valid and fill in weight = 0 */
    cpl_table_fill_column_window_double(scispec, "mweight", 0, nsci, 0.);
    for (i = 0; i < nsci; i++) {
        double err = pd_tdflux[i];
        /* Identify weight = 0 */
        if (!(err > 1. / wmin || i < pmin || i > pmax)) {
            pd_mweight[i] = 1. / err;
        }
    }

    /* Delete temporary columns */
    cpl_table_erase_column(skyspec, "tdflux");
    cpl_table_erase_column(scispec, "tdflux");

    return CPL_ERROR_NONE;
}


cpl_error_code sc_modsky_getpixelshifts(cpl_table *skyspec,
                                        const cpl_table *scispec)
{
    /*!
     * Compares the wavelength grids of science and modified sky spectrum and
     * writes the deviations in the temporary sky spectrum columns "mpix" and
     * "dpix". For a pixel in the sky spectrum "mpix" gives the number of the
     * pixel in the science spectrum with the smallest wavelength difference.
     * "dpix" provides the subpixel shift of the sky spectrum relative to the
     * science spectrum.
     *
     * \b INPUT:
     * \param skyspec  CPL table with modified sky spectrum
     * \param scispec  CPL table with science spectrum
     *
     * \b OUTPUT:
     * \param skyspec  modified sky spectrum with pixel shifts relative to
     *                 science spectrum
     *
     * \b ERRORS:
     * - none
     */

    int nsci = 0, nsky = 0, *mpix, i = 0, j = 1;
    double *skylam, *dpix, dlam = 0., shift = 0.;
    const double *scilam;

    /* Get number of data points in science and sky spectrum */
    nsci = cpl_table_get_nrow(scispec);
    nsky = cpl_table_get_nrow(skyspec);

    /* Initialise "mpix" and "dpix" columns in sky spectrum */
    cpl_table_fill_column_window_int(skyspec, "mpix", 0, nsky, 0);
    cpl_table_fill_column_window_double(skyspec, "dpix", 0, nsky, 0.);

    /* Get pointers to CPL table columns */
    scilam = cpl_table_get_data_double_const(scispec, "lambda");
    skylam = cpl_table_get_data_double(skyspec, "mlambda");
    mpix = cpl_table_get_data_int(skyspec, "mpix");
    dpix = cpl_table_get_data_double(skyspec, "dpix");

    /* Find central pixel and subpixel shift for the kernel for each pixel in
       sky spectrum */

    for (i = 0; i < nsky; i++) {

        while (j < nsci - 1 && scilam[j] < skylam[i]) {
            j++;
        }

        dlam = scilam[j] - scilam[j-1];
        shift = (skylam[i] - scilam[j]) / dlam;
        mpix[i] = (int) floor(shift + 0.5) + j;
        dpix[i] = shift - (double) mpix[i] + j;

    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_modsky_sincrebin(cpl_table *scispec, const char *scicol,
                                   const cpl_table *skyspec,
                                   const char *skycol,
                                   const cpl_vector *sinc)
{
    /*!
     * Uses a convolution with a damped sinc kernel to rebin the data of the
     * column \e skycol in the sky spectrum to the wavelength grid of the
     * science spectrum. The differences in the wavelength grids have to be
     * provided by the columns "mpix" and "dpix" (see
     * ::sc_modsky_getpixelshifts). The output is written to the column
     * \e scicol in the science spectrum. Taking an exponentially-damped sinc
     * function for the rebinning avoids an increase of the FWHM of spectral
     * lines, which would significantly deteriorate any fit. For pixels of the
     * convolution kernel for the science spectrum outside of the wavelength
     * range of the sky spectrum, the flux of the first or last pixel of the
     * sky spectrum is taken.
     *
     * \b INPUT:
     * \param scispec  CPL table with science spectrum
     * \param scicol   column name for rebinned fluxes
     * \param skyspec  CPL table with modified sky spectrum and columns "mpix"
     *                 and "dpix"
     * \param skycol   column name for input fluxes
     * \param sinc     CPL vector with damped sinc kernel
     *
     *
     * \b OUTPUT:
     * \param scispec  science spectrum table with rebinned modified sky
     *                 spectrum
     *
     * \b ERRORS:
     * - none
     */

    int nsci = 0, nsky = 0, pmin = 0, pmax = 0, i = 0, k = 0, nsinc = 0;
    const int *mpix;
    double *cflux;
    const double *dpix, *flux, *psinc;

    /* Get kernel size */
    const double radius = 5.;
    const int nkpix = 2 * radius + 1;
    double kernel[nkpix];
    /* indexing offset into precomputed kernel */
    const double k_offset = ((double)(SC_SINCRAD_PRECOMP) - (radius));
    assert(k_offset >= 0);

    /* Get number of data points in science and sky spectrum */
    nsci = cpl_table_get_nrow(scispec);
    nsky = cpl_table_get_nrow(skyspec);

    /* Initialise output column with 0 values */
    cpl_table_fill_column_window_double(scispec, scicol, 0, nsci, 0.);

    /* Get pointer to output column */
    cflux = cpl_table_get_data_double(scispec, scicol);

    /* Get pointers to required table columns of sky spectrum */
    mpix = cpl_table_get_data_int_const(skyspec, "mpix");
    dpix = cpl_table_get_data_double_const(skyspec, "dpix");
    flux = cpl_table_get_data_double_const(skyspec, skycol);

    /* Get minimum and maximum pixel in column "mpix" */
    pmin = mpix[0];
    pmax = mpix[nsky-1];

    /* Get number of data points in sinc vector */
    nsinc = cpl_vector_get_size(sinc);

    /* for interpolation of kernel we need a precomputed kernel one larger than
     * the maximum shift of 0.5 */
    assert(((nkpix - 1.) + 0.5 + k_offset) * SC_SINCNBIN + 1 < nsinc);
    assert(((0 - 0.5 + k_offset)) * SC_SINCNBIN > 0);

    /* Get pointer to sinc vector */
    psinc = cpl_vector_get_data_const(sinc);

    /* Do convolution in science spectrum table for each pixel of the modified
       sky spectrum */

    for (i = 0; i < nsky; i++) {
        double shift, rsum, sum;
        int npix, sign, h;

        if (dpix[i] > 0.5 || dpix[i] < -0.5) {
            return cpl_error_set_message(cpl_func, SC_ERROR_IIP,
                                         "%s: shift < -0.5 or > 0.5",
                                         SC_ERROR_IIP_TXT);
        }

        shift = dpix[i] - k_offset;
        /* Get kernel values for given subpixel shift from sinc vector */
        for (sum = 0., k = 0; k < nkpix; k++) {
            double x = (double) (k - shift) * SC_SINCNBIN;
            double low = psinc[(intptr_t)x];
            double high = psinc[(intptr_t)(x + 1)];
            kernel[k] = low + (high - low) * (x - floor(x));
            sum += kernel[k];
        }

        /* Normalise kernel values */
        rsum = 1. / sum;
        for (k = 0; k < nkpix; k++) {
            kernel[k] *= rsum;
        }

        /* Number of central pixels for convolution with same kernel
           (> 1 only possible for margins) */
        if (i == 0 && pmin > -radius) {
            /* Lower margin */
            npix = pmin + radius + 1;
            sign = -1;
        } else if (i == nsky - 1 && pmax < nsci - 1 + radius) {
            /* Upper margin */
            npix = nsci - pmax + radius;
            sign = +1;
        } else {
            npix = 1;
            sign = +1;
        }

        /* Convolve pixel(s) with kernel */
        for (h = 0; h < npix; h++) {
            int j;
            for (k = 0, j = mpix[i] - radius + sign * h; k < nkpix;
                 k++, j++) {
                if (j < 0 || j >= nsci) {
                    /* Skip invalid pixels */
                    continue;
                }
                cflux[j] += flux[i] * kernel[k];
            }
        }

    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_modsky_calcdev(cpl_table *scispec)
{
    /*!
     * Computes weighted deviations between modified sky line spectrum and
     * science line spectrum, which is required for the \f${\chi^2}\f$
     * derivation. The weights used are combinations of the weights of the
     * science and the rebinned sky spectrum and the \f${\sigma}\f$-clipping
     * results from ::sc_mpfit_modinitpar. Only line peaks can have non-zero
     * weights.
     *
     * \b INPUT:
     * \param scispec  CPL table with science spectrum and rebinned modified
     *                 sky spectrum
     *
     * \b OUTPUT:
     * \param scispec  science spectrum table with weighted deviations between
     *                 modified sky and science spectrum
     *
     * \b ERRORS:
     * - none
     */

    size_t i, nrow;
    const double * pdweight, * pdmweight;
    const int * pisigclip;
    double * pdcweight;

    /* Calculate combined pixel weights of science spectrum and rebinned
       modified sky spectrum under consideration of the sigma clipping
       results for science data pixels from the initial estimate */

    nrow = cpl_table_get_nrow(scispec);
    cpl_table_fill_invalid_double(scispec, "weight", 0.);
    cpl_table_fill_invalid_double(scispec, "mweight", 0.);
    pdweight = cpl_table_get_data_double_const(scispec, "weight");
    pdmweight = cpl_table_get_data_double_const(scispec, "mweight");
    pisigclip = cpl_table_get_data_int_const(scispec, "sigclip");
    pdcweight = cpl_table_get_data_double(scispec, "cweight");
    cpl_table_fill_column_window(scispec, "cweight", 0, nrow, 0.);
    for (i = 0; i < nrow; i++) {
        double weff;
        double wsci = pdweight[i];
        double wsky = pdmweight[i];
        int sigclip = pisigclip[i];
        if (wsci > 0. && wsky > 0. && sigclip == 0) {
            weff = 1. / sqrt(1. / (wsci * wsci) + 1. / (wsky * wsky));
        } else {
            weff = 0.;
        }
        pdcweight[i] = weff;
    }

    /* Calculate weighted deviations between modified sky and science
       spectrum */
    cpl_table_fill_column_window(scispec, "dev", 0, nrow, 0.);
    cpl_table_add_columns(scispec, "dev", "mlflux");
    cpl_table_subtract_columns(scispec, "dev", "lflux");
    cpl_table_multiply_columns(scispec, "dev", "cweight");

    return CPL_ERROR_NONE;
}

/**@}*/
