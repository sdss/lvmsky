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
 * \file sc_contsub.c
 *
 * Routines related to continuum identification, interpolation, and
 * subtraction
 *
 * \author Marco Barden, Stefan Noll, & ESO In-Kind Team Innsbruck
 *
 * \date   19 Sep 2013
 *
 */


/*****************************************************************************
 *                                  INCLUDES                                 *
 ****************************************************************************/

#include <sc_contsub.h>


/*****************************************************************************
 *                                    CODE                                   *
 ****************************************************************************/

cpl_error_code sc_contsub(cpl_table *spec)
{
    /*!
     * \brief
     *   Interpolate and subtract continuum from a spectrum.
     *
     * This routine uses the mask column in \em spec to identify emission
     * lines and interpolate the flux at their position using the closest
     * continuum points.
     *
     * The input spectrum has to be a CPL table containing (at least) the
     * following columns:
     *
     *      col #1: lambda
     *      col #2: flux
     *      col #3: mask
     *
     * The column "class" is is expected to contain the following values:
     *      0 = Continuum
     *      1 = Line pixel
     *      2 = Line peak
     *      3 = Isolated line peak
     *
     * On output the spectrum contains two new columns:
     *      col #4: cflux
     *      col #5: lflux
     * containing the continuum and line fluxes. Col #5, lflux, is the
     * continuum subtracted line flux. The sum of lflux and cflux corresponds
     * to col #2, flux.
     *
     * \b INPUT:
     * \param spec  spectrum
     *
     * \b OUTPUT:
     * \param spec  spectrum
     *
     * \b ERRORS:
     * - CPL_ERROR_NONE:           no error occurred
     * - CPL_ERROR_ILLEGAL_INPUT:  error in input spectra
     */

    cpl_error_code err_code = CPL_ERROR_NONE;

    cpl_table *tab;

    double *xr, *yr, *xo, *yo;              // pointers for interpolation

    long nrows = cpl_table_get_nrow(spec),  // number of rows in spec
         ncont = 0,                         // number of continuum data points
         i;

    int npix = 11;                          // pixels for continuum median

    // add new columns to spec
    if (cpl_table_has_column(spec, "cflux") == 0) {
        cpl_table_new_column(spec, "cflux", CPL_TYPE_DOUBLE);
    }
    if (cpl_table_has_column(spec, "lflux") == 0) {
        cpl_table_new_column(spec, "lflux", CPL_TYPE_DOUBLE);
    }

    // allocate space for temporary variables
    xr = (double *)calloc(nrows, sizeof(double));
    yr = (double *)calloc(nrows, sizeof(double));
    xo = (double *)calloc(nrows, sizeof(double));
    yo = (double *)calloc(nrows, sizeof(double));

    // count number of continuum and line data points
    for (i = 0; i < nrows; i++) {
        *(xo+i) = cpl_table_get_double(spec, "lambda", i, NULL);
        *(yo+i) = 0;
        if (cpl_table_get_int(spec, "class", i, NULL) == 0) {
            *(xr+ncont) = cpl_table_get_double(spec, "lambda", i, NULL);
            *(yr+ncont) = cpl_table_get_double(spec, "flux", i, NULL);
            ncont++;
        }
    }

    // smooth continuum points by median filtering
    tab = cpl_table_new(ncont);
    cpl_table_new_column(tab, "cflux", CPL_TYPE_DOUBLE);
    for (i = 0; i < ncont; i++) {
        cpl_table_set_double(tab, "cflux", i, *(yr+i));
    }
    sc_basic_filtermedian(tab, "cflux", npix);
    for (i = 0; i < ncont; i++) {
        *(yr+i) = cpl_table_get_double(tab, "cflux", i, NULL);
    }
    cpl_table_delete(tab);

    // interpolate over line regions
    err_code = sc_basic_interpollin(xo, yo, nrows, xr, yr, ncont, 0);

    // write results back to spectrum
    for (i = 0; i < nrows; i++) {
        if (cpl_table_get_int(spec, "class", i, NULL) == 0) {
            // keep original flux at continuum region
            cpl_table_set_double(spec, "cflux", i,
                             cpl_table_get_double(spec, "flux", i, NULL));
        } else {
            // set interpolated flux at line region
            cpl_table_set_double(spec, "cflux", i, *(yo + i));
        }
        // fill line flux region
        cpl_table_set_double(spec, "lflux", i,
                             cpl_table_get_double(spec, "flux", i, NULL));
    }

    // subtract continuum flux from line region
    cpl_table_subtract_columns(spec, "lflux", "cflux");

    // cleanup
    free(xr);
    free(yr);
    free(xo);
    free(yo);

    return err_code;
}


cpl_error_code sc_contsub_identcont(cpl_table *spec, const cpl_table *groups,
                                    const cpl_table *linetab,
                                    const double fluxlim,
                                    const cpl_parameterlist *parlist)
{
    /*!
     * Identifies continuum windows in a spectrum by excluding airglow lines
     * from a list with a peak flux above a threshold which is the product of
     * the input parameter 'fluxlim' and the median peak flux of the lines
     * already identified by a line finder. The remaining continuum pixels are
     * identified by a value of 0 in the 'mask' column of the spectrum table.
     * If there are no identified lines, the mean of the line list fluxes is
     * used instead of the median peak flux of the identified lines.
     *
     * \b INPUT:
     * \param spec     CPL table with spectrum
     * \param groups   CPL table with airglow line information
     * \param linetab  CPL table with information on detected lines
     * \param fluxlim  relative lower line peak flux limit
     * \param parlist  general CPL parameter list
     *
     * \b OUTPUT:
     * \param spec     spectrum with updated line/continuum flag
     *
     * \b ERRORS:
     * - see subroutines
     */

    cpl_table *selgroups, *subtab = NULL;
    const cpl_parameter *p;

    int i, j,                     // loop variables
        mid,                      // central line wavelength
        n = 0,                    // # of rows in table spec
        m = 0,                    // # of rows in table groups
        varfwhm = 0;              // wavelength-dependent FWHM?

    double fwhm = 0,              // FWHM of lines
           ratio = 0,             // flux ratio
           minflux = 0,           // minimum flux for lines from list
           nsig = 3.,    // line/cont separation extension based on line width
           lam_spec0, lam_spec1,  // neighbouring wavelengths in spec
           lam_list=0.,           // wavelength in line list
           dlam = 0,              // wavelength step
           meanlam = 0,           // mean wavelength
           flux = 0,              // flux
           lpix = 0,              // # of line pixels
           meanflux = 0.,         // mean list line flux
           meanlflux = 0.,        // mean observed line flux
           lflux = 0.;            // observed line flux

    /* Copy line group table */
    selgroups = cpl_table_duplicate(groups);

    /* Number of elements in tables */
    n = cpl_table_get_nrow(spec);
    m = cpl_table_get_nrow(selgroups);

    /* Clip groups to proper wavelength interval */
    lam_spec0 = cpl_table_get_column_min(spec, "lambda");
    lam_spec1 = cpl_table_get_column_max(spec, "lambda");
    cpl_table_unselect_all(selgroups);
    for (i = 0; i < m; i++) {
        lam_list = cpl_table_get_double(selgroups, "lambda", i, NULL);
        if (lam_list < lam_spec0 || lam_list > lam_spec1) {
            cpl_table_select_row(selgroups, i);
        }
    }
    cpl_table_erase_selected(selgroups);

    /* Update number of elements in selgroups */
    m = cpl_table_get_nrow(selgroups);

    /* delta(lambda): required for converting FWHM from [pixel] to [micron] */
    dlam = cpl_table_get_double(spec, "lambda", 1, NULL) -
           cpl_table_get_double(spec, "lambda", 0, NULL);

    /* Get FWHM */
    p = cpl_parameterlist_find_const(parlist, "fwhm");
    fwhm = cpl_parameter_get_double(p);

    /* Variable FWHM (linear change with wavelength)? */
    p = cpl_parameterlist_find_const(parlist, "varfwhm");
    varfwhm = cpl_parameter_get_int(p);

    /* Get mean wavelength */
    p = cpl_parameterlist_find_const(parlist, "meanlam");
    meanlam = cpl_parameter_get_double(p);

    /* Convert flux in group list from integrated to peak flux */
    for (i = 0; i < m; i++) {
        flux = cpl_table_get_double(selgroups, "flux", i, NULL) /
               CPL_MATH_SQRT2PI / CPL_MATH_SIG_FWHM / fwhm / dlam;
        /* Correct flux if FWHM is wavelength dependent */
        if (varfwhm == 1) {
            lam_list = cpl_table_get_double(selgroups, "lambda", i, NULL);
            flux *= meanlam / lam_list;
        }
        cpl_table_set_double(selgroups, "flux", i, flux);
    }

    /* Rebin selgroups spectrum */
    subtab = cpl_table_duplicate(spec);
    cpl_table_new_column(subtab, "flux_list", CPL_TYPE_DOUBLE);
    sc_basic_rebin(subtab, "lambda", "flux_list", selgroups, "lambda",
                   "flux");

    /* Apply scaling */
    meanflux = cpl_table_get_column_mean(subtab, "flux_list");
    ratio = cpl_table_get_column_mean(subtab, "lflux") / meanflux;
    if (ratio < 0) {
        /* Avoid negative ratios */
        for (meanlflux = 0., i = 0; i < n; i++) {
            if ((lflux = cpl_table_get(subtab, "lflux", i, NULL)) > 0) {
                meanlflux += lflux;
            }
        }
        ratio = meanlflux / meanflux;
    }
    cpl_table_multiply_scalar(selgroups, "flux", ratio);

    /* Get minimum line flux limit in spectrum */
    if (cpl_table_get_nrow(linetab) == 0) {
        /* Solution for no identified lines */
        minflux = fluxlim * cpl_table_get_column_mean(selgroups, "flux");
    } else {
        /* Fraction of median peak flux of identified lines */
        minflux = fluxlim * cpl_table_get_column_median(linetab, "peak_flux");
        /* Make sure that lines from the list are above the threshold */
        if (minflux > cpl_table_get_column_max(selgroups, "flux")) {
            minflux = fluxlim * cpl_table_get_column_mean(selgroups, "flux");
        }
    }

    /* Add line centres from line list to spectrum */

    for (i = 0, j = 0; i < n-1; ) {

        /* Check end of line list */
        if (j >= m-1) {
            break;
        }

        /* Get two consecutive wavelengths from spectrum */
        lam_spec0 = cpl_table_get_double(spec, "lambda", i, NULL);
        lam_spec1 = cpl_table_get_double(spec, "lambda", i+1, NULL);

        /* Get wavelength and flux in line list */
        lam_list = cpl_table_get_double(selgroups, "lambda", j, NULL);
        flux = cpl_table_get_double(selgroups, "flux", j, NULL);

        if (lam_list < lam_spec0) {
            /* Skip lines at low end of wavelength range in spectrum */
            j++;
            continue;
        }

        if (flux < minflux) {
            /* Use only significant lines */
            j++;
            continue;
        }

        if (lam_list < lam_spec1) {
            /* Current line in list has wavelength >= lam_spec0
               but < lam_spec1 -> find closer wavelength */
            if (lam_list-lam_spec0 < lam_spec1-lam_list) {
                /* Closer to first wavelength */
                mid = i;
            } else {
                /* Closer to second wavelength */
                mid = i+1;
            }
            /* Identify peak in mask (no change of flag for isolated lines) */
            if (cpl_table_get_int(spec, "class", mid, NULL) != 3) {
                cpl_table_set_int(spec, "class", mid, 2);
            }

            /* Next line in list */
            j++;
        }

        if (lam_list > lam_spec0) {
            /* Move to next wavelength in spectrum */
            i++;
        }

    }

    /* Loop over all pixels */

    for (i = 0; i < n-1; i++) {

        if (cpl_table_get_int(spec, "class", i, NULL) >= 2) {

            /* Peak found */

            /* Extend line regions in spectrum around peaks
               plus/minus nsig sigma [pixel] */
            lpix = nsig * fwhm / CPL_MATH_FWHM_SIG;

            /* Get peak flux and wavelength */
            flux = cpl_table_get_double(spec, "flux", i, NULL);
            lam_list = cpl_table_get_double(spec, "lambda", i, NULL);

            /* Correct peak flux if FWHM is wavelength dependent */
            if (varfwhm == 1) {
                flux *= lam_list / meanlam;
            }

            /* Strong lines are broader at the base */
            if (flux > 10 * fluxlim) {
                lpix *= 2;
            }

            /* Correct line regions if FWHM is wavelength dependent */
            if (varfwhm == 1) {
                lpix *= lam_list / meanlam;
            }

            /* Surround peak with line pixels */

            for (j = 1; j < (int) floor(lpix + 0.5); j++) {
                /* Lower wavelength end */
                if (i-j >= 0) {
                    if (cpl_table_get_int(spec, "class", i-j, NULL) < 1) {
                        /* Change continuum pixels only */
                        cpl_table_set_int(spec, "class", i-j, 1);
                    }
                }
                /* Higher wavelength end */
                if (i+j < n) {
                    if (cpl_table_get_int(spec, "class", i+j, NULL) < 1) {
                        /* Change continuum pixels only */
                        cpl_table_set_int(spec, "class", i+j, 1);
                    }
                }
            }

        }

    }

    /* Cleanup */
    cpl_table_delete(subtab);
    cpl_table_delete(selgroups);

    return CPL_ERROR_NONE;
}


int sc_contsub_check(cpl_table *spec)
{
    /*!
     * Checks distribution of continuum pixels over spectral range. If the
     * distribution allows the continuum to be interpolated in a reliable way,
     * a value of 0 is returned. The fraction of continuum pixels and the
     * fraction of the wavelength range covered by continuum pixels is
     * checked.
     *
     * \b INPUT:
     * \param spec  CPL table with spectrum
     *
     * \b OUTPUT:
     * - none
     *
     * \b ERRORS:
     * - none
     */

    cpl_table *subtab;
    int nall = 0, ncont = 0, flag = 0;
    double minlam = 0., maxlam = 0., range = 0., crange = 0.;

    /* Parameters for checking */
    int count = 5; // continuum pixel to be evaluated (at both margins)
    double minfpix = 0.2; // minimum fraction of continuum pixels
    double minfrange = 0.9; // minimum fraction of continuum wavelength range

    /* Get total number of pixels */
    nall = cpl_table_get_nrow(spec);
    if (nall == 0) {
        return 1;
    }

    /* Extract subtable of continuum pixels */
    cpl_table_unselect_all(spec);
    cpl_table_or_selected_int(spec, "class", CPL_EQUAL_TO, 0);
    subtab = cpl_table_extract_selected(spec);
    cpl_table_select_all(spec);

    /* Check fraction of continuum pixels */
    ncont = cpl_table_get_nrow(subtab);
    if ((double) ncont / (double) nall < minfpix) {
        flag = 1;
    }

    /* Get full wavelength range covered by the spectrum */
    minlam = cpl_table_get(spec, "lambda", 0, NULL);
    maxlam = cpl_table_get(spec, "lambda", nall-1, NULL);
    range = maxlam - minlam;

    /* Check wavelength range covered by continuum pixels */
    if (ncont >= 2 * count) {
        minlam = cpl_table_get(subtab, "lambda", count-1, NULL);
        maxlam = cpl_table_get(subtab, "lambda", ncont-count, NULL);
        crange = maxlam - minlam;
        if (crange / range < minfrange) {
            flag = 1;
        }
    } else {
        flag = 1;
    }

    /* Cleanup */
    cpl_table_delete(subtab);

    return flag;
}


cpl_error_code sc_contsub_skipthermir(cpl_table *spec)
{
    /*!
     * Identifies all pixels at wavelengths longer than ::SC_THERMIRLIM as
     * continuum pixels. This measure avoids the scaling of emission lines
     * originating in the lower atmosphere.
     *
     * \b INPUT:
     * \param spec  CPL table with spectrum
     *
     * \b OUTPUT:
     * \param spec  spectrum with continuum beyond thermal IR limit
     *
     * \b ERRORS:
     * - none
     */

    int n = 0, i = 0, mask = 0;
    double lam = 0.;

    /* Get number of pixles in spectrum */
    n = cpl_table_get_nrow(spec);

    /* Beyond SC_THERMIRLIM set all mask flags to 0 */
    for (i = n - 1; i >= 0; i--) {
        lam = cpl_table_get(spec, "lambda", i, NULL);
        mask = cpl_table_get(spec, "class", i, NULL);
        if (lam > SC_THERMIRLIM && mask > 0) {
            cpl_table_set(spec, "class", i, 0);
        }
    }

    return CPL_ERROR_NONE;
}

/**@}*/
