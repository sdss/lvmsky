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
 * \callgraph
 *
 * \file sc_basic.c
 *
 * Basic routines for sky correction
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  15 Feb 2011
 * \date   24 Feb 2014
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <sc_basic.h>

#include <math.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

cpl_error_code sc_basic_readline(FILE *stream, scpar par[], int *npar)
{
    /*!
     * Reads a line from a file and returns the space- or tab-separated
     * individual strings in an ::scpar structure. This structure also
     * contains each string converted into integer and double precision
     * floating point numbers. Comment lines marked by # and empty lines are
     * not considered. In this case, the next line with data is read.
     *
     * \b INPUT:
     * \param stream  name of stream (to be defined by fopen command)
     * \param par     array of ::scpar structures
     *                (to transfer parameter values)
     *
     * \b OUTPUT:
     * \param par     ::scpar array that contains the read value(s) as
     *                character string ("c"), integer ("i"), or double
     *                precision floating point number ("d").
     *                The different data types can be accessed by adding the
     *                corresponding suffixes c, i, or d to the name of the
     *                ::scpar structure.
     * \param npar    found number of parameter values (at one line)
     *
     * \b ERRORS:
     * - Invalid input parameter(s)
     * - Unexpected file structure
     * - Insufficient memory
     *
     * In the case of problems, the default values "", 0, and 0.0 are used
     * for the corresponding data types.
     */

    char errtxt[SC_MAXLEN], line[SC_MAXLEN], **str = NULL;
    cpl_boolean fl_data = CPL_FALSE;
    int npar0 = 0, i = 0, nparx = 0;

    /* Set maximum number of parameters */
    npar0 = SC_MAXPAR;

    /* Valid number of parameters? */

    if (npar0 < 1) {
        sprintf(errtxt, "%s: npar < 1 (parameters at line)",
                SC_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IIP, "%s", errtxt);
    }

    /* Fill mfpar structure with 0 */

    for (i = 0; i < npar0; i++) {
        par[i].c[0] = '\0';
        par[i].i = 0;
        par[i].d = 0.;
    }

    /* Skip comments and empty lines and check existence of data */

    while (fgets(line, SC_MAXLEN, stream) != NULL) {
        if (line[0] != '#' && strspn(line, "\n\t ") != strlen(line)) {
            fl_data = CPL_TRUE;
            break;
        }
    }

    if (fl_data == CPL_FALSE) {
        sprintf(errtxt, "%s: FILE *stream (no parameters)",
                SC_ERROR_UFS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s", errtxt);
    }

    /* Allocate memory to string array */

    str = (char **) calloc(npar0+1, sizeof(char *));
    if (str == NULL) {
        sprintf(errtxt, "%s: char **str", SC_ERROR_ISM_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_ISM, "%s", errtxt);
    }

    /* Read parameter values in strings and check number of values */

    str[0] = strtok(line, "\n\t ");
    if (npar0 > 1) {
        for (i = 1; i < npar0; i++) {
            str[i] = strtok(NULL, "\n\t ");
            if (str[i] == NULL) {
                nparx = i;
                break;
            }
        }
    }

    if (nparx == npar0) {
        str[npar0] = strtok(NULL, "\n\t ");
        if (str[npar0] != NULL) {
            free(str);
            str = NULL;
            sprintf(errtxt, "%s: FILE *stream (more parameters than "
                    "expected)", SC_ERROR_UFS_TXT);
            return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s",
                                         errtxt);
        }
    }

    /* Convert strings in integer and double numbers */

    for (i = 0; i < nparx; i++) {
        strcpy(par[i].c, str[i]);
        strtok(par[i].c, "\r");
        par[i].i = strtol(str[i], NULL, 10);
        if (sc_basic_isnumber(str[i]) == CPL_TRUE) {
            par[i].d = strtod(str[i], NULL);
        } else {
            par[i].d = 0.;
        }
    }

    /* Free memory space of string array */
    free(str);
    str = NULL;

    /* Return number of read parameters */
    *npar = nparx;

    return CPL_ERROR_NONE;
}


double sc_basic_mjd2fracyear(double mjd)
{
    /*!
     * Converts MJD into fractional year (return value, double precision).
     */

    double refjd = 0., ncycle = 0., rest = 0.;
    double frac = 0., nyears = 0., date = 0.;

    refjd = mjd - 51544.;
    ncycle = floor(refjd / 1461.);
    rest = refjd - ncycle * 1461.;

    if (rest < 366.) {
        frac = rest / 366.;
        nyears = 4 * ncycle;
    } else {
        frac = modf((rest - 366.) / 365., &nyears);
        nyears += 4 * ncycle + 1;
    }

    date = 2000. + nyears + frac;

    return date;
}


double sc_basic_fracyear2date(int *year, int *month, int *day,
                              int *hh, int *mm, double *ss,
                              const double *fracyear)
{
    /*!
     * \brief
     *   Brake up a fractional year into year, month, day.
     *
     * This routine brakes up a fractional year into year, month, day, e.g.
     * 2002.380821917808219 corresponds to 20/05/2002, optionally including
     * the time of day. The return value is the rest after braking up the
     * date, i.e. the time of day.
     *
     * \b INPUT:
     * \param fracyear  fractional year (e.g. 2002.345)
     *
     * \b OUTPUT:
     * \param year      year
     * \param month     month
     * \param day       day
     * \param hh        hour
     * \param mm        minute
     * \param ss        second
     *
     * \b RETURN:
     * fractional day (optionally broken up into hh/mm/ss)
     *
     * \b ERRORS:
     * - none
     */

    int i = 0;
    int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    int days_in_year = 365;

    double rest = 0.,
           frac_day = 0.,
           acc = 1e-10; /* corresponds to an accuracy of 0.01s */

    /* calculate the year */
    *year = floor(*fracyear);

    /* calculate remaining fractional year */
    rest = *fracyear - *year;

    /* round to desired accuracy */
    rest = floor(rest/acc + 1) * acc;

    if (*year % 4 == 0 && *year % 100 != 0) {
        /* leap year (February has 29 days) */
        days_in_month[1]++;
        days_in_year++;
    }

    /* calculate remaining fractional days */
    rest *= days_in_year;

    /* count months and subtract days from total */
    for (i = 0; i < 12 && (rest-days_in_month[i]) >= 0; i++) {
        rest -= days_in_month[i];
    }

    /* set month */
    *month = i+1;

    /* calculate days */
    *day = floor(rest) + 1;

    /* calculate remaining fractional day */
    rest -= *day - 1;
    frac_day = rest;

    /* if requested calculate time hh:mm:ss */
    if (hh != NULL) {
        rest *= 24;
        *hh = floor(rest);
        rest -= *hh;

        rest *= 60;
        *mm = floor(rest);
        rest -= *mm;

        rest *= 60;
        *ss = rest;
    }

    return frac_day;
}


cpl_error_code sc_basic_rebin(cpl_table *outspec, const char *outlam,
                              const char *outflux, const cpl_table *inspec,
                              const char *inlam, const char *influx)
{
    /*!
     * Rebins CPL table inspec with columns inlam and influx to wavelength
     * given by outlam in CPL table outspec. The resulting rebinned flux is
     * written to outflux in outspec. If outflux does not exist, it is
     * created. The routine conserves integral of flux.
     *
     * \b INPUT:
     * \param outspec  output spectrum with desired wavelength grid
     * \param outlam   output wavelength column name
     * \param outflux  output flux column name
     * \param inspec   input spectrum
     * \param inlam    input wavelength column name
     * \param influx   input flux column name
     *
     *
     * \b OUTPUT:
     * \param outspec  rebinned input spectrum
     *
     * \b ERRORS:
     * - none
     */

    int n_in = 0, n_out = 0, i = 0, jo = -1, j = 0;
    const double *inlamv = NULL, *influxv = NULL;
    double *outlamv = NULL, *outfluxv = NULL;
    double olmin = 0., olmax = 0., dol = 0., ilmin = 0., ilmax = 0., dil = 0.;
    double rdl = 0.;

    /* Get number of data points */
    n_in = cpl_table_get_nrow(inspec);
    n_out = cpl_table_get_nrow(outspec);

    /* Create flux column in output CPL table if it is not present */
    if (cpl_table_has_column(outspec, outflux) == 0) {
        cpl_table_new_column(outspec, outflux, CPL_TYPE_DOUBLE);
    }

    /* Fill flux column of output spectrum with 0. */
    cpl_table_fill_column_window(outspec, outflux, 0, n_out, 0.);

    /* No data points -> no rebinning */
    if (n_in <= 0 || n_out <= 0) {
        return CPL_ERROR_NONE;
    }

    /* Get pointers to CPL table columns */
    inlamv = cpl_table_get_data_double_const(inspec, inlam);
    influxv = cpl_table_get_data_double_const(inspec, influx);
    outlamv = cpl_table_get_data_double(outspec, outlam);
    outfluxv = cpl_table_get_data_double(outspec, outflux);

    /* One input data point only */
    if (n_in == 1) {
        cpl_table_fill_column_window(outspec, outflux, 0, n_out, influxv[0]);
        return CPL_ERROR_NONE;
    }

    for (i = 0; i < n_out; i++) {

        /* Limits of wavelength bin in output spectrum */

        if (n_out == 1) {

          /* Full range of input spectrum for one output data point */
          olmin = 1.5 * inlamv[0] - 0.5 * inlamv[1];
          olmax = 1.5 * inlamv[n_in-1] - 0.5 * inlamv[n_in-2];

        } else {

            if (i == 0) {
                olmin = 1.5 * outlamv[i] - 0.5 * outlamv[i+1];
            } else {
                olmin = olmax;
            }

            if (i == n_out - 1) {
                olmax = 1.5 * outlamv[i] - 0.5 * outlamv[i-1];
            } else {
                olmax = 0.5 * (outlamv[i] + outlamv[i+1]);
            }

        }

        dol = olmax - olmin;

        do {

            /* Limits of wavelength bin in input spectrum */

            if (j != jo) {

                if (j == 0) {
                    ilmin = 1.5 * inlamv[j] - 0.5 * inlamv[j+1];
                } else {
                    ilmin = ilmax;
                }

                if (j == n_in - 1) {
                    ilmax = 1.5 * inlamv[j] - 0.5 * inlamv[j-1];
                } else {
                    ilmax = 0.5 * (inlamv[j] + inlamv[j+1]);
                }

            }

            /* Effective range of flux value -> weight */

            if (ilmin < olmin && ilmax <= olmax) {
                dil = ilmax - olmin;
            } else if (ilmin >= olmin && ilmax > olmax) {
                dil = olmax - ilmin;
            } else if (ilmin < olmin && ilmax > olmax) {
                dil = olmax - olmin;
            } else {
                dil = ilmax - ilmin;
            }

            if (dil < 0.) {
                dil = 0.;
            }

            /* Average flux of input spectrum in output bin */

            if (dil > 0) {
                rdl = dil / dol;
                outfluxv[i] += influxv[j] * rdl;
            }

            j++;

        } while (ilmax <= olmax && j < n_in);

        j--;
        jo = j;

    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_basic_convolve(cpl_table *spec, const char *colname,
                                 const cpl_array *kernel)
{
    /*!
     * Convolution of column colname in CPL table spec with given kernel.
     *
     * \note The centre of the convolution function is shifted by -0.5 pixels
     *       for an even number of kernel pixels.
     *
     * \b INPUT:
     * \param spec     input spectrum as CPL table
     * \param colname  name of column with input data
     * \param kernel   kernel as CPL array
     *                 (required: sum of all values = 1)
     *
     * \b OUTPUT:
     * \param spec     convolved input spectrum
     *
     * \b ERRORS:
     * - Invalid object value(s)
     */

    cpl_table *tempspec;
    char errtxt[SC_MAXLEN];
    int n = 0, nkpix = 0, k = 0, kmin = 0, kmax = 0, i = 0, j = 0;
    const double *kern = NULL;
    double *flux = NULL, *tflux = NULL;
    double sum = 0., out = 0.;

    /* No data points -> no convolution */

    n = cpl_table_get_nrow(spec);
    if (n <= 0) {
        return CPL_ERROR_NONE;
    }

    /* No kernel data -> no convolution */

    nkpix = cpl_array_get_size(kernel);
    if (nkpix == 0) {
        return CPL_ERROR_NONE;
    }

    /* Get pointer to CPL array */

    kern = cpl_array_get_data_double_const(kernel);

    /* Check kernel */

    for (k = 0; k < nkpix; k++) {
        if (kern[k] < 0 || kern[k] > 1) {
            sprintf(errtxt, "%s: cpl_array *kernel "
                    "(kernel element(s) < 0 or > 1)", SC_ERROR_IOV_TXT);
            return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s",
                                         errtxt);
        }
        sum += kern[k];
    }

    if (sum < 1 - SC_TOL || sum > 1 + SC_TOL) {
        sprintf(errtxt, "%s: cpl_array *kernel (sum of kernel elements != 1)",
                SC_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }

    /* Skip convolution if number of kernel pixels is one */

    if (nkpix == 1) {
        return CPL_ERROR_NONE;
    }

    /* Create temporary spectrum and set flux to 0 in initial spectrum */

    tempspec = cpl_table_duplicate(spec);
    cpl_table_fill_column_window(spec, colname, 0, n, 0.);

    /* Get pointers to CPL table columns */

    flux = cpl_table_get_data_double(spec, colname);
    tflux = cpl_table_get_data_double(tempspec, colname);

    /* Kernel with even or odd pixel number?
       Note: centre of kernel at -0.5 pixels for even pixel number */

    if (nkpix % 2 == 0) {
        kmin = - nkpix / 2;
    } else {
        kmin = - (nkpix - 1) / 2;
    }
    kmax = kmin + nkpix - 1;

    /* Convolve spectrum with kernel */

    for (i = 0; i < n; i++) {

        for (out = 0., k = kmin; k <= kmax; k++) {

            /* Central pixel of kernel */
            j = i - k;

            /* Flux of first or last valid pixel for invalid pixels */
            if (j < 0) {
                j = 0;
            } else if (j >= n) {
                j = n - 1;
            }

            out += tflux[j] * kern[k - kmin];

        }

        flux[i] += out;

    }

    /* Delete temporary spectrum */

    cpl_table_delete(tempspec);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_basic_convolvewindow(cpl_array *convflux,
                                       const cpl_array *flux,
                                       const int range[2],
                                       const cpl_array *kernel)
{
    /*!
     * Convolution of a flux array with given kernel for a fixed range of
     * pixels. The output array has to exist and have the same number of
     * pixels as the input array. Technically, the convolution is carried out
     * by adding the convolved fluxes for each pixel in the selected window.
     *
     * \note The centre of the convolution function is shifted by -0.5 pixels
     *       for an even number of kernel pixels.
     *
     * \b INPUT:
     * \param convflux  CPL array of convolved flux values
     * \param flux      CPL array of unconvolved input flux
     * \param range     range of pixels to be considered (minimum and maximum)
     * \param kernel    kernel as CPL array
     *                  (required: sum of all values = 1)
     *
     * \b OUTPUT:
     * \param convflux  output array with added convolved flux from given
     *                  range
     *
     * \b ERRORS:
     * - Invalid object structure
     * - Invalid input parameter(s)
     * - Invalid object value(s)
     */

    char errtxt[SC_MAXLEN];
    int n = 0, i = 0, nkpix = 0, k = 0, kmin = 0, kmax = 0, jmin = 0;
    int jmax = 0, j = 0;
    const double *influx = NULL, *kern = NULL;
    double *outflux = NULL;
    double sum = 0., in0 = 0., innm1 = 0., in = 0.;

    /* No data points -> no convolution */

    n = cpl_array_get_size(flux);
    if (n <= 0) {
        return CPL_ERROR_NONE;
    }

    /* Check correspondence of data points */

    if (cpl_array_get_size(convflux) != n) {
        sprintf(errtxt, "%s: cpl_array *convflux != cpl_array *flux (size)",
                SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Check range values */

    if (range[0] > range[1]) {
        sprintf(errtxt, "%s: range[2] (min. > max.)", SC_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IIP, "%s", errtxt);
    }
    if (range[0] < 0 || range[1] >= n) {
        sprintf(errtxt, "%s: range[2] (invalid pixel)", SC_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IIP, "%s", errtxt);
    }

    /* Initialise output array if not done so far */

    if (cpl_array_has_invalid(convflux) == 1) {
        cpl_array_fill_window_double(convflux, 0, n, 0.);
    }

    /* Get pointers to CPL arrays */

    influx = cpl_array_get_data_double_const(flux);
    outflux = cpl_array_get_data_double(convflux);

    /* No kernel data -> no convolution */

    nkpix = cpl_array_get_size(kernel);
    if (nkpix == 0) {
        for (i = range[0]; i <= range[1]; i++) {
            outflux[i] += influx[i];
        }
        return CPL_ERROR_NONE;
    }

    /* Get pointer to CPL array */

    kern = cpl_array_get_data_double_const(kernel);

    /* Check kernel */

    for (k = 0; k < nkpix; k++) {
        if (kern[k] < 0 || kern[k] > 1) {
            for (i = range[0]; i <= range[1]; i++) {
                outflux[i] += influx[i];
            }
            sprintf(errtxt, "%s: cpl_array *kernel "
                    "(kernel element(s) < 0 or > 1)", SC_ERROR_IOV_TXT);
            return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s",
                                         errtxt);
        }
        sum += kern[k];
    }

    if (sum < 1 - SC_TOL || sum > 1 + SC_TOL) {
        for (i = range[0]; i <= range[1]; i++) {
            outflux[i] += influx[i];
        }
        sprintf(errtxt, "%s: cpl_array *kernel (sum of kernel elements != 1)",
                SC_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }

    /* Skip convolution if number of kernel pixels is one */

    if (nkpix == 1) {
        for (i = range[0]; i <= range[1]; i++) {
            outflux[i] += influx[i];
        }
        return CPL_ERROR_NONE;
    }

    /* Kernel with even or odd pixel number?
       Note: centre of kernel at -0.5 pixels for even pixel number */

    if (nkpix % 2 == 0) {
        kmin = - nkpix / 2;
    } else {
        kmin = - (nkpix - 1) / 2;
    }
    kmax = kmin + nkpix - 1;

    /* Set pixel range (add virtual pixels for marginal ranges) */

    if (range[0] == 0) {
        jmin = -kmax;
    } else {
        jmin = range[0];
    }
    if (range[1] == n-1) {
        jmax = n-1 - kmin;
    } else {
        jmax = range[1];
    }

    /* Set flux of virtual input pixels */

    in0 = influx[0];
    innm1 = influx[n-1];

    /* Convolve array with kernel */

    for (j = jmin; j <= jmax; j++) {

        /* Flux of real and virtual input pixels */

        if (j < 0) {
            in = in0;
        } else if (j >= n) {
            in = innm1;
        } else {
            in = influx[j];
        }

        /* Calculate output flux for each kernel element and add it to the
           corresponding output pixel */

        for (k = SC_MAX(kmin, -j); k <= SC_MIN(kmax, n-1 - j); k++) {

            i = j + k;

            outflux[i] += in * kern[k - kmin];

        }

    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_basic_filtermedian(cpl_table *spec, const char *colname,
                                     const int npix)
{
    /*!
     * Median filtering of column colname in CPL table spec. The median for
     * each pixel is calculated by considering npix-1 neighbouring pixels.
     *
     * \b INPUT:
     * \param spec     input spectrum as CPL table
     * \param colname  name of column with input data
     * \param npix     number of pixels for median
     *
     * \b OUTPUT:
     * \param spec     filtered input spectrum
     *
     * \b ERRORS:
     * - none
     */

    cpl_array *medpix;
    cpl_table *tempspec;
    int n = 0, urad = 0, j = 0, i = 0;
    double *flux = NULL, *tflux = NULL, *mflux = NULL;

    /* No data points -> no filtering */
    n = cpl_table_get_nrow(spec);
    if (n <= 0) {
        return CPL_ERROR_NONE;
    }

    /* No pixels for median -> no filtering */
    if (npix <= 0) {
        return CPL_ERROR_NONE;
    }

    /* Create temporary spectrum and set flux to 0 in initial spectrum */
    tempspec = cpl_table_duplicate(spec);
    cpl_table_fill_column_window(spec, colname, 0, n, 0.);

    /* Get pointers to CPL table columns */
    flux = cpl_table_get_data_double(spec, colname);
    tflux = cpl_table_get_data_double(tempspec, colname);

    /* Create temporary array for median calculation and initialise it with
       flux of first pixel */
    medpix = cpl_array_new(npix, CPL_TYPE_DOUBLE);
    cpl_array_fill_window(medpix, 0, n, tflux[0]);

    /* Get pointer to CPL array */
    mflux = cpl_array_get_data_double(medpix);

    /* Get upper radius of pixel range */
    urad = (int) floor(0.5 * npix);

    /* Median filtering */

    for (j = 0, i = -urad; i < n; i++) {

        /* New pixel in array for median flux */
        if (i+urad >= n) {
            mflux[j] = tflux[n-1];
        } else {
            mflux[j] = tflux[i+urad];
        }

        /* Calculate median */
        if (i >= 0) {
            flux[i] = cpl_array_get_median(medpix);
        }

        /* New position in array for median flux to be filled */
        j++;
        if (j == npix) {
            j = 0;
        }

    }

    /* Delete temporary spectrum */
    cpl_table_delete(tempspec);
    cpl_array_delete(medpix);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_basic_linecount(long *n_lines, FILE *fp)
{
    /*!
     * \brief
     *   Count the lines in a file.
     *
     * Count the number of lines in a file. File pointer is reset to beginning
     * of file.
     *
     * \b INPUT:
     * \param fp       file pointer
     *
     * \b OUTPUT:
     * \param n_lines  number of lines in file
     *
     * \b ERRORS:
     * - file not found if fseek() returns an error
     */
    char ch='\0';

    *n_lines = 0;

    if (fseek(fp, 0, 0) != 0) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_FOUND,
                                     "Problem with file pointer");
    } else {
        while (ch != EOF) {
            ch = fgetc(fp);
            if (ch == '\n'){
                (*n_lines)++;
            }
        }
        fseek(fp, 0, 0);

        return CPL_ERROR_NONE;
    }
}


void sc_basic_rhum2ppmv(double *ppmv, const double *tem, const double *p,
                        const double *hum)
{
    /*!
     * \brief
     *   Convert relative humidity to ppmv.
     *
     * A single relative humidity value is converted to a volume mixing ratio
     * using the prescription by Murphy and Koop, Review of the vapour
     * pressure of ice and supercooled water for atmospheric applications,
     * Q. J. R. Meteorol. Soc (2005), 131, pp. 1539-1565, for a given
     * temperature and pressure.
     *
     * \note
     *   Valid Temperature range is 123K - 332K.
     *
     * \b INPUT
     * \param tem   Temperature in [Kelvin].
     * \param p     Pressure in [mbar].
     * \param hum   Relative humidity in [%].
     *
     * \b OUTPUT
     * \param ppmv  Volume mixing ratio in ppmv.
     */

    double p_sat,                /* water vapour pressure */
           p_h2o;                /* partial pressure of water in [hPa] */

    double t, logt, log_ew, log_ei, log_e, log_dbl_max = log(DBL_MAX);

    /* ensure that temperature is not absolute zero */
    if (*tem <= 0.) {
        t = DBL_MIN * 2;
        cpl_msg_warning(cpl_func, "Temperatures must be larger than 0K");
    } else {
        t = *tem;
    }

    logt = log(t);

    /* vapour pressure over water */
    log_ew = 54.842763 - 6763.22 / t - 4.210 * logt + 0.000367 * t +
          tanh(0.0415 * (t - 218.8)) * (53.878 - 1331.22 / t - 9.44523 *
          logt + 0.014025 * t);

    /* vapour pressure over hexagonal ice */
    log_ei = 9.550426 - 5723.265 / t + 3.53068 * logt - 0.00728332 * t;

    /* vapour pressure */
    log_e = SC_MIN(log_ew, log_ei);

    /* avoid overflow */
    if (log_e < log_dbl_max) {
        p_sat = exp(log_e) / 100.;
    } else {
        p_sat = DBL_MAX / 100.;
        cpl_msg_warning(cpl_func, "Temperature value results in overflow");
    }

    p_h2o = SC_MIN(*hum, 100) / 100. * p_sat;

    *ppmv = SC_MAX(p_h2o / *p * 1e6, 0);
}


cpl_error_code sc_basic_rhum2ppmv_old(const cpl_array *temp,
                                      const cpl_array *pres,
                                      const cpl_array *rhum,
                                      cpl_array *ppmv)
{
    /*!
     * Convert relative humidity to ppmv.
     *
     * \b INPUT:
     * \param temp  array of temperatures in [K]
     * \param pres  array of pressure values in [mB] or [hPa]
     * \param rhum  array of relative humidity values in [%]
     *
     * \b OUTPUT:
     * \param ppmv  array of ppmv values (absolute humidity)
     *
     * Reference:  Flatau,P.J., R.L.Walko, and W.R.Cotton, 1992:
     *             "Polynomial fits to saturation vapor pressure",
     *             J.Appl.Met., v31, pp1507-1513
     *
     * Valid Temperature range is 183K - 273K (-90C - 0C).
     *
     * \b ERRORS:
     * - for pressure values = 0 a warning is printed and
     *   pressure is assumed = 1e-30.
     */

    cpl_error_code err_code;

    cpl_array *ewater,
              *eice,
              *e,
              *tc,            /* temperature in Deg C */
              *a;

//    double m_h2o = 18.01534,  /* g mol-1 molar mass of water */
//           m_dry =  28.9644,  /* g mol-1 molar mass of dry air */
    double aw[] = { 6.11583699e+00, 4.44606896e-01, 1.43177157e-02,
                    2.64224321e-04, 2.99291081e-06, 2.03154182e-08,
                    7.02620698e-11, 3.79534310e-14,-3.21582393e-16 },
           ai[] = { 6.09868993e+00, 4.99320233e-01, 1.84672631e-02,
                    4.02737184e-04, 5.65392987e-06, 5.21693933e-08,
                    3.07839583e-10, 1.05785160e-12, 1.61444444e-15 },
           val, p;

    long i, j;

    int sz = cpl_array_get_size(temp);

    ewater = cpl_array_new(sz, CPL_TYPE_DOUBLE);
    eice = cpl_array_new(sz, CPL_TYPE_DOUBLE);
    e = cpl_array_new(sz, CPL_TYPE_DOUBLE);
    tc = cpl_array_new(sz, CPL_TYPE_DOUBLE);

    a = cpl_array_new(9, CPL_TYPE_DOUBLE);            /* size of aw[], ai[] */

    for (i = 0; i < sz; i++) {
        cpl_array_set_double(tc, i,
                             cpl_array_get_double(temp, i, NULL) - 273.15);
        cpl_array_set_double(ewater, i, 0);
        cpl_array_set_double(eice, i, 0);
        cpl_array_set_double(e, i, 0);
    }

    for (i = 0; i < 9; i++) {
        cpl_array_set_double(a, i, aw[i]);
    }

    for (i = 8; i >= 0; i--) {
        for (j = 0; j < sz; j++) {
            val = cpl_array_get_double(a, i, NULL) +
                  cpl_array_get_double(ewater, j, NULL) *
                  cpl_array_get_double(tc, j, NULL);
            cpl_array_set_double(ewater, j, val);
        }
    }

    for (i = 0; i < 9; i++) {
        cpl_array_set_double(a, i, ai[i]);
    }

    for (i = 8; i >= 0; i--) {
        for (j = 0; j < sz; j++) {
            val = cpl_array_get_double(a, i, NULL) +
                  cpl_array_get_double(eice, j, NULL) *
                  cpl_array_get_double(tc, j, NULL);
            cpl_array_set_double(eice, j, val);
        }
    }

    for (i = 0; i < sz; i++) {
        val = SC_MIN(cpl_array_get_double(ewater, i, NULL),
                     cpl_array_get_double(eice, i, NULL));
        cpl_array_set_double(e, i, val);
    }

    err_code = CPL_ERROR_NONE;
    for (i = 0; i < sz; i++) {
        p = cpl_array_get_double(pres, i, NULL);
        if (p == 0) {
            p = 1e-30;
            cpl_msg_warning(cpl_func, "sc_rhum2ppmv() encountered a problem:"
                            " Division by zero");
            err_code = CPL_ERROR_DIVISION_BY_ZERO;
        }
        val = cpl_array_get_double(rhum, i, NULL) / 100. *
              cpl_array_get_double(e, i, NULL) / p * 1e6;
        if (val < 0) {
            val = 0.;
        }
        cpl_array_set_double(ppmv, i, val);
    }

    cpl_array_delete(ewater);
    cpl_array_delete(eice);
    cpl_array_delete(e);
    cpl_array_delete(tc);
    cpl_array_delete(a);

    return err_code;
}


cpl_error_code sc_basic_planck(cpl_array *bb, const cpl_array *wavelength,
                               const double temp)
{
    /*!
     * Calculate blackbody radiation.
     *
     * \b INPUT:
     * \param wavelength  wavelength array in micron (1e-6 metres).
     * \param temp        temperature in [K]
     *
     * \b OUTPUT:
     * \param bb          blackbody flux (i.e. PI*Intensity) in [ergs/cm2/s/a]
     *
     * \b ERRORS:
     * - wavelengths cannot be zero.
     */

    long sz, i;

    double c1 = 3.7417749e-5,  /* = 2 * !DPI * h * c * c in cgs units */
           c2 = 1.4387687,     /* = h * c / k            in cgs units */
           log_dbl_max = log(DBL_MAX);

    double wave, val;

    /* ensure same size of wavelength and bb arrays */
    if ((sz = cpl_array_get_size(wavelength)) != cpl_array_get_size(bb)) {
        cpl_array_set_size(bb, sz);
    }

    for (i = 0; i < sz; i++) {
        /* avoid division by zero */
        if ((wave = cpl_array_get_double(wavelength, i, NULL)) == 0) {
            return cpl_error_set_message(cpl_func, CPL_ERROR_DIVISION_BY_ZERO,
                                         "Zero wavelength value");
        }
        wave /= 1e4;
        val = c2 / wave / temp;
        /* avoid overflow */
        if (val < log_dbl_max) {
            cpl_array_set_double(bb, i, c1 / (pow(wave, 5) * (expm1(val)))
                                 * 1e-8);
        } else {
            cpl_array_set_double(bb, i, 0);
        }
    }

    return CPL_ERROR_NONE;
}


void sc_basic_dirslash(char *dir)
{
    /*!
     * Append '/' to path.
     *
     * \b INPUT:
     * \param dir  directory.
     *
     * \b OUTPUT:
     * \param dir  directory.
     *
     * \b ERRORS:
     * - none.
     *
     * \note The input directory gets overwritten.
     */

    if ((unsigned) (strrchr(dir, '/')-dir) != (unsigned) strlen(dir)-1) {
        strncat(dir, "/", 1);
    }
}


void sc_basic_abspath(char *out, const char *dir, const char *cwd)
{
    /*!
     * Expand relative path.
     *
     * \b INPUT:
     * \param dir  relative directory.
     * \param cwd  current working directory.
     *
     * \b OUTPUT:
     * \param out  output directory
     *
     * \b ERRORS:
     * - none.
     *
     * \note If the input directory is already an absolute path, the output
     *       string contains a copy of the input.
     *
     * \note Input strings cannot be longer than SC_MAXLEN.
     */

    char *ptr, *p_out, tmp_dir[SC_MAXLEN] = "";

    // ensure out is empty
    sc_basic_initstring(out, SC_MAXLEN);

    // copy input for using strtok
    strcpy(tmp_dir, dir);

    if (strncmp(tmp_dir, "/", 1) == 0) {
        // if dir is absolute path, out has to be absolute too
        strcat(out, "/");
    } else {
        // if dir is relative path, copy current working directory to out
        strcat(out, cwd);
        sc_basic_dirslash(out);
    }

    // cut up tmp_dir using strtok
    ptr = strtok(tmp_dir, "/");
    while (ptr != NULL) {
        if (strcmp(ptr, "..") == 0) {
            // if ".." move up one level
            if ((p_out = strrchr(out, '/')) != NULL) { // remove trailing "/"
                sprintf(p_out, "%c", '\0');
                if ((p_out = strrchr(out, '/')) != NULL) { // now move up
                    sprintf(p_out, "%c", '\0');
                }
                sc_basic_dirslash(out);
           }
        } else {
            // if "." do nothing
            if (strcmp(ptr, ".") != 0) {
                // else copy current level to out
                strcat(out, ptr);
                strcat(out, "/");
            }
        }
        // get next level
        ptr = strtok(NULL, "/");
    }

    sc_basic_dirslash(out);
}


cpl_error_code sc_basic_access(const char *pathname, const int mode)
{
    /*!
     * \brief
     *   Check for file or directory existence
     *
     * This function provides a wrapper for the intrinsic function access().
     * It reports back all access() errors.
     *
     * \b INPUT:
     * \param pathname  name of a file or directory
     * \param mode      mode for accessibility checks
     *                      F_OK = file existence
     *                      X_OK = file exec permission
     *                      W_OK = file write permission
     *
     * \b OUTPUT:
     * \return          CPL_ERROR_NONE on success, or on failure:
     *                  SC_ERROR_ROFS,
     *                  SC_ERROR_INVAL,
     *                  SC_ERROR_TXTBSY,
     *                  SC_ERROR_ACCES,
     *                  SC_ERROR_LOOP,
     *                  SC_ERROR_NAMETOOLONG,
     *                  SC_ERROR_NOENT,
     *                  SC_ERROR_NOTDIR,
     *                  SC_ERROR_FAULT,
     *                  SC_ERROR_IO,
     *                  SC_ERROR_NOMEM,
     *                  SC_ERROR_UNDEF
     *
     * \b ERRORS:
     * - see access manual for details
     */

    int err;

    errno = 0;

    err = access(pathname, mode);

    switch (errno) {
    case EROFS:
        return cpl_error_set_message(cpl_func, SC_ERROR_ROFS,
                                     "%s: %s", SC_ERROR_ROFS_TXT, pathname);
        break;
    case EINVAL:
        return cpl_error_set_message(cpl_func, SC_ERROR_INVAL,
                                     "%s: %s", SC_ERROR_INVAL_TXT, pathname);
        break;
    case ETXTBSY:
        return cpl_error_set_message(cpl_func, SC_ERROR_TXTBSY,
                                     "%s: %s", SC_ERROR_TXTBSY_TXT, pathname);
        break;
    case EACCES:
        return cpl_error_set_message(cpl_func, SC_ERROR_ACCES,
                                     "%s: %s", SC_ERROR_ACCES_TXT, pathname);
        break;
    case ELOOP:
        return cpl_error_set_message(cpl_func, SC_ERROR_LOOP,
                                     "%s: %s", SC_ERROR_LOOP_TXT, pathname);
        break;
    case ENAMETOOLONG:
        return cpl_error_set_message(cpl_func, SC_ERROR_NAMETOOLONG,
                                     "%s: %s", SC_ERROR_NAMETOOLONG_TXT,
                                     pathname);
        break;
    case ENOENT:
        return cpl_error_set_message(cpl_func, SC_ERROR_NOENT,
                                     "%s: %s", SC_ERROR_NOENT_TXT, pathname);
        break;
    case ENOTDIR:
        return cpl_error_set_message(cpl_func, SC_ERROR_NOTDIR,
                                     "%s: %s", SC_ERROR_NOTDIR_TXT, pathname);
        break;
    case EFAULT:
        return cpl_error_set_message(cpl_func, SC_ERROR_FAULT,
                                     "%s: %s", SC_ERROR_FAULT_TXT, pathname);
        break;
    case EIO:
        return cpl_error_set_message(cpl_func, SC_ERROR_IO,
                                     "%s: %s", SC_ERROR_IO_TXT, pathname);
        break;
    case ENOMEM:
        return cpl_error_set_message(cpl_func, SC_ERROR_NOMEM,
                                     "%s: %s", SC_ERROR_NOMEM_TXT, pathname);
        break;
    }

    if (err != 0) {
        return cpl_error_set_message(cpl_func, SC_ERROR_UNDEF,
                                     "%s: %s", SC_ERROR_UNDEF_TXT, pathname);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_basic_greg2jd(long *jd, const int year, const int month,
                                const int day)
{
    /*!
     * \brief
     *   Convert Gregorian to Julian date
     *
     * This function converts a sequence of year/month/day into a Julian date.
     *
     * \b INPUT:
     * \param year   Gregorian year
     * \param month  Gregorian month
     * \param day    Gregorian day
     *
     * \b OUTPUT:
     * \param jd     Julian date
     *
     * \b ERRORS:
     * - \em SC_ERROR_BADUSERINPUT
     */

    long jd1, jd2, jd3;

    if (year < 1 || month < 1 || month > 12 || day < 1 || day > 31) {
        return SC_ERROR_BADUSERINPUT;
    }

    jd1 = day-32075+1461*(year+4800+(month-14)/12)/4;
    jd2 = 367*(month-2-(month-14)/12*12)/12;
    jd3 = 3*((year+4900+(month-14)/12)/100)/4;

    *jd = jd1+jd2-jd3;

    return CPL_ERROR_NONE;
}


cpl_error_code sc_basic_jd2greg(int *year, int *month, int *day,
                                const long jd)
{
    /*!
     * \brief
     *   Convert Julian to Gregorian date
     *
     * This function converts a Julian date into a sequence of year/month/day.
     *
     * \b INPUT:
     * \param jd     Julian date
     *
     * \b OUTPUT:
     * \param year   Gregorian year
     * \param month  Gregorian month
     * \param day    Gregorian day
     *
     * \b ERRORS:
     * - \em SC_ERROR_BADUSERINPUT
     */

    long l, n, j, y;

    l = jd + 68569;
    n = 4*l/146097;
    l -= (146097*n + 3)/4;
    y = 4000*(l + 1)/1461001;
    l += -1461*y/4 + 31;
    j = 80*l/2447;
    *day = (int)(l - 2447*j/80);
    l = j/11;
    j += 2 - 12*l;
    *year = (int)(100*(n - 49) + y + l);
    *month = (int)j;

    if (*year < 1 || *month < 1 || *month > 12 || *day < 1 || *day > 31) {
        return SC_ERROR_BADUSERINPUT;
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_basic_absfile(char *absfilename, const char *filename)
{
    char cwd[SC_MAXLEN] = "", basedir[SC_MAXLEN] = "", *d = NULL;

    if (filename[0] != '/') {
        /* get current working directory */
        if ((d = getcwd(cwd, sizeof(cwd))) == NULL) {
            return cpl_error_set_message(cpl_func, SC_ERROR_GETCWD,
                                         "%s: %s", SC_ERROR_GETCWD_TXT, cwd);
        }
        sc_basic_dirslash(cwd);
    }

    if ((d = strrchr(filename, '/')) != NULL) {
        // filename conatins path
        strncpy(basedir, filename, (d-filename));
        sc_basic_abspath(absfilename, basedir, cwd);
        strcat(absfilename, (d+1));
    } else {
        // filename does NOT contain path
        strcpy(absfilename, cwd);
        strcat(absfilename, filename);
    }

    return CPL_ERROR_NONE;
}


cpl_boolean sc_basic_parameterlists_same(cpl_parameterlist *list1,
                                         cpl_parameterlist *list2)
{
    cpl_boolean bool = CPL_TRUE;
    cpl_type type;
    cpl_parameter *p1, *p2;

    p1 = cpl_parameterlist_get_first(list1);
    p2 = cpl_parameterlist_get_first(list2);
    type = cpl_parameter_get_type(p1);
    switch (type) {
    case CPL_TYPE_INT:
        if (cpl_parameter_get_int(p1) != cpl_parameter_get_int(p2)) {
            bool = CPL_FALSE;
        }
        break;
    case CPL_TYPE_DOUBLE:
        if (cpl_parameter_get_double(p1) != cpl_parameter_get_double(p2)) {
            bool = CPL_FALSE;
        }
        break;
    case CPL_TYPE_STRING:
        if (cpl_parameter_get_string(p1) != cpl_parameter_get_string(p2)) {
            bool = CPL_FALSE;
        }
        break;
    default:
        break;
    }

    if (bool == CPL_TRUE) {
        while ((p1 = cpl_parameterlist_get_next(list1)) != NULL &&
               (p2 = cpl_parameterlist_get_next(list2)) != NULL) {
            type = cpl_parameter_get_type(p1);
            switch (type) {
            case CPL_TYPE_INT:
                if (cpl_parameter_get_int(p1) !=
                        cpl_parameter_get_int(p2)) {
                    bool = CPL_FALSE;
                }
                break;
            case CPL_TYPE_DOUBLE:
                if (cpl_parameter_get_double(p1) !=
                        cpl_parameter_get_double(p2)) {
                    bool = CPL_FALSE;
                }
                break;
            case CPL_TYPE_STRING:
                if (cpl_parameter_get_string(p1) !=
                        cpl_parameter_get_string(p2)) {
                    bool = CPL_FALSE;
                }
                break;
            default:
                break;
            }
            if (bool == CPL_FALSE) {
                break;
            }
        }
    }
    return bool;
}

cpl_error_code sc_basic_status2txt(char *msg, const int status)
{
    // check MPFIT status
    switch (status) {
    // success
    case MP_OK_CHI:
        sprintf(msg, "Convergence in chi-square value");
        break;
    case MP_OK_PAR:
        sprintf(msg, "Convergence in parameter value");
        break;
    case MP_OK_BOTH:
        sprintf(msg, "Both MP_OK_PAR and MP_OK_CHI hold");
        break;
    case MP_OK_DIR:
        sprintf(msg, "Convergence in orthogonality");
        break;
    case MP_MAXITER:
        sprintf(msg, "Maximum number of iterations reached");
        break;
    case MP_FTOL:
        sprintf(msg, "ftol is too small; no further improvement");
        break;
    case MP_XTOL:
        sprintf(msg, "xtol is too small; no further improvement");
        break;

    case MP_GTOL:
        sprintf(msg, "gtol is too small; no further improvement");
        break;
    // failure
    case MP_ERR_INPUT:
        sprintf(msg, "General input parameter error");
        break;
    case MP_ERR_NAN:
        sprintf(msg, "User function produced non-finite values");
        break;
    case MP_ERR_FUNC:
        sprintf(msg, "No user function was supplied");
        break;
    case MP_ERR_NPOINTS:
        sprintf(msg, "No user data points were supplied");
        break;
    case MP_ERR_NFREE:
        sprintf(msg, "No free parameters");
        break;
    case MP_ERR_MEMORY:
        sprintf(msg, "Memory allocation error");
        break;
    case MP_ERR_INITBOUNDS:
        sprintf(msg, "Initial values inconsistent w constraints");
        break;
    case MP_ERR_BOUNDS:
        sprintf(msg, "Initial constraints inconsistent");
        break;
    case MP_ERR_PARAM:
        sprintf(msg, "General input parameter error");
        break;
    case MP_ERR_DOF:
        sprintf(msg, "Not enough degrees of freedom");
        break;
    default:
        sprintf(msg, "Exception: unknown exit value");
        return SC_ERROR_INVAL;
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_basic_clipmean(double *mean, double *rms,
                                 cpl_array *arr, const cpl_boolean clip) {
    /*!
     * \brief
     *   Calculate mean and rms using sigma-clipping routine
     *
     * This function calculates the mean and rms of a cpl_array first by
     * computing the median absolute difference and afterwards applying
     * Huber's method. Values identified as outliers are flagged "invalid"
     * in the input array. If a more stringent value on the mean/rms is
     * desired, /em clip can be set to 1, resulting in the outliers being
     * removed from the calculation of the mean/rms (in contrast to
     * "compressing" them to the outer boundaries of the clipping, i.e.
     * Huber's method). See: AMC Technical Briefs - ISSN 1757-5958, No.6,
     * Royal Society of Chemistry, 2001
     * http://www.rsc.org/images/brief6_tcm18-25948.pdf
     *
     * \note
     *   The input array might have invalid rows on output. Invalid rows on
     *   input are excluded from the computation.
     *
     * \b INPUT:
     * \param arr   cpl_array
     * \param clip  flag for clipping mode
     *
     * \b OUTPUT:
     * \param mean  mean
     * \param rms   root-mean-square
     *
     * \b ERRORS:
     * - none
     */

    cpl_array *mad_arr;    // median absolute difference

    double val_p, val_m, val = 0., mean_old, eps = 1e-4;

    int i, j,              // loop variables
        n,                 // size of arr
        nclip = 0, nclip_old = 0;         // total number of clipped elements

    // get size of array
    n = cpl_array_get_size(arr);

    // initialise array for median absolute difference
    mad_arr = cpl_array_new(n, CPL_TYPE_DOUBLE);

    // compute first guess for mean and rms using median absolute difference
    *mean = cpl_array_get_median(arr);
    for (i = 0; i < n; i++) {
        // compute MAD
        cpl_array_set_double(mad_arr, i,
                            fabs(cpl_array_get_double(arr, i, NULL) - *mean));
        if (cpl_array_is_valid(arr, i) == 0) {
            // remember invalid values from arr in mad_arr
            cpl_array_set_invalid(mad_arr, i);
        }
    }
    *rms = cpl_array_get_median(mad_arr) * 1.5;

    // now use Huber's method to improve estimate

    // reset mad_arr
    for (i = 0; i < n; i++) {
        cpl_array_set_double(mad_arr, i, cpl_array_get_double(arr, i, NULL));
        if (cpl_array_is_valid(arr, i) == 0) {
            // remember invalid values from arr in mad_arr
            cpl_array_set_invalid(mad_arr, i);
        }
    }

    // do at max 100 iterations (should converge much faster)
    for (j = 0, mean_old = *mean, nclip = 0, nclip_old = 0; j < 100; j++) {
        // re-compute outliers
        val_p = *mean + (1.5 * *rms);
        val_m = *mean - (1.5 * *rms);
        for (i = 0, *mean = 0; i < n; i++) {
            if (cpl_array_is_valid(mad_arr, i) == 1) {
                val = cpl_array_get_double(mad_arr, i, NULL);
                if (val > val_p) {
                    // set outliers invalid in original array
                    cpl_array_set_invalid(arr, i);
                    cpl_array_set_double(mad_arr, i, val_p);
                    nclip++;
                } else if (val < val_m) {
                    // set outliers invalid in original array
                    cpl_array_set_invalid(arr, i);
                    cpl_array_set_double(mad_arr, i, val_m);
                    nclip++;
                }
            }
        }
        *mean = cpl_array_get_mean(mad_arr);
        *rms = cpl_array_get_stdev(mad_arr)*1.134;

        if (fabs(mean_old / *mean - 1) < eps || nclip == nclip_old) {
            break;
        }
        mean_old = *mean;
        nclip_old = nclip;
    }

    if (clip == CPL_TRUE) {
        *mean = cpl_array_get_mean(arr);
        *rms = cpl_array_get_stdev(arr);
    }

    cpl_array_delete(mad_arr);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_basic_col2arr(cpl_array *arr, const cpl_table *tab,
                                const char *colname)
{
    /*!
     * \brief
     * Routine to copy a CPL table column to a CPL array.
     *
     * The target array must have size = 0, otherwise existing array data will
     * be overwritten. Data types of source array and target table column
     * must match.
     *
     * Supported data types: CPL_TYPE_INT, CPL_TYPE_FLOAT, CPL_TYPE_DOUBLE,
     *                       CPL_TYPE_STRING
     *
     * \b INPUT:
     * \param arr       (empty) CPL array to be filled
     * \param tab       CPL table containing selected column 'colname'
     * \param colname   Column name of table to be copied to array
     *
     * \b OUTPUT
     * \param arr       Filled CPL array
     *
     */

    int nrow_tab = 0;               /* # of rows in input table             */
    int nrow_arr = 0;               /* # of rows in input table             */
    int loop = 0;               /* Loop variable                            */
    char err_msg[SC_MAXLEN]; /* Error message string                     */

    cpl_type tabcoltype;        /* Type of input table column               */
    cpl_type arrcoltype;        /* Type of target array                     */

/*--------------------------------------------------------------------------*/

    /* Checking validity of input */
    nrow_tab   = cpl_table_get_nrow(tab);
    nrow_arr   = cpl_array_get_size(arr);
    tabcoltype = cpl_table_get_column_type(tab, colname);
    arrcoltype = cpl_array_get_type(arr);

    /* Input table size check */
    if (nrow_tab == 0) {
        sprintf(err_msg, "CPL input table has size = 0!\n"
                        "Cannot continue. Emergency stop.");
        return cpl_error_set_message(cpl_func,  CPL_ERROR_ILLEGAL_INPUT,
                                    "%s",  err_msg);
    }

    /* Existence of source column ? */
    if (cpl_table_has_column(tab, colname) == 0) {
        sprintf(err_msg, "CPL input table does not contain a column '%s'!\n"
                        "Cannot continue. Emergency stop.",colname);
        return cpl_error_set_message(cpl_func,  CPL_ERROR_ILLEGAL_INPUT,
                                    "%s",  err_msg);
    }

    /* Size of target array >0 ? */
    if (nrow_arr > 0) {
        cpl_msg_warning(cpl_func,"Input CPL array size > 0, existing data "
                                 "will be overwritten!");
    }

    /* CPL types match ? */
    if (tabcoltype != arrcoltype) {
        sprintf(err_msg, "CPL types of table column and array do not match!\n"
                        "Cannot continue. Emergency stop.");
        return cpl_error_set_message(cpl_func,  CPL_ERROR_ILLEGAL_INPUT,
                                    "%s",  err_msg);
    }

    /* Copying data */
    cpl_array_set_size(arr, nrow_tab);

    for (loop = 0; loop < nrow_tab; loop++)
    {
        switch (tabcoltype) {
        case CPL_TYPE_INT:
            cpl_array_set_int(arr, loop,
                              cpl_table_get_int(tab,colname,loop,NULL));
            break;
        case CPL_TYPE_FLOAT:
            cpl_array_set_float(arr, loop,
                              cpl_table_get_float(tab,colname,loop,NULL));
            break;
        case CPL_TYPE_DOUBLE:
            cpl_array_set_double(arr, loop,
                              cpl_table_get_double(tab,colname,loop,NULL));
            break;
        case CPL_TYPE_STRING:
            cpl_array_set_string(arr, loop,
                              cpl_table_get_string(tab,colname,loop));
            break;
        default:
            sprintf(err_msg, "CPL type not supported!\n"
                             "Cannot continue. Emergency stop.");
            return cpl_error_set_message(cpl_func,  CPL_ERROR_ILLEGAL_INPUT,
                                        "%s",  err_msg);
            break;
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_basic_arr2col(cpl_table *tab, const char *colname,
                                const cpl_array *arr)
{
    /*!
     * \brief
     * Routine to copy a CPL array to a column of a CPL table.
     *
     * If the target column does not exist it is created and subsequently
     * filled with array data. Also the sizes of the source array and target
     * CPL table must match, except an empty target table is used. In this
     * case, the size is set to the inout CPL array size.
     *
     * Data types of target column and source array must coincide
     *
     * Supported data types: CPL_TYPE_INT, CPL_TYPE_FLOAT, CPL_TYPE_DOUBLE,
     *                       CPL_TYPE_STRING
     *
     *
     * \b INPUT:
     * \param arr       CPL array to be copied to target column
     * \param tab       CPL table containing selected column
     * \param colname   Name of CPL table target column
     *
     * \b OUTPUT
     * \param tab       Target CPL table
     * \param colname   Target column of CPL table to be filled with array
     *                  data
     *
     */

    int nrow_arr = 0;           /* # of rows in input array                 */
    int nrow_tab = 0;           /* # of rows in target table                */
    int loop = 0;               /* Loop variable                            */
    char err_msg[SC_MAXLEN]; /* Error message string                     */

    cpl_type arrcoltype;        /* Type of target array                     */

/*--------------------------------------------------------------------------*/

    /* Checking validity of input */
    arrcoltype = cpl_array_get_type(arr);
    nrow_arr = cpl_array_get_size(arr);
    nrow_tab = cpl_table_get_nrow(tab);

    /* Empty array ? */
    if (nrow_arr == 0) {
        sprintf(err_msg, "CPL array size = 0\n"
                         "Cannot continue. Emergency stop.");
        return cpl_error_set_message(cpl_func,  CPL_ERROR_ILLEGAL_INPUT,
                                    "%s",  err_msg);
    }

    /* Sizes match? */
    if ((nrow_arr != nrow_tab) && (nrow_tab > 0)) {
        sprintf(err_msg, "CPL array and table size do not coincide!\n"
                         "Cannot continue. Emergency stop.");
        return cpl_error_set_message(cpl_func,  CPL_ERROR_ILLEGAL_INPUT,
                                    "%s",  err_msg);
    }


    /* Checking / creating column */
    if (cpl_table_has_column(tab, colname) == 0) {
        arrcoltype = cpl_array_get_type(arr);
        cpl_table_new_column(tab, colname, arrcoltype);
        if (nrow_tab == 0) {
            cpl_table_set_size(tab,nrow_arr);
        }
    }

    /* Copying info */
    for (loop = 0; loop<nrow_arr; loop++) {
        switch (arrcoltype) {
        case CPL_TYPE_INT:
            cpl_table_set_int(tab, colname, loop,
                              cpl_array_get_int(arr, loop, NULL));
            break;
        case CPL_TYPE_FLOAT:
            cpl_table_set_float(tab, colname, loop,
                              cpl_array_get_float(arr, loop, NULL));
            break;
        case CPL_TYPE_DOUBLE:
            cpl_table_set_double(tab, colname, loop,
                              cpl_array_get_double(arr, loop, NULL));
            break;
        case CPL_TYPE_STRING:
            cpl_table_set_string(tab, colname, loop,
                              cpl_array_get_string(arr, loop));
            break;
        default:
            sprintf(err_msg, "CPL type not supported!\n"
                             "Cannot continue. Emergency stop.");
            return cpl_error_set_message(cpl_func,  CPL_ERROR_ILLEGAL_INPUT,
                                        "%s", err_msg);
            break;
        }
    }

    return CPL_ERROR_NONE;
}


int sc_basic_sortarr_double(const void *p1, const void *p2)
{
    if (*(double *)p1 > *(double *)p2) {
        return 1;
    } else if (*(double *)p1 < *(double *)p2) {
        return -1;
    } else {
        return 0;
    }
}


cpl_error_code sc_basic_sortarr(cpl_array *arr, cpl_type type)
{
    /*!
     * \brief
     * Routine to sort a CPL array from lowest to largest value.
     *
     * Array data gets sorted using qsrt algorithm.
     *
     * Supported data types: CPL_TYPE_DOUBLE
     *
     * \note
     *   Array data gets changed on output (sorted)
     *
     * \b INPUT:
     * \param arr       CPL array to sort
     * \param type      CPL data type
     *
     * \b OUTPUT:
     * \param arr       sorted CPL array
     *
     * \b ERRORS:
     * - CPL_ERROR_DATA_NOT_FOUND:  arr has size zero
     * - CPL_ERROR_INVALID_TYPE:    data type not supported
     *
     */
    void *p = NULL;

    if (cpl_array_get_size(arr) == 0) {
        return CPL_ERROR_DATA_NOT_FOUND;
    }

    switch (type) {
        // Double array
        case CPL_TYPE_DOUBLE:
            // get pointer to data
            p = cpl_array_get_data_double(arr);

            // sort data
            qsort(p, cpl_array_get_size(arr),
                  cpl_type_get_sizeof(CPL_TYPE_DOUBLE),
                  sc_basic_sortarr_double);

            return CPL_ERROR_NONE;
            break;
        default:
            // type not supported
            return CPL_ERROR_INVALID_TYPE;
            break;
    }
}


cpl_error_code sc_basic_copytable_content(cpl_table *outtab,
                                          const cpl_table *intab)
{
    /*!
     * Copies the content of a table. The output table must exist and have the
     * same structure as the input table. Only columns with CPL_TYPE_STRING,
     * CPL_TYPE_INT, CPL_TYPE_FLOAT, and CPL_TYPE_DOUBLE can be handled.
     *
     * \b INPUT:
     * \param intab   input table
     *
     * \b OUTPUT:
     * \param outtab  output table with copied content of input table
     *
     * \b ERRORS:
     * - Invalid object structure
     */

    cpl_type type;
    cpl_array *colnames;
    char errtxt[SC_MAXLEN], colname[SC_LENLINE+1];
    const char **sp;
    int ncol = 0, i = 0;
    const int *ip;
    const float *fp;
    const double *dp;

    if (cpl_table_compare_structure(outtab, intab) != 0) {
        sprintf(errtxt, "%s: structure of cpl_table *outtab != *intab",
                SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }

    /* Get column labels */
    colnames = cpl_table_get_column_names(intab);
    ncol = cpl_array_get_size(colnames);

    /* Copy content of columns */

    for (i = 0; i < ncol; i++) {

        strncpy(colname, cpl_array_get_string(colnames, i), SC_LENLINE+1);
        type = cpl_table_get_column_type(intab, colname);

        if (type == CPL_TYPE_STRING) {
            sp = cpl_table_get_data_string_const(intab, colname);
            cpl_table_copy_data_string(outtab, colname, sp);
        } else if (type == CPL_TYPE_INT) {
            ip = cpl_table_get_data_int_const(intab, colname);
            cpl_table_copy_data_int(outtab, colname, ip);
        } else if (type == CPL_TYPE_FLOAT) {
            fp = cpl_table_get_data_float_const(intab, colname);
            cpl_table_copy_data_float(outtab, colname, fp);
        } else if (type == CPL_TYPE_DOUBLE) {
            dp = cpl_table_get_data_double_const(intab, colname);
            cpl_table_copy_data_double(outtab, colname, dp);
        }

    }

    /* Free memory */
    cpl_array_delete(colnames);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_basic_copytable_full(cpl_table *outtab,
                                       const cpl_table *intab)
{
    /*!
     * Copies the structure and content of a table. The output table must
     * exist. Missing columns are created. Only columns with CPL_TYPE_STRING,
     * CPL_TYPE_INT, CPL_TYPE_FLOAT, and CPL_TYPE_DOUBLE can be handled.
     *
     * \b INPUT:
     * \param outtab  empty output table
     * \param intab   input table
     *
     * \b OUTPUT:
     * \param outtab  output table with copied content of input table
     *
     * \b ERRORS:
     * - none
     */

    cpl_type type;
    cpl_array *colnames;
    char colname[SC_LENLINE+1];
    const char **sp;
    int nrow = 0, ncol = 0, i = 0;
    const int *ip;
    const float *fp;
    const double *dp;

    /* Set size of output table */
    nrow = cpl_table_get_nrow(intab);
    cpl_table_set_size(outtab, nrow);

    /* Get column labels */
    colnames = cpl_table_get_column_names(intab);
    ncol = cpl_array_get_size(colnames);

    /* Copy content of columns and create columns if necessary */

    for (i = 0; i < ncol; i++) {

        strncpy(colname, cpl_array_get_string(colnames, i), SC_LENLINE+1);
        type = cpl_table_get_column_type(intab, colname);
        if (cpl_table_has_column(outtab, colname) != 1) {
            cpl_table_new_column(outtab, colname, type);
        }

        if (type == CPL_TYPE_STRING) {
            sp = cpl_table_get_data_string_const(intab, colname);
            cpl_table_copy_data_string(outtab, colname, sp);
        } else if (type == CPL_TYPE_INT) {
            ip = cpl_table_get_data_int_const(intab, colname);
            cpl_table_copy_data_int(outtab, colname, ip);
        } else if (type == CPL_TYPE_FLOAT) {
            fp = cpl_table_get_data_float_const(intab, colname);
            cpl_table_copy_data_float(outtab, colname, fp);
        } else if (type == CPL_TYPE_DOUBLE) {
            dp = cpl_table_get_data_double_const(intab, colname);
            cpl_table_copy_data_double(outtab, colname, dp);
        }

    }

    /* Free memory */
    cpl_array_delete(colnames);

    return CPL_ERROR_NONE;
}


void sc_basic_initstring(char *str, const long n)
{
    /*!
     * \brief
     *   Initialise a string variable.
     *
     * This function initialises a given string \em str of length \em n
     * with "\0".
     *
     * \b INPUT:
     * \param n    length of string
     *
     * \b OUTPUT:
     * \param str  string
     *
     * \return     mothing
     */

    long i;

    for (i = 0; i < n; i++) {
        str[i] = '\0';
    }
}


void sc_basic_terminatestring(char *str)
{
    /*!
     * \brief
     *   Put "\0" at end of string.
     *
     * This function places a "\0" at end of the input string.
     *
     * \note
     *   No checks are performed for writing beyond the boundary of allocated
     *   memory for the input string.
     *
     * \b INPUT & OUTPUT:
     * \param str  string
     *
     * \b OUTPUT:
     * \return     nothing
     */

    char *p;                   /* character pointer */

    p = str + strlen(str);     /* set at nominal '\0' of input string */

    if (*p != '\0') {          /* if not terminated properly */
        *p = '\0';             /* place string termination character */
    }
}


char *sc_basic_strtrim(char *str)
{
    /*!
     * \brief
     *   Remove leading and trailing blanks from string.
     *
     * This function removes all leading and trailing " " characters from
     * \em str using \c isspace().
     *
     * \b INPUT:
     * \param str  string
     *
     * \b OUTPUT:
     * \return     string
     */

    if (str != NULL && *str != 0) {
        int i = 0,                  /* counter for start position of string */
            j = 0,                  /* counter for end position of string */
            len = strlen(str);      /* length of input string */

        char *p, *q,       /* pointers for looping through the input string */
            *out = (char *)malloc(len + 1);           /* new output string */

        /* set pointers to beginning of input string */
        q = str;
        p = out;

        /* skip over leading spaces */
        while (isspace(*q)) {
            i++;
            q++;
        }
        /* i now has start position of string */

        /* skip over trailing spaces */
        q = str + len - 1;
        while (isspace(*q)) {
            j++;
            q--;
        }
        /* j now has end position of string */

        j = len - j - i;            /* count number of remaining characters */
        q = str + i;          /* set pointer to beginning of trimmed string */
        for (; j>0; j--) {
            *p = *q;                                     /* copy characters */
            p++;
            q++;
        }

        *p = '\0';                                      /* terminate string */

        return(out);
    } else {
        return(NULL);
    }
}


void sc_basic_strtrim_inplace(char *str)
{
    /*!
     * \brief
     *   Remove leading and trailing blanks from string.
     *
     * This function removes all leading and trailing " " characters from
     * \em str using \c isspace(). In contrast to \c sc_basic_strtrim(), the
     * action is performed in place, i.e. the input gets overwritten.
     *
     * \b INPUT & OUTPUT:
     * \param str  string
     *
     * \b OUTPUT:
     * \return     nothing
     */

    int i = 0,                  /* counter for start position of string     */
        j = 0,                  /* counter for end position of string       */
        len = 0;                /* length of input string                   */

    char *p, *q;           /* pointers for looping through the input string */

    if (str != NULL && *str != 0) {

        len = strlen(str);
        /* set pointers to beginning of input string */
        p = str;
        q = str;

        /* skip over leading spaces */
        while (isspace(*q)) {
            i++;
            q++;
        }
        /* i now has start position of string */

        /* skip over trailing spaces */
        q = str + len - 1;
        while (isspace(*q)) {
            j++;
            if (*q == *p)   /* Avoids valgrind error                        */
            {
                break;
            }
            q--;
        }
        /* j now has end position of string */

        j = len - j - i;            /* count number of remaining characters */
        q = str + i;          /* set pointer to beginning of trimmed string */
        for (; j>0; j--) {
            *p = *q;                                     /* copy characters */
            p++;
            q++;
        }

        *p = '\0';                                      /* terminate string */
    }
}


cpl_boolean sc_basic_isnumber(char *str)
{
    /*!
     * \brief
     *   Check if string contains number.
     *
     * This function checks whether \em str contains a valid number. A leading
     * "+" or "-" is treated as a sign. A single "." is identified as a
     * decimal point. Finally, an "e" is accepted if an exponent follows this
     * letter.
     *
     * \note
     *   Surrounding spaces are not treated.
     *
     * \b INPUT:
     * \param str  string
     *
     * \b OUTPUT:
     * \return     CPL_TRUE if \em str is number, CPL_FALSE else.
     */

    cpl_boolean issign = CPL_FALSE, ispoint = CPL_FALSE, isexp = CPL_FALSE;
    int length = 0, iexp = 0, i = 0;

    if (str && *str) {

        length = strlen(str);

        for (iexp = 0, i = 0; i < length; i++) {

            if (isdigit(str[i]) != 0) {
                /* Numbers are OK */
                continue;
            } else if (str[i] == '+' || str[i] == '-') {
                /* Correct position of sign */
                if (i == 0) {
                    issign = CPL_TRUE;
                }
                if (i == length-1) {
                    return CPL_FALSE;
                } else if (i > 1 && iexp != 1) {
                    return CPL_FALSE;
                }
            } else if (str[i] == '.') {
                /* Correct position of decimal point */
                if (i == 0 && length == 1) {
                    return CPL_FALSE;
                } else if (i > 0 && ispoint == CPL_TRUE) {
                    return CPL_FALSE;
                } else if (i > 1 && iexp > 0) {
                    return CPL_FALSE;
                }
                ispoint = CPL_TRUE;
            } else if (str[i] == 'e' || str[i] == 'E') {
                /* Correct position of exponent */
                if (i == length-1) {
                    return CPL_FALSE;
                } else if (i == 1 &&
                           (issign == CPL_TRUE || ispoint == CPL_TRUE)) {
                    return CPL_FALSE;
                } else if (i > 0 && isexp == CPL_TRUE) {
                    return CPL_FALSE;
                }
                isexp = CPL_TRUE;
            } else {
                return CPL_FALSE;
            }

            /* Count number of characters after exponent sign */
            if (isexp == CPL_TRUE) {
                iexp++;
            }

        }

    } else {

        return CPL_FALSE;  // null pointer

    }

    return CPL_TRUE;
}


cpl_boolean sc_basic_isinteger(char *str)
{
    /*!
     * \brief
     *   Check if string contains an integer only.
     *
     * This function checks whether \em str contains a valid integer. A
     * leading "+" or "-" is treated as a sign.
     *
     * \note
     *   Surrounding spaces are not treated.
     *
     * \b INPUT:
     * \param str  String.
     *
     * \b OUTPUT:
     * \return     CPL_TRUE if \em str is integer, CPL_FALSE else.
     */

     int i;
     int falseflag = 0;

    /* include "-" and "+" as sign */

    if (str && *str) {
        if (str[0] == '+' || str[0] == '-') {
            i = 1;
        } else {
            i = 0;
        }

        for (;(unsigned)i < strlen(str); i++) {
            if (isdigit(str[i]) == 0) {
                falseflag = 1;
            }
        }
    }

    if (falseflag == 1) {
        return CPL_FALSE;
    } else {
        return CPL_TRUE;
    }

}


cpl_error_code sc_basic_interpollin(const double *x_out, double *y_out,
                                    const long n_out, const double *x_ref,
                                    const double *y_ref, const long n_ref,
                                    const int extrapolate)
{
    /*!
     * \brief
     *   Linear interpolation.
     *
     * This function calculates the interpolated y-values \em y_out at the
     * positions \em x_out with respect to the reference x/y pairs
     * (\em x_ref / \em y_ref).
     *
     * \note
     *   If \em extrapolate is set to 1, points outside the range of the
     *   reference vectors will be extrapolated based on the first / last
     *   values in the reference vectors for the low / high end, respectively.
     *   Otherwise, out of bounds values are replaced with the outermost
     *   points.
     *
     * \b INPUT:
     * \param x_out        Desired output x-spacing.
     * \param n_out        Length of \em x_out / \em y_out.
     * \param x_ref        Reference x-spacing.
     * \param y_ref        Reference y-values at \em x_ref.
     * \param n_ref        Length of \em x_ref / \em y_ref.
     * \param extrapolate  Extrapolation flag
     *
     * \b OUTPUT:
     * \param y_out        Requested y-values at \em x_out.
     *
     * \return             CPL_ERROR_NONE on success,
     *                     CPL_ERROR_DIVISION_BY_ZERO else.
     */

    long out, ref;     /* counter for length of output & reference array */

    const double *xr,
                 *yr;  /* pointers for looping through the reference arrays */
    const double *xo;
    double *yo;        /* pointers for looping through the output arrays */

    int divide_by_zero = 0; /* at least one value has a divide by zero */

    /* set pointers to start values */
    xr = x_ref;
    yr = y_ref;
    xo = x_out;
    yo = y_out;

    for (ref = 0, out = 0; out < n_out; out++, xo++, yo++) {
        /*
         *  find first element in ref that is larger than current out and
         *  ensure that ref does not overshoot valid range
         */
        while (*xr < *xo && ref < n_ref) {
            xr++;
            yr++;
            ref++;
        }

        /* if current out is larger than ref go one step back */
        if (xr-x_ref > 0) {
            xr--;
            yr--;
            ref--;
        }

        /* if out is larger than last ref use last two points */
        if (ref == n_ref - 1) {
            if (extrapolate == 1) {
                xr = x_ref + n_ref - 2;
                yr = y_ref + n_ref - 2;
            } else {
                xr = x_ref + n_ref - 1;
                yr = y_ref + n_ref - 1;
            }
        }

        /* calculate linear interpolation */
        if (*(xr + 1) - *xr == 0) {
            *yo = *yr;
            divide_by_zero = 1;
        } else {
            /* if out is smaller than last ref use first point */
            if (extrapolate != 1 && *xo < *xr) {
                *yo = *yr;
            } else {
                /* else use linear interpolation */
                *yo = (*xo - *xr) * (*(yr + 1) - *yr) /
                      (*(xr + 1) - *xr) + *yr;
            }
        }
    }

    if (divide_by_zero == 1) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_DIVISION_BY_ZERO,
                                     "Duplicate x-values");
    } else {
        return CPL_ERROR_NONE;
    }
}


inline int sc_basic_gaussfunc(double *fgauss, const double *xgauss,
                              const int n_data, const double *par)
{
    /*!
     * Computes Gaussian.
     *
     * \b INPUT:
     * \param xgauss  vector of x positions
     * \param n_data  number of data points
     * \param par     parameters for Gaussian
     *
     * \b OUTPUT:
     * \param fgauss  resulting Gaussian
     *
     * \b ERRORS:
     * - none
     */

    // p[0]: Normalisation
    // p[1]: position
    // p[2]: width

    int i;

    /* Compute Gaussian */
    for (i = 0; i < n_data; i++) {
        fgauss[i] = par[0] * pow(CPL_MATH_E, -pow(xgauss[i] - par[1], 2) /
                                 (2 * pow(par[2], 2)));
    }

    return 0;
}


cpl_error_code sc_basic_getfilename(char *dir, char *filename, char *suffix,
                                    const char *path)
{
    /*!
     * Splits full path into directory, file name, and file name suffix.
     * Returns NULL if a component cannot be found in the input string.
     *
     * \b INPUT:
     * \param path      path including file name
     *
     * \b OUTPUT:
     * \param dir       directory path (without slash at the end)
     * \param filename  file name
     * \param suffix    file identifier like "asc" or "fits"
     *
     * \b ERRORS:
     * - none
     */

    char str[SC_MAXLEN], *ptr;
    cpl_boolean hasdot = CPL_FALSE, hasslash = CPL_FALSE;
    size_t len;
    int i = 0;

    /* Handle NULL input */
    if (path == NULL || path[0] == '\0') {
        dir = NULL;
        filename = NULL;
        suffix = NULL;
        return CPL_ERROR_NONE;
    }

    /* Get string length */
    len = strlen(path);

    /* Get copy of input string */
    strcpy(str, path);

    /* Find separators and substitute them by spaces */
    for (i = len-1; i >= 0; i--) {
        if (str[i] == '.' && hasdot == CPL_FALSE && hasslash == CPL_FALSE) {
            str[i] = ' ';
            hasdot = CPL_TRUE;
        } else if (str[i] == '/' && hasslash == CPL_FALSE) {
            str[i] = ' ';
            hasslash = CPL_TRUE;
        }
    }

    /* Handle the case without path and file suffix */
    if (hasdot == CPL_FALSE && hasslash == CPL_FALSE) {
        if (dir != NULL) {
            sprintf(dir, ".");
        }
        if (filename != NULL) {
            strcpy(filename, str);
        }
        suffix = NULL;
        return CPL_ERROR_NONE;
    }

    /* Split string at spaces by means of strtok and return substrings */

    /* Separation of directory path and file name */
    ptr = strtok(str, " ");
    if (hasslash == CPL_FALSE) {
        if (dir != NULL) {
            sprintf(dir, ".");
        }
        if (filename != NULL) {
            sprintf(filename, "%s", ptr);
        }
    } else {
        if (dir != NULL) {
            sprintf(dir, "%s", ptr);
        }
        ptr = strtok(NULL, " ");
        if (filename != NULL) {
            sprintf(filename, "%s", ptr);
        }
    }

    /* Separation of file name and suffix */
    if (hasdot == CPL_TRUE && suffix != NULL) {
        ptr = strtok(NULL, " ");
        sprintf(suffix, "%s", ptr);
    } else {
        suffix = NULL;
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_basic_getmaskval_vector(double maskval[2],
                                          const cpl_vector *vector)
{
    /*!
     * Gets meaning of mask values provided by a CPL vector. It is assumed
     * that a good value equals 1 if only 1 and possiby 0 are present. If only
     * 0 exists, the definition is reversed. This is also true for a minimum
     * mask value of 0 and any larger number. The latter will be classified
     * as bad mask values.
     *
     * \b INPUT:
     * \param vector   input CPL vector
     *
     * \b OUTPUT:
     * \param maskval  mask values (1st: bad value, 2nd: good value)
     *
     * \b ERRORS:
     * - No data
     * - Invalid object value(s)
     */

    char errtxt[SC_MAXLEN] = "";
    int n = 0;
    double maskmin = 0., maskmax = 0.;

    /* Get size of vector */
    n = cpl_vector_get_size(vector);
    if (n <= 0) {
        sprintf(errtxt, "%s: cpl_vector *vector", SC_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_NDA, "%s", errtxt);
    }

    /* Get minimum and maximum value */
    maskmin = cpl_vector_get_min(vector);
    maskmax = cpl_vector_get_max(vector);

    /* Identify good and bad mask values */
    if ((maskmin == 0. || maskmin == 1.) && maskmax == 1.) {
        /* Integer values: 0 = bad, 1 = good */
        maskval[0] = maskmin;
        maskval[1] = maskmax;
    } else if (maskmin == 0. && maskmax != 1.) {
        /* Real values: !0 = bad, 0 = good */
        maskval[0] = maskmax;
        maskval[1] = maskmin;
    } else {
        sprintf(errtxt, "%s: cpl_vector *vector (minimum != 0 and != 1)",
                SC_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s",
                                     errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_basic_getmaskval_image(double maskval[2],
                                         const cpl_image *image)
{
    /*!
     * Gets meaning of mask values provided by a CPL image. It is assumed that
     * a good value equals 1 if only 1 and possiby 0 are present. If only 0
     * exists, the definition is reversed. This is also true for a minimum
     * mask value of 0 and any larger number. The latter will be classified
     * as bad mask values.
     *
     * \b INPUT:
     * \param image    input CPL image
     *
     * \b OUTPUT:
     * \param maskval  mask values (1st: bad value, 2nd: good value)
     *
     * \b ERRORS:
     * - No data
     * - Invalid object value(s)
     */

    char errtxt[SC_MAXLEN] = "";
    int nx = 0, ny = 0;
    double maskmin = 0., maskmax = 0.;

    /* Get size of image */
    nx = cpl_image_get_size_x(image);
    ny = cpl_image_get_size_y(image);
    if (nx <= 0 || ny <= 0) {
        sprintf(errtxt, "%s: cpl_image *image", SC_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_NDA, "%s", errtxt);
    }

    /* Get minimum and maximum value */
    maskmin = cpl_image_get_min(image);
    maskmax = cpl_image_get_max(image);

    /* Identify good and bad mask values */
    if ((maskmin == 0. || maskmin == 1.) && maskmax == 1.) {
        /* Integer values: 0 = bad, 1 = good */
        maskval[0] = 0.;
        maskval[1] = 1.;
    } else if (maskmin == 0. && maskmax != 1.) {
        /* Real values: !0 = bad, 0 = good */
        maskval[0] = maskmax;
        maskval[1] = 0.;
    } else {
        sprintf(errtxt, "%s: cpl_image *image (minimum != 0 and != 1)",
                SC_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s",
                                     errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_basic_calcsinc(cpl_vector *sinc)
{
    /*!
     * Calculates asymmetric, exponentially-damped sinc function for rebinning
     * without smoothing. The radius of the kernel, the sampling, and the sinc
     * damping constant are defined by ::SC_SINCRAD_PRECOMP, ::SC_SINCNBIN, and
     * ::SC_SINCDAMP, respectively. All constants are given in pixels. The
     * resulting kernel should be normalised to 1.
     *
     * \note The method is similar to the approach used in the IDL routine
     * "sshift2d.pro".
     *
     * \b INPUT:
     * \param sinc  empty CPL vector
     *
     * \b OUTPUT:
     * \param sinc  CPL vector with damped sinc function
     *
     * \b ERRORS:
     * - CPL_ERROR_ILLEGAL_INPUT
     */

    int nsinc, nsinch, k;
    double dx;
    double *kern;

    /* Get number of data points from header constants */
    nsinc = (2 * SC_SINCRAD_PRECOMP) * SC_SINCNBIN + 1;
    cpl_error_ensure((nsinc % 2 == 1), CPL_ERROR_ILLEGAL_INPUT,
                     return CPL_ERROR_ILLEGAL_INPUT,
                     "sc_basic_calcsinc requires odd number of points");

    /* Initialise vector for function values */
    cpl_vector_set_size(sinc, nsinc);
    cpl_vector_fill(sinc, 0.);

    /* Bisect number of data points */
    nsinch = nsinc / 2;

    /* Set radius and step size of function */
    dx = (2. * SC_SINCRAD_PRECOMP) / (nsinc - 1.);

    /* Get pointer to CPL array */
    kern = cpl_vector_get_data(sinc);

    /* Calculate left-hand side of damped sinc function */
    for (k = 0; k < nsinch; k++) {
        double x = (k - nsinch) * dx;
        if (ceil(x) == x)
            kern[k] = 0.;
        else
            kern[k] = exp(- pow(x / SC_SINCDAMP, 2))
                      * sin(CPL_MATH_PI * x) / (CPL_MATH_PI * x);
    }

    /* Set central value */
    kern[nsinch] = 1.;

    /* Mirror left-hand half to get right-hand function values */
    for (k = nsinch + 1; k < nsinc; k++) {
        kern[k] = kern[nsinc - k - 1];
    }

    return CPL_ERROR_NONE;
}

/**@}*/
