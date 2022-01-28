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
 * \file sc_fwhmest.c
 *
 * \author Marco Barden, Stefan Noll, & ESO In-Kind Team Innsbruck
 *
 * \date   28 Aug 2013
 *
 */


/*****************************************************************************
 *                                  INCLUDES                                 *
 ****************************************************************************/

#include <sc_fwhmest.h>


/*****************************************************************************
 *                                    CODE                                   *
 ****************************************************************************/

cpl_error_code sc_fwhmest(double *fwhm, double *rms,
                          cpl_table *spec, cpl_table *lines,
                          const cpl_parameterlist *parlist)
{
    /*!
     * \brief
     *   Estimate FWHM using isolated lines in spectrum.
     *
     * This routine estimates the FWHM of the lines in a spectrum based on a
     * set of identified isolated lines. Isolated lines are marked in the
     * line list. Firstly, the FWHM is approximated by the value in the line
     * list. Afterwards, fitting the data with a Gaussian is attempted. If the
     * fitting is successful, this value is preferred over the initial
     * estimate.
     *
     * The input spectrum has to be a CPL table containing (at least) the
     * following columns:
     *
     *      lambda
     *      lflux
     *      class
     *
     * The column "class" is expected to contain the following values:
     *      0 = Continuum
     *      1 = Line pixel
     *      2 = Line peak
     *      3 = Isolated line peak
     *
     * \note In the case of a wavelength-dependent line width (varfwhm = 1),
     *       the output FWHM is provided for the central wavelength of the
     *       input spectrum.
     *
     * \b INPUT:
     * \param spec     spectrum
     * \param lines    line list
     * \param parlist  general CPL parameter list
     *
     * \b OUTPUT:
     * \param fwhm     FWHM
     * \param rms      RMS error for the FWHM estimate
     *
     * \b ERRORS:
     * - CPL_ERROR_NONE:  no error occurred
     * - SC_ERROR_NDA:    no isolated lines in spectra
     */

    cpl_error_code err_code = CPL_ERROR_NONE;

    cpl_array *fwhm_arr = cpl_array_new(0, CPL_TYPE_DOUBLE),
              *iso_arr  = cpl_array_new(0, CPL_TYPE_DOUBLE);

    const cpl_parameter *p;

    int status;                         // MPFIT status
    int npar = 3;                       // Number of parameters for Gauss fit
    cpl_size maxpos;
    double pars[3];                     // User parameters

    long i, j,                          // loop variables
         n_lines,                       // number of lines in line list
         i_start = 0,                   // start of line
         i_end = 0,                     // end of line
         nline = 0,                     // number of pixels in current line
         niso = 0;                      // number of isolated lines

    int varfwhm = 0;                    // wavelength-dependent FWHM?
    double meanlam = 0;                 // mean wavelength
    double total = 0;                   // total flux in current line
    double width = 0;                   // line width

    scxydata data;                      // data for Gauss fit

    char msg[SC_MAXLEN];

    /* Variable FWHM (linear change with wavelength)? */
    p = cpl_parameterlist_find_const(parlist, "varfwhm");
    varfwhm = cpl_parameter_get_int(p);

    /* Get mean wavelength */
    p = cpl_parameterlist_find_const(parlist, "meanlam");
    meanlam = cpl_parameter_get_double(p);

    /* Initialise data array */
    data.x = cpl_array_new(0, CPL_TYPE_DOUBLE);
    data.y = cpl_array_new(0, CPL_TYPE_DOUBLE);
    data.y_err = cpl_array_new(0, CPL_TYPE_DOUBLE);

    /* Get number of lines in line list */
    n_lines = cpl_table_get_nrow(lines);

    /* Loop over all lines */

    for (i = 0; i < n_lines; i++) {

        /* Use isolated lines for FWHM estimate using Gauss fit */
        if (cpl_table_get_int(lines, "isol_flag", i, NULL) == 1) {

            i_start = cpl_table_get_int(lines, "line_px_start", i, NULL);
            i_end = cpl_table_get_int(lines, "line_px_end", i, NULL);

            /* Number of line pixels */
            nline = i_end - i_start;

            /* Copy data to xydata structure */
            cpl_array_set_size(data.x, nline);
            cpl_array_set_size(data.y_err, nline);
            cpl_array_set_size(data.y, nline);

            /* Fill private data structure to pointers with user data */
            cpl_array_copy_data_double(data.y,
                        cpl_table_get_data_double(spec, "lflux")+i_start);
            for (j = 0, total = 0; j < nline; j++) {
                cpl_array_set_double(data.x, j, (double) (i_start + j));
                cpl_array_set_double(data.y_err, j, 1.);
                total += cpl_array_get_double(data.y, j, NULL);
            }

            /* Get initial FWHM from line list */
            width = cpl_table_get_double(lines, "fwhm", i, NULL);
            if (width == 0.) {
                /* First iteration */
                width = (double) cpl_table_get_int(lines, "width", i, NULL);
            } else {
                /* Correct FWHM if line width is wavelength dependent */
                if (varfwhm == 1) {
                    width *= cpl_table_get_double(lines, "peak_lam", i,
                                                  NULL) / meanlam;
                }
            }

            /* Call mpfit() */
            pars[0] = total;
            pars[1] = cpl_table_get_int(lines, "peak_loc", i, NULL);
            pars[2] = width / CPL_MATH_FWHM_SIG;
            status = mpfit(sc_fwhmest_gaussfunc, nline, npar, pars, 0, 0,
                           (void *) &data, 0);

            sc_basic_status2txt(msg, status);
            //printf("status: %s, npar: %i, p0: %lf, p1: %lf, p2: %lf, "
            //       "FWHM: %lf\n", msg, npar, pars[0], pars[1],
            //       pars[2]*CPL_MATH_FWHM_SIG, fwhm);

            /* If fit was successful, overwrite result from crude estimate */
            if (status == MP_OK_CHI ||
                status == MP_OK_PAR ||
                status == MP_OK_BOTH ||
                status == MP_OK_DIR) {
                width = pars[2] * CPL_MATH_FWHM_SIG;
            }

            /* Count isolated lines */
            niso++;

            /* Correct FWHM if line width is wavelength dependent */
            if (varfwhm == 1) {
                width *= meanlam / cpl_table_get_double(lines, "peak_lam", i,
                                                        NULL);
            }

            /* Set FWHM in line list */
            cpl_table_set_double(lines, "fwhm", i, width);

            /* Store FWHM values in CPL array */
            cpl_array_set_size(fwhm_arr, niso);
            cpl_array_set_double(fwhm_arr, niso-1, width);

            /* Store FWHM value in array for isolated lines */
            cpl_array_set_size(iso_arr, niso);
            cpl_array_set_double(iso_arr, niso-1, width);

        }

    }

    /* Calculate mean FWHM and RMS depending on number of isolated lines */

    switch (niso) {
    case 0:
        *fwhm = 0;
        *rms = HUGE_VAL;
        break;
    case 1:
        *fwhm = cpl_array_get_double(iso_arr, 0, NULL);
        *rms = HUGE_VAL;
        break;
    case 2:
        *fwhm = cpl_array_get_min(iso_arr);
        *rms = HUGE_VAL;
        break;
    case 3:
        *fwhm = cpl_array_get_median(iso_arr);
        *rms = HUGE_VAL;
        break;
    case 4:
        cpl_array_get_maxpos(iso_arr, &maxpos);
        cpl_array_set_invalid(iso_arr, maxpos);
        *fwhm = cpl_array_get_median(iso_arr);
        *rms = HUGE_VAL;
        break;
    default:
        /* Calculate mean/RMS of FWHM array */
        sc_basic_clipmean(fwhm, rms, fwhm_arr, CPL_FALSE);
        break;
    }

    /* Clean-up */
    cpl_array_delete(fwhm_arr);
    cpl_array_delete(iso_arr);
    cpl_array_delete(data.x);
    cpl_array_delete(data.y);
    cpl_array_delete(data.y_err);

    return err_code;
}


int sc_fwhmest_gaussfunc(int n_data, int n_pars, double *p, double *deviates,
                         double **derivs, void *pdata)
{
    /*!
     * User function for CMPFIT. Provides error-weighted differences between a
     * model Gaussian and a measured line profile.
     *
     * \b INPUT:
     * \param n_data    number of data points
     * \param n_pars    number of parameters
     * \param p         array of fit parameters
     * \param derivs    derivatives (not used)
     * \param pdata     private data
     *
     * \b OUTPUT:
     * \param deviates  array of residuals
     *
     * \b ERRORS:
     * - none
     */

    int i;

    scxydata *data = (scxydata *) pdata;

    double *fgauss;

    if (n_pars) {};
    if (derivs) {};

    fgauss = (double *) malloc(n_data * sizeof(double));
    sc_basic_gaussfunc(fgauss, cpl_array_get_data_double(data->x), n_data, p);

    /* Compute function deviates */
    for (i = 0; i < n_data; i++) {
        deviates[i] = (cpl_array_get_double(data->y, i, NULL) - fgauss[i]) /
                       cpl_array_get_double(data->y_err, i, NULL);
    }

    free(fgauss);

    return 0;
}

/**@}*/
