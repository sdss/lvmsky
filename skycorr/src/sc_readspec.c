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
 * \file sc_readspec.c
 *
 * Routines for reading observing data
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  16 Feb 2011
 * \date   17 Sep 2013
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <sc_readspec.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

cpl_error_code sc_readspec(cpl_table *spec, cpl_parameterlist *parlist)
{
    /*!
     * \callgraph
     *
     * Reads a science or sky spectrum (type taken from the general parameter
     * list) in ASCII, FITS table, or FITS image format and puts the read data
     * in a CPL table consisting of the columns "lambda", "flux", and
     * "weight". Only 1D data are supported. The required observing parameters
     * are read from the FITS header (only sky spectrum) or the parameter
     * list.
     *
     * \b INPUT:
     * \param spec      empty CPL table
     * \param parlist   general CPL parameter list
     *
     * \b OUTPUT:
     * \param spec      CPL table with wavelengths, fluxes, and weights
     * \param parlist   general CPL parameter list with FITS header values
     *
     * \b ERRORS:
     * - Invalid object value(s)
     * - Error in subroutine
     * - No data
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameter *p;
    sctarr tabdat;
    char spectype[4], basedir[SC_MAXLEN], outdir[SC_MAXLEN];
    char outname[SC_MAXLEN], outfilename[SC_MAXLEN];
    char errtxt[SC_MAXLEN];
    double meanlam = 0.;

    /* Create output table columns */
    cpl_table_set_size(spec, 0);
    cpl_table_new_column(spec, "lambda", CPL_TYPE_DOUBLE);
    cpl_table_new_column(spec, "flux", CPL_TYPE_DOUBLE);
    cpl_table_new_column(spec, "dflux", CPL_TYPE_DOUBLE);
    cpl_table_new_column(spec, "mask", CPL_TYPE_INT);
    cpl_table_new_column(spec, "weight", CPL_TYPE_DOUBLE);

    /* Get type of spectrum (science or sky) */
    p = cpl_parameterlist_find(parlist, "spectype");
    strncpy(spectype, cpl_parameter_get_string(p), 4);

    /* Get name components of converted file */
    p = cpl_parameterlist_find(parlist, "inst_dir");
    strncpy(basedir, cpl_parameter_get_string(p), SC_MAXLEN);
    p = cpl_parameterlist_find(parlist, "output_dir");
    sc_basic_abspath(outdir, cpl_parameter_get_string(p), basedir);
    p = cpl_parameterlist_find(parlist, "output_name");
    strncpy(outname, cpl_parameter_get_string(p), SC_MAXLEN);

    /* Compose name of converted file ('sci' or 'sky' label) */
    if (strncmp(spectype, "SCI", 3) == 0) {
        sprintf(outfilename, "%s%s_sci.fits", outdir, outname);
    } else if (strncmp(spectype, "SKY", 3) == 0) {
        sprintf(outfilename, "%s%s_sky.fits", outdir, outname);
    } else {
        sprintf(errtxt, "%s: cpl_parameterlist *parlist (spectype neither "
                "'SCI' nor 'SKY')", SC_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }

    /* Read file and convert it to CPL table and CPL property list */
    if ((status = sc_conv_readfile(&tabdat, parlist)) != CPL_ERROR_NONE) {
        if (strncmp(spectype, "SCI", 3) == 0) {
            sprintf(errtxt, "%s: error while processing input science "
                    "spectrum", SC_ERROR_EIS_TXT);
            cpl_error_set_message(cpl_func, SC_ERROR_EIS, "%s", errtxt);
        } else if (strncmp(spectype, "SKY", 3) == 0) {
            sprintf(errtxt, "%s: error while processing input sky spectrum",
                    SC_ERROR_EIS_TXT);
            cpl_error_set_message(cpl_func, SC_ERROR_EIS, "%s", errtxt);
        }
        return status;
    }

    /* Write CPL table and CPL property list to FITS table */
    sc_conv_writetable(&tabdat, parlist);

    /* Free allocated memory */
    sc_conv_tarr_delete(&tabdat);

    /* Get spectroscopic data */
    if (cpl_error_get_code() == CPL_ERROR_NONE) {
        sc_readspec_fits(spec, parlist, outfilename);
    }

    /* Get observing parameters (for reference sky spectrum only) */
    if (cpl_error_get_code() == CPL_ERROR_NONE &&
        strncmp(spectype, "SKY", 3) == 0) {
        sc_readspec_header(parlist, outfilename);
    }

    /* Exit in case of errors */
    if ((status = cpl_error_get_code()) != CPL_ERROR_NONE) {
        if (strncmp(spectype, "SCI", 3) == 0) {
            sprintf(errtxt, "%s: error while processing input science "
                    "spectrum", SC_ERROR_EIS_TXT);
            cpl_error_set_message(cpl_func, SC_ERROR_EIS, "%s", errtxt);
        } else if (strncmp(spectype, "SKY", 3) == 0) {
            sprintf(errtxt, "%s: error while processing input sky spectrum",
                    SC_ERROR_EIS_TXT);
            cpl_error_set_message(cpl_func, SC_ERROR_EIS, "%s", errtxt);
        }
        return status;
    }

    /* Set weights for spectra without error column */
    sc_readspec_setweight(spec, parlist);

    if (cpl_table_get_nrow(spec) == 0) {
        sprintf(errtxt, "%s: cpl_table *spec", SC_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_NDA, "%s", errtxt);
    }

    /* All weights = 0? */
    if (cpl_table_get_column_max(spec, "weight") == 0) {
        sprintf(errtxt, "%s: cpl_table *spec (all weights = 0)",
                SC_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }

    /* Set mean wavelength */
    meanlam = 0.5 * (cpl_table_get(spec, "lambda", 0, NULL) +
                     cpl_table_get(spec, "lambda", cpl_table_get_nrow(spec)-1,
                                   NULL));
    p = cpl_parameterlist_find(parlist, "meanlam");
    cpl_parameter_set_double(p, meanlam);


    return CPL_ERROR_NONE;
}


cpl_error_code sc_readspec_fits(cpl_table *spec,
                                const cpl_parameterlist *parlist,
                                const char *filename)
{
    /*!
     * Reads a FITS file with tabulated spectroscopic data (wavelength, flux,
     * flux error, mask) in the 1st extension (no other extensions are
     * allowed) and puts the read data in a CPL table consisting of the
     * columns "lambda", "flux", and "weight". The names of the required FITS
     * table columns are provided by the general CPL parameter list. The
     * presence of flux error is optional. For skipping such a column, the
     * name has to be 'NONE'. The original input file could also lack a mask
     * column (also indicated by 'NONE'). However, ::sc_conv_readfile makes
     * sure that a suitable mask column indicated by the given column name or
     * ::SC_DEFMASKCOL (in the case of 'NONE') + '_I' is present. A weight of
     * zero is taken if a wavelength is excluded by the mask. If there is no
     * error column, the parameter list default error times mean flux is
     * assumed for all data points.
     *
     * \b INPUT:
     * \param spec      empty CPL table
     * \param parlist   general CPL parameter list
     * \param filename  path and name of FITS file
     *
     * \b OUTPUT:
     * \param spec      CPL table with wavelengths, fluxes, and weights
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     */

    const cpl_parameter *p;
    cpl_table *intab;
    cpl_array *colnames;
    cpl_boolean exerr = CPL_TRUE, exmask = CPL_TRUE;
    char errtxt[SC_MAXLEN], col_lam[SC_LENLINE+1];
    char col_flux[SC_LENLINE+1], col_dflux[SC_LENLINE+1];
    char col_mask[SC_LENLINE+1], col_imask[SC_LENLINE+1];
    char colname[SC_LENLINE+1];
    int coln[4] = {0, 0, 0, 0};
    int next = 0, ncolmin = 5, ncol = 0, check = 0, j = 0, nrow = 0, i = 0;
    int mask = 0;
    double wlgtomicron = 0., lam = 0., flux = 0., dflux = 0.;

    /* Write info message */
    cpl_msg_info(cpl_func, "Read %s", filename);

    /* Check file existence and number of FITS extensions */
    next = cpl_fits_count_extensions(filename);
    if (next < 0) {
        sprintf(errtxt, "%s: %s", SC_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, SC_ERROR_FOF, "%s", errtxt);
    } else if (next != 1) {
        sprintf(errtxt, "%s: %s (# of FITS extensions != 1)",
                SC_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s", errtxt);
    }

    /* Get column labels from parameter list */
    p = cpl_parameterlist_find_const(parlist, "col_lam");
    strncpy(col_lam, cpl_parameter_get_string(p), SC_LENLINE+1);
    if (strcmp(col_lam, "NONE") == 0) {
        sprintf(col_lam, "%s", SC_DEFLAMCOL);
    }
    p = cpl_parameterlist_find_const(parlist, "col_flux");
    strncpy(col_flux, cpl_parameter_get_string(p), SC_LENLINE+1);
    if (strcmp(col_flux, "NONE") == 0) {
        sprintf(col_flux, "%s", SC_DEFFLUXCOL);
    }
    p = cpl_parameterlist_find_const(parlist, "col_dflux");
    strncpy(col_dflux, cpl_parameter_get_string(p), SC_LENLINE+1);
    if (strcmp(col_dflux, "NONE") == 0) {
        cpl_table_erase_column(spec, "dflux");
        exerr = CPL_FALSE;
        ncolmin--;
    }
    p = cpl_parameterlist_find_const(parlist, "col_mask");
    strncpy(col_mask, cpl_parameter_get_string(p), SC_LENLINE+1);
    if (strcmp(col_mask, "NONE") == 0) {
        cpl_table_erase_column(spec, "mask");
        exmask = CPL_FALSE;
        sprintf(col_imask, "%s_I", SC_DEFMASKCOL);
        ncolmin--;
    } else {
        sprintf(col_imask, "%s_I", col_mask);
    }

    /* Get wavelength unit conversion factor */
    p = cpl_parameterlist_find_const(parlist, "wlgtomicron");
    wlgtomicron = cpl_parameter_get_double(p);

    /* Read FITS table */
    intab = cpl_table_load(filename, 1, 0);

    /* Get column labels */
    colnames = cpl_table_get_column_names(intab);
    ncol = cpl_array_get_size(colnames);

    /* Check existence of columns */

    for (check = 0, j = 0; j < ncol; j++) {
        strncpy(colname, cpl_array_get_string(colnames, j),
                SC_LENLINE+1);
        if (strcmp(colname, col_lam) == 0) {
            coln[0] = j;
            check++;
        } else if (strcmp(colname, col_flux) == 0) {
            coln[1] = j;
            check++;
        } else if (strcmp(colname, col_dflux) == 0) {
            coln[2] = j;
            check++;
        } else if (strcmp(colname, col_mask) == 0) {
            check++;
        } else if (strcmp(colname, col_imask) == 0) {
            coln[3] = j;
            check++;
       }
    }

    if (check < ncolmin) {
        cpl_table_set_size(spec, 0);
        cpl_array_delete(colnames);
        cpl_table_delete(intab);
        sprintf(errtxt, "%s: %s (missing FITS table column(s))",
                SC_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s",
                                     errtxt);
    }

    /* Resize output table */
    nrow = cpl_table_get_nrow(intab);
    cpl_table_set_size(spec, nrow);

    /* Transfer wavelength and flux/transmission and compute weights */

    for (i = 0; i < nrow; i++) {

        lam = cpl_table_get(intab,
                            cpl_array_get_string(colnames, coln[0]), i, NULL)
              * wlgtomicron;
        cpl_table_set(spec, "lambda", i, lam);

        flux = cpl_table_get(intab,
                             cpl_array_get_string(colnames, coln[1]), i,
                             NULL);
        cpl_table_set(spec, "flux", i, flux);

        if (exerr == CPL_TRUE) {
            dflux = cpl_table_get(intab,
                                  cpl_array_get_string(colnames, coln[2]), i,
                                  NULL);
            cpl_table_set(spec, "dflux", i, dflux);
        } else {
            dflux = 1.;
        }

        mask = cpl_table_get(intab, cpl_array_get_string(colnames, coln[3]),
                             i, NULL);
        if (exmask == CPL_TRUE) {
            cpl_table_set(spec, "mask", i, mask);
        }

        if (dflux <= 0. || mask == 0) {
            cpl_table_set(spec, "weight", i, 0.);
        } else {
            cpl_table_set(spec, "weight", i, 1. / dflux);
        }

    }

    /* Delete temporary CPL objects */
    cpl_array_delete(colnames);
    cpl_table_delete(intab);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_readspec_header(cpl_parameterlist *parlist,
                                  const char *filename)
{
    /*!
     * Reads keywords from a FITS file header and puts the read values in the
     * general CPL parameter list. The presence of keywords for date (MJD or
     * date in years), time (UTC in s), and telescope altitude angle (in deg)
     * is required. The keyword names are included in the general parameter
     * list. The default names are MJD_OBS, TM-START, and ESO TEL ALT. If an
     * expected keyword cannot be found in the FITS header, an error message
     * will be written to the CPL error structure. The values of the required
     * parameters are also part of the general parameter list. If a value is
     * set manually, this value is taken instead of the corresponding FITS
     * keyword content.
     *
     * \b INPUT:
     * \param parlist   general CPL parameter list
     * \param filename  path and name of FITS file
     *
     * \b OUTPUT:
     * \param parlist   general CPL parameter list with FITS header values
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     */

    cpl_parameter *p;
    cpl_propertylist *header = NULL;
    cpl_property *prop;
    cpl_type type;
    char key[SC_NKEY][SC_LENLINE+1] = {"", "", ""}, errtxt[SC_MAXLEN];
    int next = 0, ext = 0, i = 0, nerr = 0;
    double val = 0.;

    /* Parameter names related to FITS keywords */
    char keypar[SC_NKEY][SC_LENLINE+1] =
         {"date_key", "time_key", "telalt_key"};
    char valpar[SC_NKEY][SC_LENLINE+1] =
         {"date_val", "time_val", "telalt_val"};

    /* Find FITS extension with required keywords and write the header of this
       extension into CPL property list */
    next = cpl_fits_count_extensions(filename);
    for (ext = 0; ext <= next; ext++) {
        header = cpl_propertylist_load(filename, ext);
        if (header == NULL) {
            sprintf(errtxt, "%s: %s", SC_ERROR_FOF_TXT, filename);
            return cpl_error_set_message(cpl_func, SC_ERROR_FOF, "%s",
                                         errtxt);
        }
        /* Check first keyword */
        p = cpl_parameterlist_find(parlist, keypar[0]);
        strncpy(key[0], cpl_parameter_get_string(p), SC_LENLINE+1);
        prop = cpl_propertylist_get_property(header, key[0]);
        if (prop != NULL) {
            break;
        } else {
            if (ext < next) {
                cpl_propertylist_delete(header);
            }
        }
    }

    /* Write info message */
    cpl_msg_info(cpl_func, "Take keywords from FITS extension %d", ext);

    /* Update general parameter list with keywords from CPL property list */

    for (i = 0; i < SC_NKEY; i++) {

        /* Get value from parameter list */
        p = cpl_parameterlist_find(parlist, valpar[i]);
        val = cpl_parameter_get_double(p);

        /* Get name of required FITS keyword */
        p = cpl_parameterlist_find(parlist, keypar[i]);
        strncpy(key[i], cpl_parameter_get_string(p), SC_LENLINE+1);

        /* Write info message */
        cpl_msg_info(cpl_func, "Read keyword %s", key[i]);

        /* Get FITS keyword */
        prop = cpl_propertylist_get_property(header, key[i]);

        if (prop == NULL && val == -1.) {

            /* Set error message in the case of missing keyword */
            nerr++;
            sprintf(errtxt, "%s: %s (keyword %s not found)",
                    SC_ERROR_UFS_TXT, filename, key[i]);
            cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s", errtxt);

        } else {

            /* Take FITS keyword only if parameter value was not given
               manually */
            if (val == -1.) {
                /* Check for type */
                type = cpl_property_get_type(prop);
                if (type == CPL_TYPE_DOUBLE)  {
                    val = cpl_property_get_double(prop);
                } else if (type == CPL_TYPE_FLOAT) {
                    val = (double) cpl_property_get_float(prop);
                } else if (type == CPL_TYPE_INT) {
                    val = (double) cpl_property_get_int(prop);
                } else {
                    nerr++;
                    sprintf(errtxt, "%s: %s (non-numerical keyword %s)",
                            SC_ERROR_UFS_TXT, filename, key[i]);
                    cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s",
                                          errtxt);
                    continue;
                }
            }

            /* MJD -> date in years (if required; save MJD before) */
            if (i == 0 && val > 3000.) {
                p = cpl_parameterlist_find(parlist, "mjd");
                cpl_parameter_set_double(p, val);
                val = sc_basic_mjd2fracyear(val);
            }

            /* Write new value into parameter list */
            p = cpl_parameterlist_find(parlist, valpar[i]);
            cpl_parameter_set_double(p, val);

        }

    }

    /* Free allocated memory */
    cpl_propertylist_delete(header);

    /* Return SC_ERROR_UFS in the case of keyword mismatch */
    if (nerr > 0) {
        return SC_ERROR_UFS;
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_readspec_setweight(cpl_table *spec,
                                     const cpl_parameterlist *parlist)
{
    /*!
     * Sets the weights for spectra that do not contain an error column.
     * The default relative error provided by the general CPL parameter list
     * is multiplied by the mean flux of all wavelengths that are used by the
     * fit. The resulting absolute error is taken for all wavelengths with
     * non-zero weight.
     *
     * \b INPUT:
     * \param spec     CPL table with observed spectrum
     * \param parlist  general CPL parameter list
     *
     * \b OUTPUT:
     * \param spec     spectrum with optimised weights
     *
     * \b ERRORS:
     * - Invalid object value(s)
     */

    cpl_error_code err = CPL_ERROR_NONE;
    const cpl_parameter *p;
    char col_dflux[SC_LENLINE+1], errtxt[SC_MAXLEN];
    int nrow = 0, i = 0, nw = 0;
    double deferr = 0., weight = 0., fsum = 0., weight0 = 0.;

    /* Return if error data exist */
    p = cpl_parameterlist_find_const(parlist, "col_dflux");
    strncpy(col_dflux, cpl_parameter_get_string(p), SC_LENLINE+1);
    if (strcmp(col_dflux, "NONE") != 0) {
        return CPL_ERROR_NONE;
    }

    /* Get default relative error */
    p = cpl_parameterlist_find_const(parlist, "default_error");
    deferr = cpl_parameter_get_double(p);

    /* Count data points with non-zero weight and sum up their fluxes */
    nrow = cpl_table_get_nrow(spec);
    for (i = 0; i < nrow; i++) {
        weight = cpl_table_get(spec, "weight", i, NULL);
        if (weight > 0) {
            nw++;
            fsum += cpl_table_get(spec, "flux", i, NULL);
        }
    }

    /* Set errors to default error * mean if error column is missing */
    if (fsum > 0) {
        weight0 = (double) nw / (fsum * deferr);
        for (i = 0; i < nrow; i++) {
            weight = cpl_table_get(spec, "weight", i, NULL);
            if (weight > 0) {
                cpl_table_set(spec, "weight", i, weight0);
            }
        }
    }

    /* Check for errors */

    if (fsum == 0) {
        sprintf(errtxt, "%s: cpl_table *spec (all fluxes = 0)",
                SC_ERROR_IOV_TXT);
        err = cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }

    if (nw == 0) {
        sprintf(errtxt, "%s: cpl_table *spec (all weights = 0)",
                SC_ERROR_IOV_TXT);
        err = cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }

    return err;
}

/**@}*/
