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
 * \file sc_conv.c
 *
 * Routines for conversion of files
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  16 Apr 2013
 * \date   31 Jan 2014
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <sc_conv.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

cpl_error_code sc_conv_preptable(const char *parfile_in)
{
    /*!
     * \callgraph
     *
     * Converts input file into a FITS table for SKYCORR.
     *
     * \b INPUT:
     * \param parfile_in  name of parameter file
     *
     * \b ERRORS:
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameterlist *parlist;
    sctarr tabdat;
    char parfile[SC_MAXLEN] = "";

    /* Read driver file */
    sc_basic_absfile(parfile, parfile_in);
    parlist = cpl_parameterlist_new();
    if ((status = sc_par_readfile(parlist, parfile)) != CPL_ERROR_NONE) {
        cpl_parameterlist_delete(parlist);
        return status;
    }

    /* Read file and convert it to CPL table and CPL property list */
    if ((status = sc_conv_readfile(&tabdat, parlist)) !=
        CPL_ERROR_NONE) {
        cpl_parameterlist_delete(parlist);
        return status;
    }

    /* Write CPL table and CPL property list to FITS table */
    sc_conv_writetable(&tabdat, parlist);

    /* Free allocated memory */
    cpl_parameterlist_delete(parlist);
    sc_conv_tarr_delete(&tabdat);

    /* Return error code of last error */
    return cpl_error_get_code();
}


cpl_error_code sc_conv_readfile(sctarr *tabdat,
                                const cpl_parameterlist *parlist)
{
    /*!
     * \callgraph
     *
     * Reads an ASCII file or FITS file (either FITS table or 1D fits image)
     * and write its/their data into a ::sctarr structure, which contains an
     * array of CPL tables (spectra) and CPL property lists (header keywords).
     *
     * \b INPUT:
     * \param parlist  general CPL parameter list
     *
     * \b OUTPUT:
     * \param tabdat   ::sctarr structure with data of input file
     *
     * \b ERRORS:
     * - Invalid object values
     * - File opening failed
     * - Unexpected file structure
     * - see subroutines
     */

    const cpl_parameter *p;
    cpl_array *colnames = NULL;
    cpl_table *extnames = NULL;
    scvarr vecdat;
    char spectype[4] = "", filename[SC_MAXLEN] = "";
    char basedir[SC_MAXLEN] = "";
    int fitsformat = 0;

    /* Get type of input spectrum (science or sky) */
    p = cpl_parameterlist_find_const(parlist, "spectype");
    strncpy(spectype, cpl_parameter_get_string(p), 4);

    /* Take science (spectype = "SCI") or sky spectrum (spectype = "SKY")
       and write info message */
    if (strncmp(spectype, "SCI", 3) == 0) {
        p = cpl_parameterlist_find_const(parlist, "scispec");
        strncpy(filename, cpl_parameter_get_string(p), SC_MAXLEN);
        cpl_msg_info(cpl_func, "Input science data file: %s", filename);
    } else if (strncmp(spectype, "SKY", 3) == 0) {
        p = cpl_parameterlist_find_const(parlist, "skyspec");
        strncpy(filename, cpl_parameter_get_string(p), SC_MAXLEN);
        cpl_msg_info(cpl_func, "Input sky data file: %s", filename);
    } else {
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV,
                         "%s: cpl_parameterlist *parlist (spectype neither "
                         "'SCI' nor 'SKY')", SC_ERROR_IOV_TXT);
    }

    /* Modify path of input file if not absolute path */
    p = cpl_parameterlist_find_const(parlist, "inst_dir");
    strncpy(basedir, cpl_parameter_get_string(p), SC_MAXLEN);
    if (filename[0] != '/') {
        char relfilename[SC_MAXLEN];
        if(snprintf(relfilename, SC_MAXLEN, "%s/%s", basedir, filename) >=
           SC_MAXLEN) {
            return cpl_error_set_message(cpl_func, SC_ERROR_NAMETOOLONG,
                                         "Relative filename too long: %s/%s",
                                         basedir, filename);
        }
        sc_basic_absfile(filename, relfilename);
    }

    /* Get file type */
    sc_conv_checkfitsformat(&fitsformat, filename);

    /* Read file dependent on type */
    if (fitsformat == 0) {
        /* Read non-FITS file (ASCII format assumed) */
        colnames = cpl_array_new(0, CPL_TYPE_STRING);
        sc_conv_setcolnames(colnames, parlist);
        sc_conv_ascii_read(tabdat, filename, colnames);
        cpl_array_delete(colnames);
    } else if (fitsformat == 1) {
        /* Read FITS table */
        sc_conv_tarr_read(tabdat, filename);
    } else if (fitsformat == 2) {
        /* Read 1D FITS image */
        sc_conv_varr_read(&vecdat, filename);
        extnames = cpl_table_new(0);
        if (cpl_error_get_code() == CPL_ERROR_NONE) {
            sc_conv_setextnames_varr(extnames, &vecdat, parlist);
        }
        if (cpl_error_get_code() == CPL_ERROR_NONE) {
            sc_conv_varr2tarr(tabdat, &vecdat, extnames);
        }
        sc_conv_varr_delete(&vecdat);
        cpl_table_delete(extnames);
    } else if (fitsformat > 2) {
        return cpl_error_set_message(cpl_func, SC_ERROR_UFS,
                                     "%s: %s (multidimensional FITS image)",
                                     SC_ERROR_UFS_TXT, filename);
    }

    /* Adapt data for use by SKYCORR */
    if (cpl_error_get_code() == CPL_ERROR_NONE) {
        sc_conv_modtable(tabdat, parlist);
    }

    /* Return error code of last error */
    return cpl_error_get_code();
}


cpl_error_code sc_conv_checkfitsformat(int *fitsformat, const char *filename)
{
    /*!
     * Identifies the FITS type (table or image) by means of the header
     * keyword XTENSION of the first FITS extension. If an extension is not
     * present, the image format is assumed. The keyword NAXIS is used to
     * distinguish between 1D, 2D, and 3D images. The routine returns 0 for
     * non-FITS format (e.g. ASCII), 1 for FITS table, 2 for 1D FITS image, 3
     * for 2D FITS image, 4 for 3D FITS image, and -1 in the case of errors.
     * It is assumed that all extensions have the same FITS type.
     *
     * \b INPUT:
     * \param filename    path and name of input FITS file
     *
     * \b OUTPUT:
     * \param fitsformat  flag for FITS format
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     */

    FILE *stream;
    cpl_propertylist *header = NULL;
    cpl_property *prop;
    int next = 0;

    /* Check file existence */
    if ((stream = fopen(filename, "r")) == NULL) {
        *fitsformat = -1;
        return cpl_error_set_message(cpl_func, SC_ERROR_FOF, "%s: %s",
                                     SC_ERROR_FOF_TXT, filename);
    }
    fclose(stream);

    /* Get number of extensions */
    next = cpl_fits_count_extensions(filename);

    /* FITS table or FITS image? */

    if (next == -1) {
        /* Not a FITS file */
        cpl_errorstate_set(CPL_ERROR_NONE);
        *fitsformat = 0;
    } else if (next == 0) {
        /* File header from zeroth extension */
        header = cpl_propertylist_load(filename, 0);
        /* FITS tables require at least one extension */
        *fitsformat = 2;
    } else {
        /* File header from first extension */
        header = cpl_propertylist_load(filename, 1);
        /* Read FITS keyword XTENSION */
        prop = cpl_propertylist_get_property(header, "XTENSION");
        if (prop == NULL) {
            *fitsformat = -1;
            cpl_propertylist_delete(header);
            return cpl_error_set_message(cpl_func, SC_ERROR_UFS,
                                         "%s: %s (keyword XTENSION not found)",
                                         SC_ERROR_UFS_TXT, filename);
        } else {
            if (strcmp(cpl_property_get_string(prop), "IMAGE") == 0) {
                /* FITS image */
                *fitsformat = 2;
            } else if (strcmp(cpl_property_get_string(prop), "BINTABLE")
                       == 0) {
                /* FITS table */
                *fitsformat = 1;
            } else {
                *fitsformat = -1;
                cpl_propertylist_delete(header);
                return cpl_error_set_message(cpl_func, SC_ERROR_UFS,
                         "%s: %s (keyword XTENSION != IMAGE or BINTABLE)",
                         SC_ERROR_UFS_TXT, filename);
            }
        }
    }

    /* 1D or 2D FITS image? */

    if (*fitsformat == 2) {
        prop = cpl_propertylist_get_property(header, "NAXIS");
        if (prop == NULL) {
            *fitsformat = -1;
            cpl_propertylist_delete(header);
            return cpl_error_set_message(cpl_func, SC_ERROR_UFS,
                                         "%s: %s (keyword NAXIS not found)",
                                         SC_ERROR_UFS_TXT, filename);
        } else {
            if (cpl_property_get_int(prop) == 1) {
                /* 1D FITS image */
                *fitsformat = 2;
            } else if (cpl_property_get_int(prop) == 2) {
                /* 2D FITS image */
                *fitsformat = 3;
            } else if (cpl_property_get_int(prop) == 3) {
                /* 3D FITS image */
                *fitsformat = 4;
            } else {
                *fitsformat = -1;
                cpl_propertylist_delete(header);
                return cpl_error_set_message(cpl_func, SC_ERROR_UFS,
                                     "%s: %s (keyword NAXIS != 1, 2, or 3)",
                                     SC_ERROR_UFS_TXT, filename);
            }
        }

    }

    /* Free allocated memory */
    cpl_propertylist_delete(header);

    /* Return error code of last error */
    return cpl_error_get_code();
}


cpl_error_code sc_conv_setcolnames(cpl_array *colnames,
                                   const cpl_parameterlist *parlist)
{
    /*!
     * Gets column names for ASCII file from the general parameter list. The
     * mandatory columns are wavelength and flux. Optionally, error and mask
     * can also be provided. If this is not the case, this is indicated by
     * 'NULL' in the parameter structure.
     *
     * \b INPUT:
     * \param colnames  empty CPL array for strings
     * \param parlist   general CPL parameter list
     *
     * \b OUTPUT:
     * \param colnames  CPL array of column names for reading of ASCII file
     *
     * \b ERRORS:
     * - none
     */

    const cpl_parameter *p;
    char col_lam[SC_LENLINE+1] = "", col_flux[SC_LENLINE+1] = "";
    char col_dflux[SC_LENLINE+1] = "", col_mask[SC_LENLINE+1] = "";
    int j = 0;

    /* Set default size of array for column names */
    cpl_array_set_size(colnames, 4);

    /* Get names of wavelength and flux column */
    p = cpl_parameterlist_find_const(parlist, "col_lam");
    strncpy(col_lam, cpl_parameter_get_string(p), SC_LENLINE+1);
    cpl_array_set_string(colnames, j, col_lam);
    j++;
    p = cpl_parameterlist_find_const(parlist, "col_flux");
    strncpy(col_flux, cpl_parameter_get_string(p), SC_LENLINE+1);
    cpl_array_set_string(colnames, j, col_flux);
    j++;

    /* Get information on existence and names of error and mask column */
    p = cpl_parameterlist_find_const(parlist, "col_dflux");
    strncpy(col_dflux, cpl_parameter_get_string(p), SC_LENLINE+1);
    if (strcmp(col_dflux, "NONE") != 0) {
        cpl_array_set_string(colnames, j, col_dflux);
        j++;
    }
    p = cpl_parameterlist_find_const(parlist, "col_mask");
    strncpy(col_mask, cpl_parameter_get_string(p), SC_LENLINE+1);
    if (strcmp(col_mask, "NONE") != 0) {
        cpl_array_set_string(colnames, j, col_mask);
        j++;
    }

    /* Resize array */
    cpl_array_set_size(colnames, j);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_setextnames_varr(cpl_table *extnames,
                                        const scvarr *vecdat,
                                        const cpl_parameterlist *parlist)
{
    /*!
     * Prepares CPL table with FITS extension numbers and names of a ::scvarr
     * structure. Extension names are provided by the \e columns parameter of
     * the driver file. The name "NONE" will be replaced by the default label
     * for the data type.
     *
     * \b INPUT:
     * \param extnames  empty CPL table
     * \param vecdat    ::scvarr structure with data of 1D FITS image
     * \param parlist   general CPL parameter list
     *
     * \b OUTPUT:
     * \param extnames  CPL table with FITS extension numbers and names
     *
     * \b ERRORS:
     * - Invalid object structure
     * - see subroutines
     */

    const cpl_parameter *p;
    cpl_property *prop;
    cpl_boolean none_lam = CPL_FALSE, none_flux = CPL_FALSE;
    char errtxt[SC_MAXLEN] = "", colname[SC_LENLINE+1] = "";
    char extname[SC_LENLINE+1] = "";
    char **col;
    int ncol0 = 4, ncol = 0, next = 0, check = 0, j = 0, h = 0;
    int *extn;

    /* Set default size of output table */
    cpl_table_set_size(extnames, ncol0);

    /* Create columns for extension numbers and names */
    if (cpl_table_has_column(extnames, "col") != 1) {
        cpl_table_new_column(extnames, "col", CPL_TYPE_STRING);
    }
    if (cpl_table_has_column(extnames, "extn") != 1) {
        cpl_table_new_column(extnames, "extn", CPL_TYPE_INT);
    }
    cpl_table_fill_column_window_int(extnames, "extn", 0, ncol0, -1);

    /* Get column/extension names from parameter list */

    ncol = ncol0;

    p = cpl_parameterlist_find_const(parlist, "col_lam");
    sprintf(colname, "%s", cpl_parameter_get_string(p));
    cpl_table_set_string(extnames, "col", 0, colname);
    if (strcmp(colname, "NONE") == 0) {
        cpl_table_set_string(extnames, "col", 0, SC_DEFLAMCOL);
        ncol--;
        none_lam = CPL_TRUE;
    }

    p = cpl_parameterlist_find_const(parlist, "col_flux");
    sprintf(colname, "%s", cpl_parameter_get_string(p));
    cpl_table_set_string(extnames, "col", 1, colname);
    if (strcmp(colname, "NONE") == 0) {
        cpl_table_set_string(extnames, "col", 1, SC_DEFFLUXCOL);
        none_flux = CPL_TRUE;
    }

    p = cpl_parameterlist_find_const(parlist, "col_dflux");
    sprintf(colname, "%s", cpl_parameter_get_string(p));
    cpl_table_set_string(extnames, "col", 2, colname);
    if (strcmp(colname, "NONE") == 0) {
        ncol--;
    }

    p = cpl_parameterlist_find_const(parlist, "col_mask");
    sprintf(colname, "%s", cpl_parameter_get_string(p));
    cpl_table_set_string(extnames, "col", 3, colname);
    if (strcmp(colname, "NONE") == 0) {
        ncol--;
    }

    /* Get pointers to table columns */
    col = cpl_table_get_data_string(extnames, "col");
    extn = cpl_table_get_data_int(extnames, "extn");

    /* Get number of extensions */
    next = vecdat->next;

    /* Check existence of extensions with expected names and set extension
       numbers */

    for (check = 0, j = 0; j <= next; j++) {
        prop = cpl_propertylist_get_property(vecdat->head[j], "EXTNAME");
        if (prop == NULL && j > 0) {
            continue;
        } else if (prop == NULL && j == 0) {
            /* Make sure that flux vector is counted independent of existence
               of keyword EXTNAME (requires 'NONE' as extension name)  */
            if (none_flux == CPL_TRUE) {
                extn[1] = 0;
                check++;
            }
        } else {
            strncpy(extname, cpl_property_get_string(prop), SC_LENLINE+1);
            for (h = 0; h < ncol0; h++) {
                if (strncmp(extname, col[h], SC_LENLINE+1) == 0) {
                    extn[h] = j;
                    check++;
                }
            }
        }
    }

    if (check != ncol) {
        if (extn[0] < 0 && none_lam == CPL_FALSE) {
            sprintf(errtxt, "%s: scvarr *vecdat (wavelength extension '%s' "
                    "not found)", SC_ERROR_IOS_TXT, col[0]);
            return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s",
                                         errtxt);
        }
        if (extn[1] < 0 && none_flux == CPL_FALSE) {
            sprintf(errtxt, "%s: scvarr *vecdat (flux extension '%s' not "
                    "found)", SC_ERROR_IOS_TXT, col[1]);
            return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s",
                                         errtxt);
        }
        if (extn[2] < 0 && strcmp(col[2], "NONE") != 0) {
            sprintf(errtxt, "%s: scvarr *vecdat (flux error extension '%s' "
                    "not found)", SC_ERROR_IOS_TXT, col[2]);
            return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s",
                                         errtxt);
        }
        if (extn[3] < 0 && strcmp(col[3], "NONE") != 0) {
            sprintf(errtxt, "%s: scvarr *vecdat (mask extension '%s' not "
                    "found)", SC_ERROR_IOS_TXT, col[3]);
            return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s",
                                         errtxt);
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_setextnames_iarr(cpl_table *extnames,
                                        const sciarr *imadat,
                                        const cpl_parameterlist *parlist)
{
    /*!
     * Prepares CPL table with FITS extension numbers and names of a ::sciarr
     * structure. Extension names are provided by the \e columns parameter of
     * the driver file. The name "NONE" will be replaced by the default label
     * for the data type.
     *
     * \b INPUT:
     * \param extnames  empty CPL table
     * \param imadat    ::sciarr structure with data of 1D FITS image
     * \param parlist   general CPL parameter list
     *
     * \b OUTPUT:
     * \param extnames  CPL table with FITS extension numbers and names
     *
     * \b ERRORS:
     * - Invalid object structure
     * - see subroutines
     */

    const cpl_parameter *p;
    cpl_property *prop;
    cpl_boolean none_lam = CPL_FALSE, none_flux = CPL_FALSE;
    char errtxt[SC_MAXLEN] = "", colname[SC_LENLINE+1] = "";
    char extname[SC_LENLINE+1] = "";
    char **col;
    int ncol0 = 4, ncol = 0, next = 0, check = 0, j = 0, h = 0;
    int *extn;

    /* Set default size of output table */
    cpl_table_set_size(extnames, ncol0);

    /* Create columns for extension numbers and names */
    if (cpl_table_has_column(extnames, "col") != 1) {
        cpl_table_new_column(extnames, "col", CPL_TYPE_STRING);
    }
    if (cpl_table_has_column(extnames, "extn") != 1) {
        cpl_table_new_column(extnames, "extn", CPL_TYPE_INT);
    }
    cpl_table_fill_column_window_int(extnames, "extn", 0, ncol0, -1);

    /* Get column/extension names from parameter list */

    ncol = ncol0;

    p = cpl_parameterlist_find_const(parlist, "col_lam");
    sprintf(colname, "%s", cpl_parameter_get_string(p));
    cpl_table_set_string(extnames, "col", 0, colname);
    if (strcmp(colname, "NONE") == 0) {
        cpl_table_set_string(extnames, "col", 0, SC_DEFLAMCOL);
        ncol--;
        none_lam = CPL_TRUE;
    }

    p = cpl_parameterlist_find_const(parlist, "col_flux");
    sprintf(colname, "%s", cpl_parameter_get_string(p));
    cpl_table_set_string(extnames, "col", 1, colname);
    if (strcmp(colname, "NONE") == 0) {
        cpl_table_set_string(extnames, "col", 1, SC_DEFFLUXCOL);
        none_flux = CPL_TRUE;
    }

    p = cpl_parameterlist_find_const(parlist, "col_dflux");
    sprintf(colname, "%s", cpl_parameter_get_string(p));
    cpl_table_set_string(extnames, "col", 2, colname);
    if (strcmp(colname, "NONE") == 0) {
        ncol--;
    }

    p = cpl_parameterlist_find_const(parlist, "col_mask");
    sprintf(colname, "%s", cpl_parameter_get_string(p));
    cpl_table_set_string(extnames, "col", 3, colname);
    if (strcmp(colname, "NONE") == 0) {
        ncol--;
    }

    /* Get pointers to table columns */
    col = cpl_table_get_data_string(extnames, "col");
    extn = cpl_table_get_data_int(extnames, "extn");

    /* Get number of extensions */
    next = imadat->next;

    /* Check existence of extensions with expected names and set extension
       numbers */

    for (check = 0, j = 0; j <= next; j++) {
        prop = cpl_propertylist_get_property(imadat->head[j], "EXTNAME");
        if (prop == NULL && j > 0) {
            continue;
        } else if (prop == NULL && j == 0) {
            /* Make sure that flux vector is counted independent of existence
               of keyword EXTNAME (requires 'NONE' as extension name)  */
            if (none_flux == CPL_TRUE) {
                extn[1] = 0;
                check++;
            }
        } else {
            strncpy(extname, cpl_property_get_string(prop), SC_LENLINE+1);
            for (h = 0; h < ncol0; h++) {
                if (strncmp(extname, col[h], SC_LENLINE+1) == 0) {
                    extn[h] = j;
                    check++;
                }
            }
        }
    }

    if (check != ncol) {
        if (extn[0] < 0 && none_lam == CPL_FALSE) {
            sprintf(errtxt, "%s: sciarr *imadat (wavelength extension '%s' "
                    "not found)", SC_ERROR_IOS_TXT, col[0]);
            return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s",
                                         errtxt);
        }
        if (extn[1] < 0 && none_flux == CPL_FALSE) {
            sprintf(errtxt, "%s: sciarr *imadat (flux extension '%s' not "
                    "found)", SC_ERROR_IOS_TXT, col[1]);
            return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s",
                                         errtxt);
        }
        if (extn[2] < 0 && strcmp(col[2], "NONE") != 0) {
            sprintf(errtxt, "%s: sciarr *imadat (flux error extension '%s' "
                    "not found)", SC_ERROR_IOS_TXT, col[2]);
            return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s",
                                         errtxt);
        }
        if (extn[3] < 0 && strcmp(col[3], "NONE") != 0) {
            sprintf(errtxt, "%s: sciarr *imadat (mask extension '%s' not "
                    "found)", SC_ERROR_IOS_TXT, col[3]);
            return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s",
                                         errtxt);
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_varr2tarr(sctarr *tabdat, const scvarr *vecdat,
                                 const cpl_table *extnames)
{
    /*!
     * Fills CPL tables and CPL property lists of a ::sctarr structure by data
     * read from a 1D FITS image with optional extensions for error, mask, and
     * wavelength. The image data is provided by a ::scvarr structure. The
     * extension numbers and names are given by a CPL table. Extensions are
     * not read if names are missing or wrong. Images for different chips have
     * to be converted separately.
     *
     * \b INPUT:
     * \param vecdat    ::scvarr structure with data of 1D FITS image
     * \param extnames  CPL table with FITS extension numbers and names
     *
     * \b OUTPUT:
     * \param tabdat    ::sctarr structure with data of FITS table
     *
     * \b ERRORS:
     * - No data
     * - Invalid object structure
     * - Invalid object value(s)
     */

    char errtxt[SC_MAXLEN] = "";
    const char **col;
    int ncol = 0, next = 0, h = 0, nrow = 0, j = 0, i = 0;
    const int *extn;
    double crpix = 0., crval = 0., cdelt = 0.;

    /* Initialise CPL tables and CPL property lists for content of scvarr
       structure (memory allocation) */
    sc_conv_tarr_init(tabdat, 1);
    tabdat->tab[0] = cpl_table_new(0);
    tabdat->tab[1] = cpl_table_new(0);

    /* Transfer FITS header data */
    tabdat->head[0] = cpl_propertylist_duplicate(vecdat->head[0]);
    tabdat->head[1] = cpl_propertylist_new();

    /* Check table for extension numbers and names */
    ncol = cpl_table_get_nrow(extnames);
    if (ncol <= 0) {
        sc_conv_tarr_delete(tabdat);
        sprintf(errtxt, "%s: cpl_table *extnames", SC_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_NDA, "%s", errtxt);
    }

    /* Check column names */
    if (cpl_table_has_column(extnames, "col") != 1 ||
        cpl_table_has_column(extnames, "extn") != 1) {
        sc_conv_tarr_delete(tabdat);
        sprintf(errtxt, "%s: cpl_table *extnames (column names)",
                SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Check existence of data */
    if (cpl_table_get_column_max(extnames, "extn") == -1) {
        sc_conv_tarr_delete(tabdat);
        sprintf(errtxt, "%s: cpl_table *extnames (no extension with valid "
                "data)", SC_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }

    /* Get number of extensions */
    next = vecdat->next;

    /* Check extension numbers */
    if (cpl_table_get_column_max(extnames, "extn") > next) {
        sc_conv_tarr_delete(tabdat);
        sprintf(errtxt, "%s: cpl_table *extnames (incorrect extension "
                "number)", SC_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }

    /* Get pointers to table columns */
    col = cpl_table_get_data_string_const(extnames, "col");
    extn = cpl_table_get_data_int_const(extnames, "extn");

    /* Create table columns */
    cpl_table_new_column(tabdat->tab[1], col[0], CPL_TYPE_DOUBLE);
    for (h = 1; h < ncol; h++) {
        if (extn[h] >= 0) {
            cpl_table_new_column(tabdat->tab[1], col[h], CPL_TYPE_DOUBLE);
        }
    }

    /* Set size of data table */
    nrow = cpl_vector_get_size(vecdat->vec[0]);
    if (nrow <= 0) {
        sc_conv_tarr_delete(tabdat);
        sprintf(errtxt, "%s: scvarr *vecdat", SC_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_NDA, "%s", errtxt);
    } else if (nrow > 0 && next > 0) {
        for (j = 1; j <= next; j++) {
            if (cpl_vector_get_size(vecdat->vec[j]) != nrow) {
                sc_conv_tarr_delete(tabdat);
                sprintf(errtxt, "%s: scvarr *vecdat (vector size differs for "
                        "different extensions", SC_ERROR_IOS_TXT);
                return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s",
                                             errtxt);
            }
        }
    }
    cpl_table_set_size(tabdat->tab[1], nrow);

    /* Get wavelength grid from FITS header if a wavelength vector is not
       provided */
    if (extn[0] < 0) {
         crpix = sc_conv_getwcskey(vecdat->head[0], "CRPIX1");
         if (cpl_error_get_code() == CPL_ERROR_DATA_NOT_FOUND) {
             sc_conv_tarr_delete(tabdat);
             sprintf(errtxt, "%s: scvarr *vecdat (keyword CRPIX1 not found)",
                     SC_ERROR_IOS_TXT);
             return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s",
                                          errtxt);
         }
         crval = sc_conv_getwcskey(vecdat->head[0], "CRVAL1");
         if (cpl_error_get_code() == CPL_ERROR_DATA_NOT_FOUND) {
             sc_conv_tarr_delete(tabdat);
             sprintf(errtxt, "%s: scvarr *vecdat (keyword CRVAL1 not found)",
                     SC_ERROR_IOS_TXT);
             return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s",
                                          errtxt);
         }
         cdelt = sc_conv_getwcskey(vecdat->head[0], "CDELT1");
         if (cpl_error_get_code() == CPL_ERROR_DATA_NOT_FOUND) {
             sc_conv_tarr_delete(tabdat);
             sprintf(errtxt, "%s: scvarr *vecdat (keyword CDELT1 not found)",
                     SC_ERROR_IOS_TXT);
             return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s",
                                          errtxt);
         }
    }

    /* Convert vectors to table columns */
    for (h = 0; h < ncol; h++) {
        for (i = 0; i < nrow; i++) {
            if (extn[h] >= 0) {
                cpl_table_set(tabdat->tab[1], col[h], i,
                              cpl_vector_get(vecdat->vec[extn[h]], i));
            } else if (extn[0] < 0 && h == 0) {
                cpl_table_set(tabdat->tab[1], col[h], i,
                              crval + (i - crpix + 1) * cdelt);
            }
        }
    }

    return CPL_ERROR_NONE;
}


double sc_conv_getwcskey(const cpl_propertylist *plist, const char *key)
{
    /*!
     * Retrieves a FITS WCS key from a property list
     *
     * \b INPUT:
     * \param plist  cpl_propertylist
     * \param key    key to retrieve
     *
     * \b RETURN:
     * - FITS WCS key
     *
     * \b ERRORS:
     * - CPL_ERROR_DATA_NOT_FOUND: key does not exist
     */

    double d;
    const cpl_errorstate cleanstate = cpl_errorstate_get();

    d = cpl_propertylist_get_double(plist, key);

    /* WCS keys can be written as ints but should be interpreted as doubles */
    if (cpl_error_get_code() == CPL_ERROR_TYPE_MISMATCH) {
        cpl_errorstate_set(cleanstate);
        d = (double) cpl_propertylist_get_int(plist, key);
    }

    return d;
}


static cpl_table * sc_convert_sdp_table(const cpl_table * inptable,
                                const char * col_lam, const char * col_flux,
                                const char * col_dflux, const char * col_mask)
{
    /*!
     * Convert SDP formatted VOCLASS Spectrum data to internal format.
     * This format consists of one row containing arrays with the data in the
     * columns
     *
     * \b INPUT:
     * \param inptable   input Spectrum table with one row containing arrays
     * \param col_lam    name of wavelength column
     * \param col_flux   name of flux column
     * \param col_dflux  name of flux error column or NULL
     * \param col_mask   name of mask column or NULL
     *
     * \b OUTPUT:
     * new table in with arrays expanded to rows
     *
     * \b ERRORS:
     * - Missing or inconsistent columns
     */
    const cpl_array * alam = cpl_table_get_array(inptable, col_lam, 0);
    if (alam == NULL) {
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, "Expected %s "
                              "array column in input", col_lam);;
        return NULL;
    }
    const cpl_array * aflux = cpl_table_get_array(inptable, col_flux, 0);
    if (aflux == NULL) {
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, "Expected %s "
                              "array column in input", col_flux);
        return NULL;
    }
    const cpl_array * adflux = NULL, * amask = NULL;
    if (col_dflux) {
        adflux = cpl_table_get_array(inptable, col_dflux, 0);
        if (adflux == NULL) {
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, "Expected %s "
                                  "array column in input", col_dflux);
            return NULL;
        }
    }
    if (col_mask) {
        amask = cpl_table_get_array(inptable, col_mask, 0);
        if (adflux == NULL) {
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, "Expected %s "
                                  "array column in input", col_mask);
            return NULL;
        }
    }
    if (cpl_array_get_size(alam) != cpl_array_get_size(aflux) ||
        (adflux && cpl_array_get_size(alam) != cpl_array_get_size(adflux)) ||
        (amask && cpl_array_get_size(alam) != cpl_array_get_size(amask))) {
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Table columns do not all have the same size");
        return NULL;
    }

    cpl_table * ntab = cpl_table_new(cpl_array_get_size(alam));

    /* const casts due to PIPE-6075 */
    cpl_array * dalam = cpl_array_cast((cpl_array*)alam, CPL_TYPE_DOUBLE);
    cpl_table_new_column(ntab, col_lam, CPL_TYPE_DOUBLE);
    cpl_table_copy_data_double(ntab, col_lam, cpl_array_get_data_double_const(dalam));
    cpl_array_delete(dalam);

    cpl_array * daflux = cpl_array_cast((cpl_array*)aflux, CPL_TYPE_DOUBLE);
    cpl_table_new_column(ntab, col_flux, CPL_TYPE_DOUBLE);
    cpl_table_copy_data_double(ntab, col_flux, cpl_array_get_data_double_const(daflux));
    cpl_array_delete(daflux);

    if(adflux) {
        cpl_array * dadflux = cpl_array_cast((cpl_array*)adflux, CPL_TYPE_DOUBLE);
        cpl_table_new_column(ntab, col_dflux, CPL_TYPE_DOUBLE);
        cpl_table_copy_data_double(ntab, col_dflux, cpl_array_get_data_double_const(dadflux));
        cpl_array_delete(dadflux);
    }
    if (amask) {
        cpl_table_new_column(ntab, col_mask, CPL_TYPE_INT);
        cpl_table_copy_data_int(ntab, col_mask, cpl_array_get_data_int_const(amask));
    }

    return ntab;
}


cpl_error_code sc_conv_modtable(sctarr *tabdat,
                                const cpl_parameterlist *parlist)
{
    /*!
     * Modifies data in ::sctarr structure in order to be consistent with
     * requirements of SKYCORR. Only the first table extension is considered.
     * The routine creates a column with the integer mask values 0 (rejected)
     * and 1 (ok). By default, it is expected that input mask values agree
     * with this definition. If other values are found, 0 is assumed to be ok
     * and all other values cause pixel rejection. Possible nan in the flux
     * and error columns are substituted by zero flux and indicated by
     * mask = 0. The same is performed for negative errors. Finally,
     * unreliable fluxes at the edges (defined by ::SC_EDGEFRAC) are masked.
     *
     * \b INPUT:
     * \param tabdat   ::sctarr structure with read data
     * \param parlist  general CPL parameter list
     *
     * \b OUTPUT:
     * \param tabdat   ::sctarr structure with modified data
     *
     * \b ERRORS:
     * - No data
     * - Invalid object structure
     */

    const cpl_parameter *p;
    char col_lam[SC_LENLINE+1] = "", col_flux[SC_LENLINE+1] = "";
    char col_dflux[SC_LENLINE+1] = "", col_mask[SC_LENLINE+1] = "";
    char errtxt[SC_MAXLEN] = "", col_imask[SC_LENLINE+1] = "";
    char col_renamed[SC_LENLINE+1] = "";
    cpl_boolean exerr = CPL_TRUE, exmask = CPL_TRUE, isnanflux = CPL_FALSE;
    cpl_boolean isnanum = CPL_FALSE, isedgeflux = CPL_FALSE;
    cpl_boolean isedge = CPL_FALSE, isnandflux = CPL_FALSE;
    cpl_boolean isoutrange = CPL_FALSE, is0 = CPL_FALSE, isnomask = CPL_FALSE;
    cpl_boolean ismask0 = CPL_FALSE;
    int next = 0, nrow = 0, edgepix = 0, i = 0, nmask0 = 0, imask = 0;
    double minflux = 0., maxflux = 0., flux = 0., dflux = 0., mask = 0.;

    /* Get names of wavelength and flux column */
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

    /* Get information on existence and names of error and mask column */
    p = cpl_parameterlist_find_const(parlist, "col_dflux");
    strncpy(col_dflux, cpl_parameter_get_string(p), SC_LENLINE+1);
    if (strcmp(col_dflux, "NONE") == 0) {
        exerr = CPL_FALSE;
    }
    p = cpl_parameterlist_find_const(parlist, "col_mask");
    strncpy(col_mask, cpl_parameter_get_string(p), SC_LENLINE+1);
    if (strcmp(col_mask, "NONE") == 0) {
        exmask = CPL_FALSE;
    }

    /* Check number of extensions */
    next = tabdat->next;
    if (next <= 0) {
        sprintf(errtxt, "%s: sctarr *tabdat", SC_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_NDA, "%s", errtxt);
    } else if (next > 1) {
        sprintf(errtxt, "%s: sctarr *tabdat (number of extensions != 1)",
                SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Check existence of expected columns */
    if (cpl_table_has_column(tabdat->tab[1], col_lam) != 1) {
        sprintf(errtxt, "%s: sctarr *tabdat (wavelength column '%s' not "
                "found)", SC_ERROR_IOS_TXT, col_lam);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }
    if (cpl_table_has_column(tabdat->tab[1], col_flux) != 1) {
        sprintf(errtxt, "%s: sctarr *tabdat (flux column '%s' not found)",
                SC_ERROR_IOS_TXT, col_flux);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }
    if (exerr == TRUE) {
        if (cpl_table_has_column(tabdat->tab[1], col_dflux) != 1) {
            sprintf(errtxt, "%s: sctarr *tabdat (flux error column '%s' not "
                    "found)", SC_ERROR_IOS_TXT, col_dflux);
            return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s",
                                         errtxt);
        }
    }
    if (exmask == TRUE) {
        if (cpl_table_has_column(tabdat->tab[1], col_mask) != 1) {
            sprintf(errtxt, "%s: sctarr *tabdat (mask column '%s' not found)",
                    SC_ERROR_IOS_TXT, col_mask);
            return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s",
                                         errtxt);
        }
    }

    /* check for and convert sdp data format */
    if (cpl_propertylist_has(tabdat->head[1], "VOCLASS") &&
        strncmp(cpl_propertylist_get_string(tabdat->head[1],
                                            "VOCLASS"),
                "SPECTRUM", strlen("SPECTRUM")) == 0) {
        cpl_table * ntab = sc_convert_sdp_table(tabdat->tab[1], col_lam, col_flux,
                                                exerr ? col_dflux : NULL,
                                                exmask ? col_mask : NULL);
        if (cpl_error_get_code()) {
            return cpl_error_get_code();
        }
        cpl_table_delete(tabdat->tab[1]);
        tabdat->tab[1] = ntab;
    }

    /* Create integer mask column (and rename already existing column of the
       same name) */
    if (exmask == CPL_TRUE) {
        sprintf(col_imask, "%s_I", col_mask);
    } else {
        sprintf(col_imask, "%s_I", SC_DEFMASKCOL);
    }
    if (cpl_table_has_column(tabdat->tab[1], col_imask) == 1) {
        sprintf(col_renamed, "%s_orig", col_imask);
        if (cpl_table_has_column(tabdat->tab[1], col_renamed) != 1) {
            cpl_msg_info(cpl_func, "Name of internal integer mask already "
                         "used: Rename %s in %s", col_imask, col_renamed);
            cpl_table_name_column(tabdat->tab[1], col_imask, col_renamed);
        } else {
            cpl_msg_info(cpl_func, "Use of reserved mask column names: "
                         "Erase %s, keep %s", col_imask, col_renamed);
            cpl_table_erase_column(tabdat->tab[1], col_imask);
        }
    }
    cpl_table_new_column(tabdat->tab[1], col_imask, CPL_TYPE_INT);

    /* Get number of rows */
    nrow = cpl_table_get_nrow(tabdat->tab[1]);

    /* Get number of edge pixels */
    edgepix = (int) ceil(nrow * SC_EDGEFRAC);

    /* Get minimum and maximum flux excluding the edge pixels */
    for (minflux = 0., maxflux = 0., i = edgepix; i < nrow - edgepix; i++) {
        flux = cpl_table_get(tabdat->tab[1], col_flux, i, NULL);
        if (flux < minflux) {
            minflux = flux;
        } else if (flux > maxflux) {
            maxflux = flux;
        }
    }

    /* Check for nan values and negative errors and correct them */

    for (nmask0 = 0, i = 0; i < nrow; i++) {

        flux = cpl_table_get(tabdat->tab[1], col_flux, i, NULL);
        if (isnan(flux) != 0) {
            cpl_table_set(tabdat->tab[1], col_flux, i, 0.);
            isnanflux = CPL_TRUE;
            isnanum = CPL_TRUE;
            isedgeflux = CPL_FALSE;
        } else if (flux < SC_MAXFLUXFAC * minflux ||
                   flux > SC_MAXFLUXFAC * maxflux) {
            isnanflux = CPL_FALSE;
            isedgeflux = CPL_TRUE;
            isedge = CPL_TRUE;
        } else {
            isnanflux = CPL_FALSE;
            isedgeflux = CPL_FALSE;
        }

        if (exerr == CPL_TRUE) {
            dflux = cpl_table_get(tabdat->tab[1], col_dflux, i, NULL);
            if (dflux <= 0 || isnan(dflux) != 0) {
                cpl_table_set(tabdat->tab[1], col_dflux, i, 0.);
                isnandflux = CPL_TRUE;
                isoutrange = CPL_TRUE;
            } else {
                isnandflux = CPL_FALSE;
            }
        }

        if (exmask == CPL_TRUE) {
            mask = cpl_table_get(tabdat->tab[1], col_mask, i, NULL);
        } else {
            mask = 1;
        }
        if (isnanflux == CPL_TRUE || isnandflux == CPL_TRUE ||
            isedgeflux == CPL_TRUE) {
            cpl_table_set(tabdat->tab[1], col_imask, i, -1);
            if (mask != 0 && mask != 1) {
                isnomask = CPL_TRUE;
            } else if (mask == 0) {
                is0 = CPL_TRUE;
                nmask0++;
            }
        } else {
            if (mask == 0) {
                cpl_table_set(tabdat->tab[1], col_imask, i, 0);
                is0 = CPL_TRUE;
                nmask0++;
            } else {
                cpl_table_set(tabdat->tab[1], col_imask, i, 1);
                if (mask != 1) {
                    isnomask = CPL_TRUE;
                }
            }
        }

    }

    /* Return if mask cannot be interpreted */
    if (isnomask == CPL_TRUE && is0 == CPL_FALSE) {
        sprintf(errtxt, "%s: sctarr *tabdat (all mask value(s) != 0 or 1)",
                SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Correct integer mask values if required (e.g. reverse definition) */

    for (i = 0; i < nrow; i++) {
        imask = cpl_table_get(tabdat->tab[1], col_imask, i, NULL);
        if (isnomask == CPL_TRUE) {
            if (imask == 0) {
                cpl_table_set(tabdat->tab[1], col_imask, i, 1);
            } else {
                cpl_table_set(tabdat->tab[1], col_imask, i, 0);
            }
        } else {
            if (imask == -1) {
                cpl_table_set(tabdat->tab[1], col_imask, i, 0);
            } else if (nmask0 == nrow) {
               cpl_table_set(tabdat->tab[1], col_imask, i, 1);
               ismask0 = CPL_TRUE;
            }
        }
    }

    /* Print info message in the case of bad fluxes, errors, or mask
       values */

    if (isnanum == CPL_TRUE) {
        cpl_msg_info(cpl_func, "Input data: flux(es) = 'nan' "
                     "-> set mask = 0");
    }

    if (isedge == CPL_TRUE) {
        cpl_msg_info(cpl_func, "Input data: unreliable edge flux(es) "
                     "-> set mask = 0");
    }

    if (isoutrange == CPL_TRUE) {
        cpl_msg_info(cpl_func, "Input data: error(s) <= 0 or 'nan' "
                     "-> set mask = 0");
    }

    if (isnomask == CPL_TRUE) {
        cpl_msg_info(cpl_func, "Input data: mask value(s) != 0 or 1 "
                     "-> reverse definition (0 -> 1; != 0 -> 0)");
    }

    if (ismask0 == CPL_TRUE) {
        cpl_msg_info(cpl_func, "Input data: all mask values = 0 "
                     "-> reverse definition (0 -> 1)");
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_writetable(const sctarr *tabdat,
                                  const cpl_parameterlist *parlist)
{
    /*!
     * Writes CPL table (spectrum) and CPL property list (header keywords)
     * from an ::sctarr structure into a FITS table. Depending on science or
     * sky spectrum type, the suffix 'sci' or 'sky' is added to the output
     * file name.
     *
     * \b INPUT:
     * \param tabdat   ::sctarr structure with FITS data
     * \param parlist  general CPL parameter list
     *
     * \b ERRORS:
     * - Invalid object value(s)
     */

    const cpl_parameter *p;
    char basedir[SC_MAXLEN] = "", outdir[SC_MAXLEN] = "";
    char outname[SC_MAXLEN] = "", spectype[4] = "";
    char outfile[SC_MAXLEN] = "", errtxt[SC_MAXLEN] = "";

    /* Get output file name and write info message */
    p = cpl_parameterlist_find_const(parlist, "inst_dir");
    strncpy(basedir, cpl_parameter_get_string(p), SC_MAXLEN);
    p = cpl_parameterlist_find_const(parlist, "output_dir");
    sc_basic_abspath(outdir, cpl_parameter_get_string(p), basedir);
    p = cpl_parameterlist_find_const(parlist, "output_name");
    strncpy(outname, cpl_parameter_get_string(p), SC_MAXLEN);
    p = cpl_parameterlist_find_const(parlist, "spectype");
    strncpy(spectype, cpl_parameter_get_string(p), 4);
    if (strncmp(spectype, "SCI", 3) == 0) {
        sprintf(outfile, "%s%s_sci.fits", outdir, outname);
        cpl_msg_info(cpl_func, "Convert input science data file into %s",
                     outfile);
    } else if (strncmp(spectype, "SKY", 3) == 0) {
        sprintf(outfile, "%s%s_sky.fits", outdir, outname);
        cpl_msg_info(cpl_func, "Convert input sky data file into %s",
                     outfile);
    } else {
        sprintf(errtxt, "%s: cpl_parameterlist *parlist (spectype neither "
                "'SCI' nor 'SKY')", SC_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }

    /* Write new FITS table for first chip and new extension for all chips */
    sc_conv_tarr_write(outfile, tabdat);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_writeresults(const char *parfile_in)
{
    /*!
     * \callgraph
     *
     * Writes results of SKYCORR into a file of the same format as the input
     * data file.
     *
     * \b INPUT:
     * \param parfile_in  name of parameter file
     *
     * \b ERRORS:
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_table *results = NULL;
    cpl_parameterlist *parlist;
    char parfile[SC_MAXLEN] = "";

    /* Read driver file */
    sc_basic_absfile(parfile, parfile_in);
    parlist = cpl_parameterlist_new();
    if ((status = sc_par_readfile(parlist, parfile)) != CPL_ERROR_NONE) {
        cpl_parameterlist_delete(parlist);
        return status;
    }

    /* Read file with CALCTRANS results and convert it into CPL table */
    results = cpl_table_new(0);
    if ((status = sc_conv_readresults(results, parlist)) != CPL_ERROR_NONE) {
        cpl_parameterlist_delete(parlist);
        cpl_table_delete(results);
        return status;
    }

    /* Write CPL table and CPL property list to file with format of input data
       file */
    sc_conv_writefile(results, parlist);

    /* Free allocated memory */
    cpl_parameterlist_delete(parlist);
    cpl_table_delete(results);

    /* Return error code of last error */
    return cpl_error_get_code();
}


cpl_error_code sc_conv_readresults(cpl_table *results,
                                   const cpl_parameterlist *parlist)
{
    /*!
     * Fills CPL table with the results of SKYCORR saved as a FITS table.
     * Converts micron into initial wavelength units if necessary.
     *
     * \b INPUT:
     * \param results  empty CPL table
     * \param parlist  general CPL parameter list
     *
     * \b OUTPUT:
     * \param results  CPL table with results of SKYCORR
     *
     * \b ERRORS:
     * - File opening failed
     */

    const cpl_parameter *p;
    cpl_table *tmptab = NULL;
    char basedir[SC_MAXLEN] = "", outdir[SC_MAXLEN] = "";
    char outname[SC_MAXLEN] = "", filename[SC_MAXLEN] = "";
    char errtxt[SC_MAXLEN] = "";
    double wlgtomicron = 0.;

    /* Get path and name of results file */
    p = cpl_parameterlist_find_const(parlist, "inst_dir");
    strncpy(basedir, cpl_parameter_get_string(p), SC_MAXLEN);
    p = cpl_parameterlist_find_const(parlist, "output_dir");
    sc_basic_abspath(outdir, cpl_parameter_get_string(p), basedir);
    p = cpl_parameterlist_find_const(parlist, "output_name");
    strncpy(outname, cpl_parameter_get_string(p), SC_MAXLEN);
    sprintf(filename, "%s%s_fit.fits", outdir, outname);

    /* Load FITS table into CPL table */
    tmptab = cpl_table_load(filename, 1, 0);
    if (tmptab == NULL) {
        sprintf(errtxt, "%s: %s", SC_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, SC_ERROR_FOF, "%s", errtxt);
    }

    /* Copy content of temporary table into results table */
    sc_basic_copytable_full(results, tmptab);

    /* Get wavelength unit conversion factor */
    p = cpl_parameterlist_find_const(parlist, "wlgtomicron");
    wlgtomicron = cpl_parameter_get_double(p);

    /* Change wavelength units if necessary */
    cpl_table_divide_scalar(results, "lambda", wlgtomicron);

    /* Free allocated memory */
    cpl_table_delete(tmptab);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_writefile(cpl_table *results,
                                 const cpl_parameterlist *parlist)
{
    /*!
     * \callgraph
     *
     * Writes results from SKYCORR into an ASCII file or FITS file (FITS
     * table or 1D FITS image) dependent on the format of the input data file.
     * Apart from sky-subtracted fluxes, several additional SKYCORR-related
     * keywords are written.
     *
     * \b INPUT:
     * \param results  CPL table with results of SKYCORR
     * \param parlist  general CPL parameter list
     *
     * \b OUTPUT:
     * \param results  results of SKYCORR with intial wavelength grid
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     * - see subroutines
     */

    const cpl_parameter *p;
    cpl_table *extnames = NULL;
    sctarr tabdat;
    scvarr vecdat;
    char spectype[4] = "", filename[SC_MAXLEN] = "";
    char errtxt[SC_MAXLEN] = "", basedir[SC_MAXLEN] = "";
    char relfilename[SC_MAXLEN] = "", outdir[SC_MAXLEN] = "";
    char indir[SC_MAXLEN] = "", infilename[SC_MAXLEN] = "";
    char suffix[SC_MAXLEN] = "", outfilename[SC_MAXLEN] = "";
    int fitsformat = 0;

    /* Create sctarr structure and read FITS table with data of initial input
       file */
    if (sc_conv_readprepfits(&tabdat, parlist) != CPL_ERROR_NONE) {
        return cpl_error_get_code();
    }

    /* Write results into sctarr structure */
    if (sc_conv_results2tarr(&tabdat, results) != CPL_ERROR_NONE) {
        sc_conv_tarr_delete(&tabdat);
        return cpl_error_get_code();
    }

    /* Get type of input spectrum (science or sky) */
    p = cpl_parameterlist_find_const(parlist, "spectype");
    strncpy(spectype, cpl_parameter_get_string(p), 4);

    /* Take science (spectype = "SCI") or sky spectrum (spectype = "SKY") */
    if (strncmp(spectype, "SCI", 3) == 0) {
        p = cpl_parameterlist_find_const(parlist, "scispec");
        strncpy(filename, cpl_parameter_get_string(p), SC_MAXLEN);
    } else if (strncmp(spectype, "SKY", 3) == 0) {
        p = cpl_parameterlist_find_const(parlist, "skyspec");
        strncpy(filename, cpl_parameter_get_string(p), SC_MAXLEN);
    } else {
        sc_conv_tarr_delete(&tabdat);
        sprintf(errtxt, "%s: cpl_parameterlist *parlist (spectype neither "
                "'SCI' nor 'SKY')", SC_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }

    /* Modify path of input file if not absolute path */
    p = cpl_parameterlist_find_const(parlist, "inst_dir");
    strncpy(basedir, cpl_parameter_get_string(p), SC_MAXLEN);
    if (filename[0] != '/') {
        sprintf(relfilename, "%s/%s", basedir, filename);
        sc_basic_absfile(filename, relfilename);
    }

    /* Get output file name */
    p = cpl_parameterlist_find_const(parlist, "output_dir");
    sc_basic_abspath(outdir, cpl_parameter_get_string(p), basedir);
    sc_basic_getfilename(indir, infilename, suffix, filename);
    sprintf(outfilename, "%s%s_SC.%s", outdir, infilename, suffix);

    /* Write info message */
    cpl_msg_info(cpl_func, "Corrected input science data file: %s",
                 outfilename);

    /* Get file type */
    sc_conv_checkfitsformat(&fitsformat, filename);

    /* Write file dependent on input type */
    if (fitsformat == 0) {
        /* Write non-FITS file (ASCII format assumed) */
        sc_conv_erasemaskcol(&tabdat, parlist);
        sc_conv_ascii_write(outfilename, &tabdat);
    } else if (fitsformat == 1) {
        /* Write FITS table */
        sc_conv_erasemaskcol(&tabdat, parlist);
        sc_par_addkeywords(tabdat.head[0], parlist); // add. FITS keywords
        sc_conv_tarr_write(outfilename, &tabdat);
    } else if (fitsformat == 2) {
        /* Write 1D FITS images */
        sc_conv_varr_read(&vecdat, filename);
        extnames = cpl_table_new(0);
        if (cpl_error_get_code() == CPL_ERROR_NONE) {
            sc_conv_setextnames_varr(extnames, &vecdat, parlist);
        }
        if (cpl_error_get_code() == CPL_ERROR_NONE) {
            sc_conv_resultstarr2varr(&vecdat, extnames, &tabdat);
            sc_par_addkeywords(vecdat.head[0], parlist); // add. FITS keywords
        }
        if (cpl_error_get_code() == CPL_ERROR_NONE) {
            sc_conv_varr_write(outfilename, &vecdat);
        }
        sc_conv_varr_delete(&vecdat);
        cpl_table_delete(extnames);
    } else if (fitsformat > 2) {
        sc_conv_tarr_delete(&tabdat);
        sprintf(errtxt, "%s: %s (multidimensional FITS image)",
                SC_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s", errtxt);
    }

    /* Free allocated memory */
    sc_conv_tarr_delete(&tabdat);

    /* Return error code of last error */
    return cpl_error_get_code();
}


cpl_error_code sc_conv_readprepfits(sctarr *tabdat,
                                    const cpl_parameterlist *parlist)
{
    /*!
     * Fills CPL table of a new ::sctarr structure by data read from a FITS
     * table prepared for SKYCORR.
     *
     * \b INPUT:
     * \param parlist  general CPL parameter list
     *
     * \b OUTPUT:
     * \param tabdat   ::sctarr structure with data of FITS file
     *
     * \b ERRORS:
     * - Invalid object value(s)
     * - see ::sc_conv_tarr_read
     */

    const cpl_parameter *p;
    char basedir[SC_MAXLEN] = "", outdir[SC_MAXLEN] = "";
    char outname[SC_MAXLEN] = "", spectype[4] = "";
    char filename[SC_MAXLEN] = "", errtxt[SC_MAXLEN] = "";

    /* Get path and name of results file */
    p = cpl_parameterlist_find_const(parlist, "inst_dir");
    strncpy(basedir, cpl_parameter_get_string(p), SC_MAXLEN);
    p = cpl_parameterlist_find_const(parlist, "output_dir");
    sc_basic_abspath(outdir, cpl_parameter_get_string(p), basedir);
    p = cpl_parameterlist_find_const(parlist, "output_name");
    strncpy(outname, cpl_parameter_get_string(p), SC_MAXLEN);
    p = cpl_parameterlist_find_const(parlist, "spectype");
    strncpy(spectype, cpl_parameter_get_string(p), 4);
    if (strncmp(spectype, "SCI", 3) == 0) {
        sprintf(filename, "%s%s_sci.fits", outdir, outname);
    } else if (strncmp(spectype, "SKY", 3) == 0) {
        sprintf(filename, "%s%s_sky.fits", outdir, outname);
    } else {
        sprintf(errtxt, "%s: cpl_parameterlist *parlist (spectype neither "
                "'SCI' nor 'SKY')", SC_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }

    /* Read FITS table */
    sc_conv_tarr_read(tabdat, filename);

    /* Return error code of last error */
    return cpl_error_get_code();
}


cpl_error_code sc_conv_results2tarr(sctarr *tabdat, cpl_table *results)
{
    /*!
     * Writes results from SKYCORR into an ::sctarr structure that already
     * contains data from the input file. A new table column "scflux" is
     * created, which contains the data from the "scflux" column of the
     * results data table. If columns named "scdflux" and "scmask" exist in
     * the results data table, they are also copied by using the same names
     * for the output table.
     *
     * \b INPUT:
     * \param tabdat   ::sctarr structure with data of FITS file prepared for
     *                 SKYCORR
     * \param results  CPL table with results of SKYCORR
     *
     * \b OUTPUT:
     * \param tabdat   input ::sctarr structure extended for SKYCORR results
     *
     * \b ERRORS:
     * - No data
     * - Invalid object structure
     * - Inconsistent data grids
     */

    char errtxt[SC_MAXLEN] = "";
    int next = 0, nrow = 0;

    /* Check existence of required column in results data table */
    if (cpl_table_has_column(results, "scflux") != 1) {
        sprintf(errtxt, "%s: cpl_table *results (required column 'scflux' "
                "not present)", SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Check number of extensions */
    next = tabdat->next;
    if (next <= 0) {
        sprintf(errtxt, "%s: sctarr *tabdat", SC_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_NDA, "%s", errtxt);
    } else if (next > 1) {
        sprintf(errtxt, "%s: sctarr *tabdat (number of extensions != 1)",
                SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Check non-existence of 'sc' results columns in input data table */
    if (cpl_table_has_column(tabdat->tab[1], "scflux") != 0 ||
        cpl_table_has_column(tabdat->tab[1], "scdflux") != 0 ||
        cpl_table_has_column(tabdat->tab[1], "scmask") != 0) {
        sprintf(errtxt, "%s: sctarr *tabdat ('sc' results column(s) already "
                "exist(s))", SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Get and check number of data points */
    nrow = cpl_table_get_nrow(results);
    if (nrow == 0) {
        sprintf(errtxt, "%s: cpl_table *results", SC_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_NDA, "%s", errtxt);
    } else if (nrow != cpl_table_get_nrow(tabdat->tab[1])) {
        sprintf(errtxt, "%s: cpl_table *results != sctarr tabdat->tab[1]",
                SC_ERROR_IDG_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IDG, "%s", errtxt);
    }

    /* Copy scflux column with SKYCORR results into table with input data */
    cpl_table_duplicate_column(tabdat->tab[1], "scflux", results, "scflux");

    /* Copy scdflux column if present */
    if (cpl_table_has_column(results, "scdflux") == 1) {
        cpl_table_duplicate_column(tabdat->tab[1], "scdflux",
                                   results, "scdflux");
    }

    /* Copy scmask column if present */
    if (cpl_table_has_column(results, "scmask") == 1) {
        cpl_table_duplicate_column(tabdat->tab[1], "scmask",
                                   results, "scmask");
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_erasemaskcol(sctarr *tabdat,
                                    const cpl_parameterlist *parlist)
{
    /*!
     * Erases the integer mask column for SKYCORR in an ::sctarr structure.
     * The column was created by ::sc_conv_modtable.
     *
     * \b INPUT:
     * \param tabdat   ::sctarr structure with read data
     * \param parlist  general CPL parameter list
     *
     * \b OUTPUT:
     * \param tabdat   ::sctarr structure without integer mask column
     *
     * \b ERRORS:
     * - No data
     * - Invalid object structure
     */

    const cpl_parameter *p;
    char col_mask[SC_LENLINE+1] = "", col_imask[SC_LENLINE+1] = "";
    char errtxt[SC_MAXLEN] = "", col_renamed[SC_LENLINE+1] = "";
    int next = 0;

    /* Get name of integer mask column */
    p = cpl_parameterlist_find_const(parlist, "col_mask");
    strncpy(col_mask, cpl_parameter_get_string(p), SC_LENLINE+1);
    if (strcmp(col_mask, "NONE") == 0) {
        sprintf(col_imask, "%s_I", SC_DEFMASKCOL);
    } else {
        sprintf(col_imask, "%s_I", col_mask);
    }

    /* Check number of extensions */
    next = tabdat->next;
    if (next <= 0) {
        sprintf(errtxt, "%s: sctarr *tabdat", SC_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_NDA, "%s", errtxt);
    } else if (next > 1) {
        sprintf(errtxt, "%s: sctarr *tabdat (number of extensions != 1)",
                SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Erase integer mask column if it exists */
    if (cpl_table_has_column(tabdat->tab[1], col_imask) == 1) {
        cpl_table_erase_column(tabdat->tab[1], col_imask);
    }

    /* Rename a column with the same original name as the internal integer
       mask column (temporary suffix: "_orig") */
    sprintf(col_renamed, "%s_orig", col_imask);
    if (cpl_table_has_column(tabdat->tab[1], col_renamed) == 1) {
        cpl_table_name_column(tabdat->tab[1], col_renamed, col_imask);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_resultstarr2varr(scvarr *vecdat,
                                        const cpl_table *extnames,
                                        const sctarr *tabdat)
{
    /*!
     * Writes sky subtracted spectrum from SKYCORR into an ::scvarr structure
     * representing a FITS image of the same format as the input data file.
     *
     * \b INPUT:
     * \param vecdat    ::scvarr structure with data of input 1D FITS image
     * \param extnames  CPL table with FITS extension numbers and names
     * \param tabdat    ::sctarr structure with input data and SKYCORR
     *                  results
     *
     * \b OUTPUT:
     * \param vecdat    ::scvarr structure with applied sky subtraction
     *
     * \b ERRORS:
     * - No data
     * - Invalid object structure
     * - Invalid object value(s)
     * - Inconsistent data grids
     */

    char errtxt[SC_MAXLEN] = "";
    cpl_boolean hascol = CPL_TRUE;
    int next = 0, extn_flux = 0, extn_dflux = 0, extn_mask = 0, nrow = 0;
    int j = 0, i = 0;
    int *scmask = NULL;
    double maskval[2] = {0., 1.};
    double *flux = NULL, *dflux = NULL, *mask = NULL, *scflux = NULL;
    double *scdflux = NULL;

    /* Check table for extension numbers and names */
    if (cpl_table_get_nrow(extnames) != 4) {
        sprintf(errtxt, "%s: cpl_table *extnames (number of rows != 4)",
                SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Check column names */
    if (cpl_table_has_column(extnames, "col") != 1 ||
        cpl_table_has_column(extnames, "extn") != 1) {
        sprintf(errtxt, "%s: cpl_table *extnames (column names)",
                SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Check existence of data */
    if (cpl_table_get_column_max(extnames, "extn") == -1) {
        sprintf(errtxt, "%s: cpl_table *extnames (no extension with valid "
                "data)", SC_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }

    /* Get number of extensions */
    next = vecdat->next;

    /* Check extension numbers */
    if (cpl_table_get_column_max(extnames, "extn") > next) {
        sprintf(errtxt, "%s: cpl_table *extnames (incorrect extension "
                "number)", SC_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }

    /* Get extension numbers */
    extn_flux = cpl_table_get(extnames, "extn", 1, NULL);
    extn_dflux = cpl_table_get(extnames, "extn", 2, NULL);
    extn_mask = cpl_table_get(extnames, "extn", 3, NULL);

    /* Check number of extensions in results table */
    if (tabdat->next != 1) {
        sprintf(errtxt, "%s: sctarr *tabdat (number of extensions != 1)",
                SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Get number of data points */
    nrow = cpl_table_get_nrow(tabdat->tab[1]);
    if (nrow == 0) {
        sprintf(errtxt, "%s: sctarr *tabdat", SC_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_NDA, "%s", errtxt);
    }

    /* Check number of data points */
    for (j = 0; j <= next; j++) {
        if (cpl_vector_get_size(vecdat->vec[j]) != nrow) {
                sprintf(errtxt, "%s: scvarr vecdat->vec[%d] != "
                        "sctarr tabdat->tab[1]", SC_ERROR_IDG_TXT, j);
                return cpl_error_set_message(cpl_func, SC_ERROR_IDG, "%s",
                                             errtxt);
        }
    }

    /* Get the two (extreme) original mask values */
    if (sc_conv_getmaskval(maskval, tabdat, extnames) != CPL_ERROR_NONE) {
        return cpl_error_get_code();
    }

    /* Get pointers to image data */
    flux = cpl_vector_get_data(vecdat->vec[extn_flux]);
    if (extn_dflux >= 0) {
        dflux = cpl_vector_get_data(vecdat->vec[extn_dflux]);
    }
    if (extn_mask >= 0) {
        mask = cpl_vector_get_data(vecdat->vec[extn_mask]);
    }

    /* Get pointers to columns of results data table */
    scflux = cpl_table_get_data_double(tabdat->tab[1], "scflux");
    if (scflux == NULL) hascol = CPL_FALSE;
    if (extn_dflux >= 0) {
        scdflux = cpl_table_get_data_double(tabdat->tab[1], "scdflux");
        if (scdflux == NULL) hascol = CPL_FALSE;
    }
    if (extn_mask >= 0) {
        scmask = cpl_table_get_data_int(tabdat->tab[1], "scmask");
        if (scmask == NULL) hascol = CPL_FALSE;
    }
    if (hascol == CPL_FALSE) {
      cpl_table_dump_structure(tabdat->tab[1], NULL);
        sprintf(errtxt, "%s: sctarr *tabdat (required 'sc' column(s) not "
                "present)", SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Write sky-corrected flux, errors, and mask into flux image vectors */
    for (i = 0; i < nrow; i++) {
        flux[i] = scflux[i];
        if (extn_dflux >= 0) {
            dflux[i] = scdflux[i];
        }
        if (extn_mask >= 0) {
            /* Correct mask value if necessary (good -> bad) */
            if (mask[i] == maskval[1] && scmask[i] == 0) {
                mask[i] = maskval[0];
            }
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_getmaskval(double maskval[2], const sctarr *tabdat,
                                  const cpl_table *extnames)
{
    /*!
     * Gets the bad (maskval[0]) and good (maskval[1]) mask values from an
     * :sctarr structure containing the data of the FITS table prepared for
     * SKYCORR.
     *
     * \b INPUT:
     * \param tabdat    ::sctarr structure with input data
     * \param extnames  CPL table with FITS extension numbers and names
     *
     * \b OUTPUT:
     * \param maskval   array for bad and good mask value
     *
     * \b ERRORS:
     * - Invalid object structure
     * - No data
     */

    char errtxt[SC_MAXLEN] = "", col_mask[SC_LENLINE+1] = "";
    char col_imask[SC_LENLINE+1] = "";
    cpl_boolean hascol = CPL_TRUE;
    cpl_size imaskmaxpos = 0;
    int extn_mask = 0;
    double maskmin = 0., maskmax = 0.;

    /* Default mask values */
    maskval[0] = 0;
    maskval[1] = 1;

    /* Check table for extension numbers and names */
    if (cpl_table_get_nrow(extnames) != 4) {
        sprintf(errtxt, "%s: cpl_table *extnames (number of rows != 4)",
                SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Check column names */
    if (cpl_table_has_column(extnames, "col") != 1 ||
        cpl_table_has_column(extnames, "extn") != 1) {
        sprintf(errtxt, "%s: cpl_table *extnames (column names)",
                SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Get name and extension of input mask column */
    sprintf(col_mask, "%s", cpl_table_get_string(extnames, "col", 3));
    extn_mask = cpl_table_get_int(extnames, "extn", 3, NULL);

    /* Get name of integer mask column */
    if (extn_mask >= 0) {
        sprintf(col_imask, "%s_I", col_mask);
    } else {
        sprintf(col_imask, "%s_I", SC_DEFMASKCOL);
    }

    /* Check number of extensions in results table */
    if (tabdat->next != 1) {
        sprintf(errtxt, "%s: sctarr *tabdat (number of extensions != 1)",
                SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Check number of data points */
    if (cpl_table_get_nrow(tabdat->tab[1]) == 0) {
        sprintf(errtxt, "%s: sctarr *tabdat", SC_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_NDA, "%s", errtxt);
    }

    /* Check existence of mask columns in input data table */
    if (extn_mask >= 0) {
        if (cpl_table_has_column(tabdat->tab[1], col_mask) != 1) {
            hascol = CPL_FALSE;
        }
    }
    if (cpl_table_has_column(tabdat->tab[1], col_imask) != 1) {
        hascol = CPL_FALSE;
    }
    if (hascol == CPL_FALSE) {
        sprintf(errtxt, "%s: sctarr *tabdat (required input column(s) not "
                "present)", SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Get the two (extreme) original mask values */
    if (extn_mask >= 0) {
        maskmin = cpl_table_get_column_min(tabdat->tab[1], col_mask);
        maskmax = cpl_table_get_column_max(tabdat->tab[1], col_mask);
        cpl_table_get_column_maxpos(tabdat->tab[1], col_imask, &imaskmaxpos);
        if (cpl_table_get_int(tabdat->tab[1], col_imask, imaskmaxpos, NULL)
            == 1) {
            maskval[0] = maskmax;
            maskval[1] = maskmin;
        } else {
            maskval[0] = maskmin;
            maskval[1] = maskmax;
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_extract1d(const char *parfile_in,
                                 const char *scifilename,
                                 const char *skyfilename,
                                 const double minrelprofflux,
                                 const double lowpixfrac)
{
    /*!
     * \callgraph
     *
     * Converts two input 2D FITS images into compatible object and sky
     * reference 1D FITS images. Instead of two different FITS images, a
     * single file could also be used to produce object and sky spectrum. The
     * output file names and the FITS extension labels are taken from the
     * provided parameter file. The instrumental set-ups and the file
     * structures have to be identical for both input frames. For the input
     * parameters \e minrelprofflux and \e lowpixfrac, see
     * ::sc_conv_iarr2varr_sci and ::sc_conv_iarr2varr_sky, repectively.
     *
     * \b INPUT:
     * \param parfile_in      name of parameter file
     * \param scifilename     name of input 2D FITS image for object
     *                        extraction
     * \param skyfilename     name of input 2D FITS image for sky extraction
     * \param minrelprofflux  minimum relative profile flux (science spectrum)
     * \param lowpixfrac      column-specific fraction of pixels with flux
     *                        below selected pixel (sky spectrum)
     *
     * \b ERRORS:
     * - Unexpected file structure
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameterlist *parlist;
    cpl_parameter *p;
    cpl_table *extnames1, *extnames2;
    cpl_vector *selpix;
    sciarr imadat1, imadat2;
    scvarr vecdat1, vecdat2;
    char parfile[SC_MAXLEN] = "", errtxt[SC_MAXLEN] = "";
    char outscifilename[SC_MAXLEN] = "", outskyfilename[SC_MAXLEN] = "";
    int fitsformat = 0, i = 0;

    /* Read driver file */
    sc_basic_absfile(parfile, parfile_in);
    parlist = cpl_parameterlist_new();
    if ((status = sc_par_readfile(parlist, parfile)) != CPL_ERROR_NONE) {
        cpl_parameterlist_delete(parlist);
        return status;
    }

    /* Get type of science file */
    if ((status = sc_conv_checkfitsformat(&fitsformat, scifilename)) !=
        CPL_ERROR_NONE) {
        cpl_parameterlist_delete(parlist);
        return status;
    }

    /* Check for correct FITS format (2D image) of science file */
    if (fitsformat != 3) {
        cpl_parameterlist_delete(parlist);
        sprintf(errtxt, "%s: %s (no 2D FITS image)", SC_ERROR_UFS_TXT,
                scifilename);
        return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s", errtxt);
    }

    /* Get file type of sky file */
    if ((status = sc_conv_checkfitsformat(&fitsformat, skyfilename)) !=
        CPL_ERROR_NONE) {
        cpl_parameterlist_delete(parlist);
        return status;
    }

    /* Check for correct FITS format (2D image) of sky file */
    if (fitsformat != 3) {
        cpl_parameterlist_delete(parlist);
        sprintf(errtxt, "%s: %s (no 2D FITS image)", SC_ERROR_UFS_TXT,
                skyfilename);
        return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s", errtxt);
    }

    /* Read 2D FITS science image */
    if ((status = sc_conv_iarr_read(&imadat1, scifilename)) !=
        CPL_ERROR_NONE) {
        cpl_parameterlist_delete(parlist);
        return status;
    }

    /* Read 2D FITS sky image */
    if ((status = sc_conv_iarr_read(&imadat2, skyfilename)) !=
        CPL_ERROR_NONE) {
        cpl_parameterlist_delete(parlist);
        sc_conv_iarr_delete(&imadat1);
        return status;
    }

    /* Check correspondence of pixel number in both images */
    if (cpl_image_get_size_x(imadat1.ima[0]) !=
        cpl_image_get_size_x(imadat2.ima[0]) ||
        cpl_image_get_size_y(imadat1.ima[0]) !=
        cpl_image_get_size_y(imadat2.ima[0])) {
        cpl_parameterlist_delete(parlist);
        sc_conv_iarr_delete(&imadat1);
        sc_conv_iarr_delete(&imadat2);
        sprintf(errtxt, "%s: %s != %s (different pixel number)",
                SC_ERROR_UFS_TXT, scifilename, skyfilename);
        return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s", errtxt);
    }

    /* Get names and positions of FITS extensions for science image */
    extnames1 = cpl_table_new(0);
    if ((status = sc_conv_setextnames_iarr(extnames1, &imadat1, parlist)) !=
         CPL_ERROR_NONE) {
        cpl_parameterlist_delete(parlist);
        cpl_table_delete(extnames1);
        sc_conv_iarr_delete(&imadat1);
        sc_conv_iarr_delete(&imadat2);
        return status;
    }

    /* Get names and positions of FITS extensions for sky image */
    extnames2 = cpl_table_new(0);
    if ((status = sc_conv_setextnames_iarr(extnames2, &imadat2, parlist)) !=
         CPL_ERROR_NONE) {
        cpl_parameterlist_delete(parlist);
        cpl_table_delete(extnames1);
        cpl_table_delete(extnames2);
        sc_conv_iarr_delete(&imadat1);
        sc_conv_iarr_delete(&imadat2);
        return status;
    }

    /* Compare table data */
    for (i = 0; i < cpl_table_get_nrow(extnames1); i++) {
        if (strcmp(cpl_table_get_string(extnames1, "col", i),
                   cpl_table_get_string(extnames2, "col", i)) != 0 ||
            cpl_table_get_int(extnames1, "extn", i, NULL) !=
            cpl_table_get_int(extnames2, "extn", i, NULL)) {
            cpl_parameterlist_delete(parlist);
            cpl_table_delete(extnames1);
            cpl_table_delete(extnames2);
            sc_conv_iarr_delete(&imadat1);
            sc_conv_iarr_delete(&imadat2);
            sprintf(errtxt, "%s: %s != %s (extension names and/or positions)",
                    SC_ERROR_UFS_TXT, scifilename, skyfilename);
            return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s",
                                         errtxt);
        }
    }

    /* Convert sciarr structure (2D) into science scvarr structure (1D) */
    selpix = cpl_vector_new(1);
    if ((status = sc_conv_iarr2varr_sci(&vecdat1, selpix, &imadat1, extnames1,
                                        minrelprofflux)) != CPL_ERROR_NONE) {
        cpl_parameterlist_delete(parlist);
        cpl_table_delete(extnames1);
        cpl_table_delete(extnames2);
        cpl_vector_delete(selpix);
        sc_conv_iarr_delete(&imadat1);
        sc_conv_iarr_delete(&imadat2);
        return status;
    }

    /* Convert sciarr structure (2D) into sky scvarr structure (1D) */
    if ((status = sc_conv_iarr2varr_sky(&vecdat2, &imadat2, extnames2,
                                        lowpixfrac, selpix)) !=
        CPL_ERROR_NONE) {
        cpl_parameterlist_delete(parlist);
        cpl_table_delete(extnames1);
        cpl_table_delete(extnames2);
        cpl_vector_delete(selpix);
        sc_conv_iarr_delete(&imadat1);
        sc_conv_iarr_delete(&imadat2);
        sc_conv_varr_delete(&vecdat1);
        return status;
    }

    /* Free allocated memory */
    cpl_table_delete(extnames1);
    cpl_table_delete(extnames2);
    cpl_vector_delete(selpix);
    sc_conv_iarr_delete(&imadat1);
    sc_conv_iarr_delete(&imadat2);

    /* Get output file names from parameter file */
    p = cpl_parameterlist_find(parlist, "scispec");
    strncpy(outscifilename, cpl_parameter_get_string(p), SC_MAXLEN);
    p = cpl_parameterlist_find(parlist, "skyspec");
    strncpy(outskyfilename, cpl_parameter_get_string(p), SC_MAXLEN);

    /* Write scvarr structure into 1D FITS science image */
    if ((status = sc_conv_varr_write(outscifilename, &vecdat1)) !=
        CPL_ERROR_NONE) {
        cpl_parameterlist_delete(parlist);
        sc_conv_varr_delete(&vecdat1);
        sc_conv_varr_delete(&vecdat2);
        return status;
    }

    /* Write scvarr structure into 1D FITS sky image */
    if ((status = sc_conv_varr_write(outskyfilename, &vecdat2)) !=
        CPL_ERROR_NONE) {
        cpl_parameterlist_delete(parlist);
        sc_conv_varr_delete(&vecdat1);
        sc_conv_varr_delete(&vecdat2);
        return status;
    }

    /* Free allocated memory */
    cpl_parameterlist_delete(parlist);
    sc_conv_varr_delete(&vecdat1);
    sc_conv_varr_delete(&vecdat2);

    /* Return error code of last error */
    return cpl_error_get_code();
}


cpl_error_code sc_conv_iarr2varr_sci(scvarr *vecdat, cpl_vector *selpix,
                                     const sciarr *imadat,
                                     const cpl_table *extnames,
                                     const double minrelprofflux)
{
    /*!
     * Fills CPL vectors and CPL property lists of a ::scvarr structure by
     * data read from a 2D FITS image with optional extensions for error,
     * mask, and wavelength, and added in spatial direction. The image data is
     * provided by a ::sciarr structure. The extension numbers and names are
     * given by a CPL table. Extensions are not read if names are missing or
     * wrong.
     *
     * The flux summation only considers pixels with a flux of the spatial
     * slit profile relative to the maximum value above a given minimum.
     * Morover, only good pixels are added if a mask extension is provided.
     * If pixels along the spatial direction are not considered, the summed
     * flux is corrected by means of a model of the spatial flux distribution.
     * The latter is derived from a projection of the 2D spectrum. If errors
     * exist, they are added quadratically. Possible correlations of the
     * pixels are not considered.
     *
     * The effective number of pixels that are added up for each wavelength is
     * output by means of an CPL vector. These data can be used to extract a
     * corresponding 1D sky spectrum.
     *
     * \b INPUT:
     * \param selpix          empty CPL vector
     * \param imadat          ::sciarr structure with data of 2D FITS image
     * \param extnames        CPL table with FITS extension numbers and names
     * \param minrelprofflux  minimum relative profile flux
     *
     * \b OUTPUT:
     * \param vecdat          ::scvarr structure with summed data of 2D FITS
     *                        image
     * \param selpix          vector of effective pixel numbers
     *
     * \b ERRORS:
     * - Invalid input parameter(s)
     * - No data
     * - Invalid object structure
     * - Invalid object value(s)
     * - see subroutines
     */

    cpl_vector *row, *prof, *relsum;
    char errtxt[SC_MAXLEN] = "";
    int next = 0, h = 0, nx = 0, ny = 0, j = 0, i = 0, qual = 0, spix = 0;
    const int *extn;
    double maxerr = 0., maskval[2] = {0., 0.}, profmin = 0., profmax = 0.;
    double proflim = 0., profval = 0., sum = 0., scale = 0., mask = 0.;
    double val = 0.;

    /* Check input parameter */
    if (minrelprofflux < 0. || minrelprofflux > 1.) {
        sprintf(errtxt, "%s: double minrelprofflux (value outside [0,1])",
                SC_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IIP, "%s", errtxt);
    }

    /* Get number of extensions */
    next = imadat->next;

    /* Initialise CPL vectors and CPL property lists for content of sciarr
       structure (memory allocation) */
    sc_conv_varr_init(vecdat, next);
    for (h = 0; h <= next; h++) {
        vecdat->vec[h] = cpl_vector_new(1);
    }

    /* Transfer FITS header data and delete keywords specific for 2D images */
    for (h = 0; h <= next; h++) {
        vecdat->head[h] = cpl_propertylist_duplicate(imadat->head[h]);
        cpl_propertylist_erase(vecdat->head[h], "CRPIX2");
        cpl_propertylist_erase(vecdat->head[h], "CRVAL2");
        cpl_propertylist_erase(vecdat->head[h], "CDELT2");
        cpl_propertylist_erase(vecdat->head[h], "CTYPE2");
    }

    /* Check table for extension numbers and names */
    if (cpl_table_get_nrow(extnames) != 4) {
        sc_conv_varr_delete(vecdat);
        sprintf(errtxt, "%s: cpl_table *extnames (number of rows != 4)",
                SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Check column names */
    if (cpl_table_has_column(extnames, "col") != 1 ||
        cpl_table_has_column(extnames, "extn") != 1) {
        sc_conv_varr_delete(vecdat);
        sprintf(errtxt, "%s: cpl_table *extnames (column names)",
                SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Check existence of data */
    if (cpl_table_get_column_max(extnames, "extn") == -1) {
        sc_conv_varr_delete(vecdat);
        sprintf(errtxt, "%s: cpl_table *extnames (no extension with valid "
                "data)", SC_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }

    /* Check extension numbers */
    if (cpl_table_get_column_max(extnames, "extn") > next) {
        sc_conv_varr_delete(vecdat);
        sprintf(errtxt, "%s: cpl_table *extnames (incorrect extension "
                "number)", SC_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }

    /* Get pointer to extension numbers */
    extn = cpl_table_get_data_int_const(extnames, "extn");

    /* Check image size */
    nx = cpl_image_get_size_x(imadat->ima[0]);
    ny = cpl_image_get_size_y(imadat->ima[0]);
    if (nx <= 0 || ny <= 0) {
        sc_conv_varr_delete(vecdat);
        sprintf(errtxt, "%s: sciarr *imadat", SC_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_NDA, "%s", errtxt);
    } else if (nx > 0 && ny > 0 && next > 0) {
        for (h = 1; h <= next; h++) {
            if (cpl_image_get_size_x(imadat->ima[h]) != nx ||
                cpl_image_get_size_y(imadat->ima[h]) != ny) {
                sc_conv_varr_delete(vecdat);
                sprintf(errtxt, "%s: sciarr *imadat (image size differs for "
                        "different extensions", SC_ERROR_IOS_TXT);
                return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s",
                                             errtxt);
            }
        }
    }

    /* Get maximum error if available */
    if (extn[2] >= 0) {
        maxerr = cpl_image_get_max(imadat->ima[extn[2]]);
    }

    /* Get meaning of mask values (identification of good and bad pixels) */
    if (extn[3] >= 0) {
        if (sc_basic_getmaskval_image(maskval, imadat->ima[extn[3]]) !=
            CPL_ERROR_NONE) {
            return SC_ERROR_IOV;
        }
    }

    /* Create row (x) and profile (y) vector */
    row = cpl_vector_new(nx);
    prof = cpl_vector_new(ny);

    /* Get profile along slit */
    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            cpl_vector_set(row, i,
                           cpl_image_get(imadat->ima[extn[1]], i+1, j+1,
                                         &qual));
        }
        cpl_vector_set(prof, j, cpl_vector_get_median(row));
    }

    /* Delete temporary row vector */
    cpl_vector_delete(row);

    /* Subtract minimum flux from profile function */
    profmin = cpl_vector_get_min(prof);
    cpl_vector_subtract_scalar(prof, profmin);

    /* Avoid pixels with negligible flux contribution */
    profmax = cpl_vector_get_max(prof);
    proflim = minrelprofflux * profmax;
    for (j = 0; j < ny; j++) {
         profval = cpl_vector_get(prof, j);
         if (profval < proflim) {
             cpl_vector_set(prof, j, 0.);
         }
    }

    /* Normalise profile function to get weights */
    for (sum = 0., j = 0; j < ny; j++) {
        sum += cpl_vector_get(prof, j);
    }
    if (sum == 0.) {
        sc_conv_varr_delete(vecdat);
        cpl_vector_delete(prof);
        sprintf(errtxt, "%s: sciarr *imadat (flux sum = 0)",
                SC_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }
    cpl_vector_divide_scalar(prof, sum);

    /* Calculate correction function for columns with bad pixels */
    relsum = cpl_vector_new(nx);
    if (extn[3] >= 0) {
        for (i = 0; i < nx; i++) {
            for (sum = 0., j = 0; j < ny; j++) {
                if (cpl_image_get(imadat->ima[extn[3]], i+1, j+1, &qual) <
                    maskval[1] + SC_TOL) {
                    sum += cpl_vector_get(prof, j);
                }
            }
            cpl_vector_set(relsum, i, sum);
        }
    } else {
        cpl_vector_fill(relsum, 1.);
    }

    /* Set size of vecdat vectors */
    for (h = 0; h <= next; h++) {
        cpl_vector_set_size(vecdat->vec[h], nx);
    }

    /* Set size of selpix vector and initialise all elements */
    cpl_vector_set_size(selpix, nx);
    cpl_vector_fill(selpix, 0.);

    /* Convert sciarr to scvarr structure (unidentified extensions are
       handled like flux) */

    for (h = 0; h <= next; h++) {
        for (i = 0; i < nx; i++) {

            /* Get scaling factor for pixel flux */
            scale = cpl_vector_get(relsum, i);

            if (h == extn[3]) {

                /* Set mask values (bad value only if no good pixel in image
                   column) */
                if (scale < 2 * SC_TOL) {
                    cpl_vector_set(vecdat->vec[extn[3]], i, maskval[0]);
                } else {
                    cpl_vector_set(vecdat->vec[extn[3]], i, maskval[1]);
                }

            } else {

                /* Add values of good pixels */

                for (sum = 0., spix = 0, j = 0; j < ny; j++) {

                    /* Skip pixels with negligible flux contribution */
                    profval = cpl_vector_get(prof, j);
                    if (profval == 0.) {
                        continue;
                    }

                    /* Get mask value for pixel */
                    if (extn[3] >= 0) {
                        mask = cpl_image_get(imadat->ima[extn[3]], i+1, j+1,
                                             &qual);
                    } else {
                        mask = maskval[1];
                    }

                    /* Get pixel value */
                    val = cpl_image_get(imadat->ima[h], i+1, j+1, &qual);
                    if (val <= 0. && h == extn[0]) {
                        sc_conv_varr_delete(vecdat);
                        cpl_vector_delete(prof);
                        cpl_vector_delete(relsum);
                        sprintf(errtxt, "%s: sciarr *imadat "
                                "(wavelength <= 0)", SC_ERROR_IOV_TXT);
                        return cpl_error_set_message(cpl_func, SC_ERROR_IOV,
                                                     "%s", errtxt);
                    }

                    /* Consider good pixels only if not wavelength */
                    if (mask == maskval[1] || h == extn[0]) {
                        /* Add pixel values */
                        if (h == extn[2]) {
                            /* Squared summation of error pixels */
                            sum += val * val;
                        } else {
                            sum += val;
                        }
                        /* Number of good pixels */
                        if (h == extn[1]) {
                            spix++;
                        }
                    }

                }

                /* Save effective pixel number */
                if (h == extn[1]) {
                    if (scale > 0.) {
                        cpl_vector_set(selpix, i, (double) spix / scale);
                    } else {
                        cpl_vector_set(selpix, i, (double) spix);
                    }
                }

                /* Get resulting value depending on kind of data */
                if (h == extn[0]) {
                    /* Wavelength */
                    sum /= (double) ny;
                } else if (h == extn[2]) {
                    /* Error */
                    if (sum > 0. && scale > 0.) {
                        sum = sqrt(sum) / scale;
                    } else {
                        sum = maxerr;
                    }
                } else {
                    /* Flux */
                    if (scale > 0.) {
                        sum /= scale;
                    }
                }

                /* Write resulting value into output vector */
                cpl_vector_set(vecdat->vec[h], i, sum);

            }

        }
    }

    /* Free allocated memory */
    cpl_vector_delete(prof);
    cpl_vector_delete(relsum);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_iarr2varr_sky(scvarr *vecdat, const sciarr *imadat,
                                     const cpl_table *extnames,
                                     const double lowpixfrac,
                                     const cpl_vector *selpix)
{
    /*!
     * Fills CPL vectors and CPL property lists of a ::scvarr structure by
     * data read from a 2D FITS image with optional extensions for error,
     * mask, and wavelength, and median-averaged in spatial direction. The
     * image data is provided by a ::sciarr structure. The extension numbers
     * and names are given by a CPL table. Extensions are not read if names
     * are missing or wrong.
     *
     * The median calculation only considers good pixels provided that a mask
     * extension is provided. The input parameter \e lowpixfrac allows one to
     * take a pixel in the sorted pixel list which deviates from the middle
     * one (0.5). Values below 0.5 (full available range: 0 to 1) reduce the
     * influence of spectra of bright objects on the calculation of the
     * background median flux. The output flux is scaled to equal the sky flux
     * for the input effective number of pixels provided by a CPL vector for
     * each wavelength.
     *
     * \b INPUT:
     * \param imadat      ::sciarr structure with data of 2D FITS image
     * \param extnames    CPL table with FITS extension numbers and names
     * \param lowpixfrac  column-specific fraction of pixels with flux below
     *                    selected pixel
     * \param selpix      CPL vector of effective pixel numbers
     *
     * \b OUTPUT:
     * \param vecdat      ::scvarr structure with median-averaged data of 2D
     *                    FITS image
     *
     * \b ERRORS:
     * - Invalid input parameter(s)
     * - Invalid object value(s)
     * - No data
     * - Invalid object structure
     * - see subroutines
     */

    cpl_vector *yvec, *ysel;
    char errtxt[SC_MAXLEN] = "";
    int next = 0, h = 0, nx = 0, ny = 0, j = 0, i = 0, ngood = 0, qual = 0;
    int jmed = 0;
    const int *extn;
    double maskval[2] = {0., 0.};
    double mask = 0., val = 0., med = 0., nmed = 0.;

    /* Check input parameter */
    if (lowpixfrac < 0. || lowpixfrac > 1.) {
        sprintf(errtxt, "%s: double lowfracpix (value outside [0,1])",
                SC_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IIP, "%s", errtxt);
    }

    /* Check validity of CPL vector of effective pixel numbers */
    if (cpl_vector_get_max(selpix) > 0.) {
    } else {
        sprintf(errtxt, "%s: cpl_vector *selpix (invalid elements or all 0)",
                SC_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }

    /* Get number of extensions */
    next = imadat->next;

    /* Initialise CPL vectors and CPL property lists for content of sciarr
       structure (memory allocation) */
    sc_conv_varr_init(vecdat, next);
    for (h = 0; h <= next; h++) {
        vecdat->vec[h] = cpl_vector_new(1);
    }

    /* Transfer FITS header data and delete keywords specific for 2D images */
    for (h = 0; h <= next; h++) {
        vecdat->head[h] = cpl_propertylist_duplicate(imadat->head[h]);
        cpl_propertylist_erase(vecdat->head[h], "CRPIX2");
        cpl_propertylist_erase(vecdat->head[h], "CRVAL2");
        cpl_propertylist_erase(vecdat->head[h], "CDELT2");
        cpl_propertylist_erase(vecdat->head[h], "CTYPE2");
    }

    /* Check table for extension numbers and names */
    if (cpl_table_get_nrow(extnames) != 4) {
        sc_conv_varr_delete(vecdat);
        sprintf(errtxt, "%s: cpl_table *extnames (number of rows != 4)",
                SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Check column names */
    if (cpl_table_has_column(extnames, "col") != 1 ||
        cpl_table_has_column(extnames, "extn") != 1) {
        sc_conv_varr_delete(vecdat);
        sprintf(errtxt, "%s: cpl_table *extnames (column names)",
                SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Check existence of data */
    if (cpl_table_get_column_max(extnames, "extn") == -1) {
        sc_conv_varr_delete(vecdat);
        sprintf(errtxt, "%s: cpl_table *extnames (no extension with valid "
                "data)", SC_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }

    /* Check extension numbers */
    if (cpl_table_get_column_max(extnames, "extn") > next) {
        sc_conv_varr_delete(vecdat);
        sprintf(errtxt, "%s: cpl_table *extnames (incorrect extension "
                "number)", SC_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }

    /* Get pointer to extension numbers */
    extn = cpl_table_get_data_int_const(extnames, "extn");

    /* Check image size */
    nx = cpl_image_get_size_x(imadat->ima[0]);
    ny = cpl_image_get_size_y(imadat->ima[0]);
    if (nx <= 0 || ny <= 0) {
        sc_conv_varr_delete(vecdat);
        sprintf(errtxt, "%s: sciarr *imadat", SC_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_NDA, "%s", errtxt);
    } else if (nx > 0 && ny > 0 && next > 0) {
        for (h = 1; h <= next; h++) {
            if (cpl_image_get_size_x(imadat->ima[h]) != nx ||
                cpl_image_get_size_y(imadat->ima[h]) != ny) {
                sc_conv_varr_delete(vecdat);
                sprintf(errtxt, "%s: sciarr *imadat (image size differs for "
                        "different extensions", SC_ERROR_IOS_TXT);
                return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s",
                                             errtxt);
            }
        }
    }

    /* Compare size of image in x direction and size of selpix vector */
    if (cpl_vector_get_size(selpix) != nx) {
        sc_conv_varr_delete(vecdat);
        sprintf(errtxt, "%s: cpl_vector *selpix (size) != sciarr *imadat "
                "(x size)", SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Get meaning of mask values (identification of good and bad pixels) */
    if (extn[3] >= 0) {
        if (sc_basic_getmaskval_image(maskval, imadat->ima[extn[3]]) !=
            CPL_ERROR_NONE) {
           sc_conv_varr_delete(vecdat);
           return SC_ERROR_IOV;
        }
    }

    /* Create temporary y column vectors */
    yvec = cpl_vector_new(ny);
    ysel = cpl_vector_new(ny);

    /* Set size of output vectors */
    for (h = 0; h <= next; h++) {
        cpl_vector_set_size(vecdat->vec[h], nx);
    }

    /* Convert sciarr to scvarr structure (unidentified extensions are
       handled like flux) */

    for (h = 0; h <= next; h++) {
        for (i = 0; i < nx; i++) {

            for (ngood = 0, j = 0; j < ny; j++) {

                /* Get mask value for pixel */
                if (extn[3] >= 0) {
                    mask = cpl_image_get(imadat->ima[extn[3]], i+1, j+1,
                                         &qual);
                } else {
                    mask = maskval[1];
                }

                /* Count good pixels */
                if (mask == maskval[1]) {
                    ngood++;
                }

                /* Skip rest of loop in the case of mask extension */
                if (h == extn[3]) {
                    continue;
                }

                /* Get pixel value */
                val = cpl_image_get(imadat->ima[h], i+1, j+1, &qual);
                if (val <= 0. && h == extn[0]) {
                    sc_conv_varr_delete(vecdat);
                    cpl_vector_delete(yvec);
                    cpl_vector_delete(ysel);
                    sprintf(errtxt, "%s: sciarr *imadat (wavelength <= 0)",
                            SC_ERROR_IOV_TXT);
                    return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s",
                                                 errtxt);
                }

                /* Set value for good and bad pixels */
                cpl_vector_set(yvec, j, val);
                if (mask == maskval[1]) {
                    /* good pixel */
                    cpl_vector_set(ysel, j, val);
                } else {
                    /* bad pixel */
                    cpl_vector_set(ysel, j, HUGE_VAL);
                }

            }

            /* Set mask values (bad value only if no good pixel in image
               column) */
            if (h == extn[3]) {
                if (ngood == 0) {
                    cpl_vector_set(vecdat->vec[extn[3]], i, maskval[0]);
                } else {
                    cpl_vector_set(vecdat->vec[extn[3]], i, maskval[1]);
                }
                continue;
            }

            /* Get resulting value depending on kind of data */
            if (ngood > 0) {
                cpl_vector_sort(ysel, CPL_SORT_ASCENDING);
                jmed = ceil(lowpixfrac * ngood) - 1;
                if (jmed < 0) {
                    jmed = 0;
                }
                med = cpl_vector_get(ysel, jmed);
            } else {
                cpl_vector_sort(yvec, CPL_SORT_ASCENDING);
                jmed = ceil(lowpixfrac * ny) - 1;
                if (jmed < 0) {
                    jmed = 0;
                }
                med = cpl_vector_get(yvec, jmed);
            }

            /* Rough error estimate */
            if (h == extn[2]) {
                nmed = 2. * ngood * SC_MIN(lowpixfrac, 1. - lowpixfrac);
                if (nmed > 0.) {
                    med /= sqrt(nmed);
                }
            }

            /* Provide fluxes and errors for given effective pixel number */
            if (h == extn[1] || h == extn[2]) {
                med *= cpl_vector_get(selpix, i);
            }


            /* Write resulting value into output vector */
            cpl_vector_set(vecdat->vec[h], i, med);

        }
    }

    /* Free allocated memory */
    cpl_vector_delete(yvec);
    cpl_vector_delete(ysel);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_ascii_read(sctarr *tabdat, const char *filename,
                                  const cpl_array *colnames)
{
    /*!
     * Fills CPL table of a ::sctarr structure by data read from an ASCII
     * file. The number and names of the expected columns are provided by a
     * CPL array. Header lines in the ASCII file are allowed if they are
     * marked by '#'.
     *
     * \b INPUT:
     * \param filename  path and name of input ASCII file
     * \param colnames  CPL array of column names for reading of ASCII file
     *
     * \b OUTPUT:
     * \param tabdat    ::sctarr structure with data of ASCII file
     *
     * \b ERRORS:
     * - No data
     * - File opening failed
     * - Unexpected file structure
     */

    FILE *stream;
    scpar val[SC_MAXPAR];
    char errtxt[SC_MAXLEN] = "", str[SC_LENLINE+2] = "";
    const char **col;
    int ncolmin = 0, nrec = 0, j = 0, i = 0;
    int ncol0 = SC_MAXPAR, ncol = SC_MAXPAR;

    /* Get size of column name array */
    ncolmin = cpl_array_get_size(colnames);
    if (ncolmin == 0 || colnames == NULL) {
        sprintf(errtxt, "%s: cpl_array *colnames", SC_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_NDA, "%s", errtxt);
    }

    /* Get pointer for column name array */
    col = cpl_array_get_data_string_const(colnames);

    /* Check file existence and open ASCII file */
    if ((stream = fopen(filename, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", SC_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, SC_ERROR_FOF, "%s", errtxt);
    }

    /* Find number of data lines */
    while (fgets(str, SC_LENLINE+2, stream) != NULL) {
        if (str[0] == '#' || str[0] == '\n') {
        } else if (isdigit(str[0]) || isspace(str[0]) || str[0] == '-') {
            nrec++;
        } else {
            fclose(stream);
            sprintf(errtxt, "%s: %s (unexpected first character at line)",
                    SC_ERROR_UFS_TXT, filename);
            return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s",
                                         errtxt);
        }
    }
    rewind(stream);

    /* No data points */
    if (nrec == 0) {
        fclose(stream);
        sprintf(errtxt, "%s: %s (no data)",
                SC_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s", errtxt);
    }

    /* Number of values per line */
    sc_basic_readline(stream, val, &ncol0);
    if (ncol0 < ncolmin) {
        fclose(stream);
        sprintf(errtxt, "%s: %s (too low number of columns)",
                SC_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s", errtxt);
    }
    rewind(stream);

    /* Initialise CPL table and CPL property list for content of ASCII file
       (memory allocation) */
    sc_conv_tarr_init(tabdat, 1);

    /* Create empty CPL property lists */
    tabdat->head[0] = cpl_propertylist_new();
    tabdat->head[1] = cpl_propertylist_new();

    /* Create CPL tables of required size (put data in first extension of
       output file) */
    tabdat->tab[0] = cpl_table_new(0);
    tabdat->tab[1] = cpl_table_new(nrec);

    /* Create required table columns */
    for (j = 0; j < ncolmin; j++) {
        cpl_table_new_column(tabdat->tab[1], col[j], CPL_TYPE_DOUBLE);
    }

    /* Read spectral data from file and write it to CPL table */

    for (i = 0; i < nrec; i++) {

        sc_basic_readline(stream, val, &ncol);

        if (ncol != ncol0) {
            sc_conv_tarr_delete(tabdat);
            fclose(stream);
            sprintf(errtxt, "%s: %s (unexpected number of values at line)",
                    SC_ERROR_UFS_TXT, filename);
            return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s",
                                         errtxt);
        }

        for (j = 0; j < ncolmin; j++) {
            cpl_table_set(tabdat->tab[1], col[j], i, val[j].d);
        }

    }

    /* Close ASCII file */
    fclose(stream);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_ascii_write(const char *filename, const sctarr *tabdat)
{
    /*!
     * Writes the table data of an ::sctarr structure into an ASCII file. Only
     * one extension is allowed. The column names are written into a header
     * line starting with '#'. The routine supports the column types STRING,
     * INT, FLOAT, and DOUBLE.
     *
     * \b INPUT:
     * \param filename  path and name of output ASCII file
     * \param tabdat    ::sctarr structure with tabulated data
     *
     * \b ERRORS:
     * - Invalid object structure
     * - No data
     * - File opening failed
     * - see subroutines
     */

    FILE *stream;
    cpl_type coltype;
    cpl_array *colnames = NULL;
    char errtxt[SC_MAXLEN] = "", str[SC_MAXLEN] = "";
    char strcomp[SC_LENLINE+2] = "";
    char **col;
    int nrow = 0, ncol = 0, i = 0, j = 0;

    /* Check number of extensions */
    if (tabdat->next != 1) {
        sprintf(errtxt, "%s: sctarr *tabdat (number of extensions != 1)",
                SC_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s", errtxt);
    }

    /* Get and check number of data points */
    nrow = cpl_table_get_nrow(tabdat->tab[1]);
    if (nrow == 0) {
        sprintf(errtxt, "%s: sctarr *tabdat", SC_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_NDA, "%s", errtxt);
    }

    /* Get column names in input table */
    colnames = cpl_table_get_column_names(tabdat->tab[1]);
    ncol = cpl_array_get_size(colnames);

    /* Get pointer to array */
    col = cpl_array_get_data_string(colnames);

    /* Open output ASCII file */
    if ((stream = fopen(filename, "w+")) == NULL) {
        cpl_array_delete(colnames);
        sprintf(errtxt, "%s: %s", SC_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, SC_ERROR_FOF, "%s", errtxt);
    }

    /* Build string of column names */
    sprintf(str, "#");
    for (j = 0; j < ncol; j++) {
        sprintf(strcomp, " %s", col[j]);
        strncat(str, strcomp, strlen(strcomp));
    }

    /* Write column names into ASCII file */
    fprintf(stream, "%s\n", str);

    /* Write data into ASCII file */

    for (i = 0; i < nrow; i++) {

        /* Build string of data values */
        str[0] = '\0';
        for (j = 0; j < ncol; j++) {
            coltype = cpl_table_get_column_type(tabdat->tab[1], col[j]);
            if (coltype == CPL_TYPE_STRING) {
                sprintf(strcomp, "%s", cpl_table_get_string(tabdat->tab[1],
                                                            col[j], i));
            } else if (coltype == CPL_TYPE_INT) {
                sprintf(strcomp, "%d", cpl_table_get_int(tabdat->tab[1],
                                                         col[j], i, NULL));
            } else if (coltype == CPL_TYPE_FLOAT) {
                sprintf(strcomp, "%f", cpl_table_get_float(tabdat->tab[1],
                                                           col[j], i, NULL));
            } else if (coltype == CPL_TYPE_DOUBLE) {
                sprintf(strcomp, "%e", cpl_table_get_double(tabdat->tab[1],
                                                            col[j], i, NULL));
            } else {
                fclose(stream);
                cpl_array_delete(colnames);
                sprintf(errtxt, "%s: sctarr tabdat->tab[1] (unsupported "
                        "column type)", SC_ERROR_IOS_TXT);
                return cpl_error_set_message(cpl_func, SC_ERROR_IOS, "%s",
                                             errtxt);
            }
            strncat(str, strcomp, strlen(strcomp));
            if (j != ncol-1) {
                sprintf(strcomp, " ");
                strncat(str, strcomp, strlen(strcomp));
            }
        }

        /* Write data string into ASCII file */
        fprintf(stream, "%s\n", str);

   }

    /* Close ASCII file */
    fclose(stream);

    /* Free allocated memory */
    cpl_array_delete(colnames);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_tarr_init(sctarr *tabdat, const int next)
{
    /*!
     * Initialises array of CPL tables and CPL property lists as ::sctarr
     * structure depending on the input number of FITS extensions.
     *
     * \b INPUT:
     * \param next    number of FITS extensions
     *
     * \b OUTPUT:
     * \param tabdat  empty ::sctarr structure for FITS table data
     *
     * \b ERRORS:
     * - none
     */

    tabdat->next = next;
    tabdat->head = cpl_calloc(next+1, sizeof(cpl_propertylist *));
    tabdat->tab = cpl_calloc(next+1, sizeof(cpl_table *));

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_tarr_read(sctarr *tabdat, const char *filename)
{
    /*!
     * Reads a FITS table with an arbritrary number of extensions and puts the
     * data into an ::sctarr structure consisting of an array of CPL tables
     * and CPL property lists.
     *
     * \b INPUT:
     * \param filename  path and name of input FITS table
     *
     * \b OUTPUT:
     * \param tabdat    ::sctarr structure with FITS table data
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     */

    char errtxt[SC_MAXLEN] = "";
    int next = 0, fitsformat = 0, i = 0;

    /* Get number of extensions and check file */
    next = cpl_fits_count_extensions(filename);
    if (next < 0) {
        sprintf(errtxt, "%s: %s", SC_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, SC_ERROR_FOF, "%s", errtxt);
    } else if (next == 0) {
        sprintf(errtxt, "%s: %s (number of FITS extensions = 0)",
                SC_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s", errtxt);
    }

    /* Check FITS format */
    sc_conv_checkfitsformat(&fitsformat, filename);
    if (fitsformat != 1) {
        sprintf(errtxt, "%s: %s (no FITS table)", SC_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s", errtxt);
    }

    /* Initialise CPL tables and CPL property lists for content of FITS
       file (memory allocation) */
    sc_conv_tarr_init(tabdat, next);

    /* Read FITS file */
    for (i = 0; i <= next; i++) {
        tabdat->head[i] = cpl_propertylist_load(filename, i);
        if (i == 0) {
            tabdat->tab[i] = cpl_table_new(0);
        } else {
            tabdat->tab[i] = cpl_table_load(filename, i, 0);
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_tarr_write(const char *filename, const sctarr *tabdat)
{
    /*!
     * Writes data of an ::sctarr structure into a FITS table. The number of
     * FITS extensions depends on the content of the ::sctarr structure.
     *
     * \b INPUT:
     * \param filename  path and name of output FITS table
     * \param tabdat    ::sctarr structure with FITS table data
     *
     * \b ERRORS:
     * - none
     */

    int next = 0, i = 0;

    /* Get number of FITS extensions */
    next = tabdat->next;

    /* Write FITS table with next extensions */
    for (i = 0; i < next; i++) {
        if (i == 0) {
            cpl_table_save(tabdat->tab[i+1], tabdat->head[0],
                           tabdat->head[i+1], filename, CPL_IO_CREATE);
        } else {
            cpl_table_save(tabdat->tab[i+1], NULL,
                           tabdat->head[i+1], filename, CPL_IO_EXTEND);
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_tarr_delete(sctarr *tabdat)
{
    /*!
     * Deletes an ::sctarr structure, which contains an array of CPL tables
     * and CPL property lists.
     *
     * \b INPUT:
     * \param tabdat  ::sctarr structure with FITS table data
     *
     * \b ERRORS:
     * - none
     */

    int next = 0, i = 0;

    if (tabdat == NULL) {
        return CPL_ERROR_NONE;
    }

    next = tabdat->next;

    for (i = 0; i <= next; i++) {
        cpl_propertylist_delete(tabdat->head[i]);
        cpl_table_delete(tabdat->tab[i]);
    }

    cpl_free(tabdat->head);
    cpl_free(tabdat->tab);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_varr_init(scvarr *vecdat, const int next)
{
    /*!
     * Initialises array of CPL vectors and CPL property lists as ::scvarr
     * structure depending on the input number of FITS extensions.
     *
     * \b INPUT:
     * \param next    number of FITS extensions
     *
     * \b OUTPUT:
     * \param vecdat  empty ::scvarr structure for 1D FITS image data
     *
     * \b ERRORS:
     * - none
     */

    vecdat->next = next;
    vecdat->head = cpl_calloc(next+1, sizeof(cpl_propertylist *));
    vecdat->vec = cpl_calloc(next+1, sizeof(cpl_vector *));

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_varr_read(scvarr *vecdat, const char *filename)
{
    /*!
     * Reads a 1D FITS image with an arbritrary number of extensions and puts
     * the data into an ::scvarr structure consisting of an array of CPL
     * vectors and CPL property lists.
     *
     * \b INPUT:
     * \param filename  path and name of input 1D FITS image
     *
     * \b OUTPUT:
     * \param vecdat    ::scvarr structure with 1D FITS image data
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     */

    char errtxt[SC_MAXLEN] = "";
    int next = 0, fitsformat = 0, i = 0;

    /* Get number of extensions and check file */
    next = cpl_fits_count_extensions(filename);
    if (next < 0) {
        sprintf(errtxt, "%s: %s", SC_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, SC_ERROR_FOF, "%s", errtxt);
    }

    /* Check FITS format */
    sc_conv_checkfitsformat(&fitsformat, filename);
    if (fitsformat != 2) {
        sprintf(errtxt, "%s: %s (no 1D FITS image)", SC_ERROR_UFS_TXT,
                filename);
        return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s", errtxt);
    }

    /* Initialise CPL tables and CPL property lists for content of FITS
       file (memory allocation) */
    sc_conv_varr_init(vecdat, next);

    /* Read FITS file */
    for (i = 0; i <= next; i++) {
        vecdat->head[i] = cpl_propertylist_load(filename, i);
        vecdat->vec[i] = cpl_vector_load(filename, i);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_varr_write(const char *filename, const scvarr *vecdat)
{
    /*!
     * Writes data of an ::scvarr structure into a 1D FITS image. The number
     * of FITS extensions depends on the content of the ::scvarr structure.
     *
     * \b INPUT:
     * \param filename  path and name of output 1D FITS image
     * \param vecdat    ::scvarr structure with 1D FITS image data
     *
     * \b ERRORS:
     * - none
     */

    int next = 0, i = 0;

    /* Get number of FITS extensions */
    next = vecdat->next;

    /* Write 1D FITS image with next extensions */
    for (i = 0; i <= next; i++) {
        if (i == 0) {
            cpl_vector_save(vecdat->vec[i], filename, CPL_TYPE_FLOAT,
                            vecdat->head[i], CPL_IO_CREATE);
        } else {
            cpl_vector_save(vecdat->vec[i], filename, CPL_TYPE_FLOAT,
                            vecdat->head[i], CPL_IO_EXTEND);
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_varr_delete(scvarr *vecdat)
{
    /*!
     * Deletes an ::scvarr structure, which contains an array of CPL vectors
     * and CPL property lists.
     *
     * \b INPUT:
     * \param vecdat  ::scvarr structure with 1D FITS image data
     *
     * \b ERRORS:
     * - none
     */

    int next = 0, i = 0;

    if (vecdat == NULL) {
        return CPL_ERROR_NONE;
    }

    next = vecdat->next;

    for (i = 0; i <= next; i++) {
        cpl_propertylist_delete(vecdat->head[i]);
        cpl_vector_delete(vecdat->vec[i]);
    }

    cpl_free(vecdat->head);
    cpl_free(vecdat->vec);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_iarr_init(sciarr *imadat, const int next)
{
    /*!
     * Initialises array of CPL images and CPL property lists as ::sciarr
     * structure depending on the input number of FITS extensions.
     *
     * \b INPUT:
     * \param next    number of FITS extensions
     *
     * \b OUTPUT:
     * \param imadat  empty ::sciarr structure for 2D FITS image data
     *
     * \b ERRORS:
     * - none
     */

    imadat->next = next;
    imadat->head = cpl_calloc(next+1, sizeof(cpl_propertylist *));
    imadat->ima = cpl_calloc(next+1, sizeof(cpl_image *));

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_iarr_read(sciarr *imadat, const char *filename)
{
    /*!
     * Reads a 2D FITS image with an arbritrary number of extensions and puts
     * the data into an ::sciarr structure consisting of an array of CPL
     * images and CPL property lists.
     *
     * \b INPUT:
     * \param filename  path and name of input 2D FITS image
     *
     * \b OUTPUT:
     * \param imadat    ::sciarr structure with 2D FITS image data
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     */

    char errtxt[SC_MAXLEN] = "";
    int next = 0, fitsformat = 0, i = 0;

    /* Get number of extensions and check file */
    next = cpl_fits_count_extensions(filename);
    if (next < 0) {
        sprintf(errtxt, "%s: %s", SC_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, SC_ERROR_FOF, "%s", errtxt);
    }

    /* Check FITS format */
    sc_conv_checkfitsformat(&fitsformat, filename);
    if (fitsformat != 3) {
        sprintf(errtxt, "%s: %s (no 2D FITS image)",SC_ERROR_UFS_TXT,
                filename);
        return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s", errtxt);
    }

    /* Initialise CPL tables and CPL property lists for content of FITS
       file (memory allocation) */
    sc_conv_iarr_init(imadat, next);

    /* Read FITS file */
    for (i = 0; i <= next; i++) {
        imadat->head[i] = cpl_propertylist_load(filename, i);
        imadat->ima[i] = cpl_image_load(filename, CPL_TYPE_UNSPECIFIED, 0, i);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_iarr_write(const char *filename, const sciarr *imadat)
{
    /*!
     * Writes data of an ::sciarr structure into a 2D FITS image. The number
     * of FITS extensions depends on the content of the ::sciarr structure.
     *
     * \b INPUT:
     * \param filename  path and name of output 2D FITS image
     * \param imadat    ::sciarr structure with 2D FITS image data
     *
     * \b ERRORS:
     * - none
     */

    int next = 0, i = 0;

    /* Get number of FITS extensions */
    next = imadat->next;

    /* Write 2D FITS image with next extensions */
    for (i = 0; i <= next; i++) {
        if (i == 0) {
            cpl_image_save(imadat->ima[i], filename, CPL_TYPE_UNSPECIFIED,
                           imadat->head[i], CPL_IO_CREATE);
        } else {
            cpl_image_save(imadat->ima[i], filename, CPL_TYPE_UNSPECIFIED,
                           imadat->head[i], CPL_IO_EXTEND);
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_conv_iarr_delete(sciarr *imadat)
{
    /*!
     * Deletes an ::sciarr structure, which contains an array of CPL images
     * and CPL property lists.
     *
     * \b INPUT:
     * \param imadat  ::sciarr structure with 2D FITS image data
     *
     * \b ERRORS:
     * - none
     */

    int next = 0, i = 0;

    if (imadat == NULL) {
        return CPL_ERROR_NONE;
    }

    next = imadat->next;

    for (i = 0; i <= next; i++) {
        cpl_propertylist_delete(imadat->head[i]);
        cpl_image_delete(imadat->ima[i]);
    }

    cpl_free(imadat->head);
    cpl_free(imadat->ima);

    return CPL_ERROR_NONE;
}

/**@}*/
