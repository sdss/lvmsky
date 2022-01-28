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
 * \ingroup sm_general
 */

/**@{*/

/*!
 * \file sm_general.c
 *
 * Basic routines used for the sky model
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  22 Sep 2009
 * \date   06 Oct 2015
 */


/*****************************************************************************
 *                                INCLUDES                                   *
 ****************************************************************************/

#include <sm_general.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

cpl_error_code sm_spec_malloc(smspec *spec, const int size)
{
    /*!
     * Allocates memory for an ::smspec structure of given size.
     * Assumes symmetric errors (see ::smspec).
     *
     * \b INPUT:
     * \param spec  ::smspec structure (no memory allocated)
     * \param size  number of data points
     *
     * \b OUTPUT:
     * \param spec  ::smspec structure with allocated memory (all values = 0)
     *
     * \b ERRORS:
     * - ISM: Insufficient memory
     */

    char errtxt[SM_MAXLEN+1];

    /* Existence of symmetric errors by default */
    spec->type = 3;

    /* Allocate memory */

    spec->n = size;
    spec->dat = (smdat *) calloc(spec->n, sizeof(smdat));

    if (spec->dat == NULL && spec->n > 0) {
        spec->n = 0;
        sprintf(errtxt, "%s: smspec *spec", SM_ERROR_ISM_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_ISM, "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_create(smspec *outspec, const double limlam[2],
                              const double dlam)
{
    /*!
     * Initialisation of a wavelength grid that is characterised by constant
     * bin size. The spectrum type 3 is set (see ::smspec).
     *
     * \b INPUT:
     * \param limlam   lower and upper limit of wavelength range
     * \param dlam     bin size
     *
     * \b OUTPUT:
     * \param outspec  ::smspec structure that provides wavelengths, fluxes
     *                 (= 0.), and flux errors (= 0.) of a spectrum
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     */

    char errtxt[SM_MAXLEN+1];
    int i;

    /* Test input parameters */

    if (limlam[1] <= limlam[0] || dlam <= 0.) {
        /* Return spectrum with zero data points */
        outspec->type = 3;
        outspec->n = 0;
        sprintf(errtxt, "%s: limlam[1] <= limlam[0] || dlam <= 0 "
                "(wavelength grid)", SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Get number of data points and allocate memory */
    if ((int) sm_spec_malloc(outspec, ((limlam[1] - limlam[0]) /
            dlam + 0.5) + 1) == (int) SM_ERROR_ISM) {
        return SM_ERROR_ISM;
    }

    /* Write wavelength grid and flux = 0 in structure smspec */

    outspec->dat[0].lam = limlam[0];
    outspec->dat[0].flux = 0.;
    outspec->dat[0].dflux1 = 0.;
    outspec->dat[0].dflux2 = 0.;
    for (i = 1; i < outspec->n; i++) {
        outspec->dat[i].lam = outspec->dat[i-1].lam + dlam;
        outspec->dat[i].flux = 0.;
        outspec->dat[i].dflux1 = 0.;
        outspec->dat[i].dflux2 = 0.;
   }

   return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_read(smspec *spec, const char *filename)
{
    /*!
     * Fills ::smspec structure by data read from a file.
     * Each line of the file must consist of wavelength and flux.
     * The presence of one or two optional error columns is detected
     * automatically.
     * Header lines are allowed if they are marked by #.
     *
     * \b INPUT:
     * \param filename  name of data file
     *
     * \b OUTPUT:
     * \param spec      output spectrum with read values
     *
     * \b ERRORS:
     * - FOF: File opening failed
     * - UFS: Unexpected file structure
     * - CPL: Access out of range [warning]
     */

    FILE *stream;
    char errtxt[SM_MAXLEN+1], str[SM_LENLINE+2];
    smbool isoutrange = F;
    int nhead = 0, nrec = 0, ncol0 = 3, i = 0, ncol = 0;
    double lam = 0, flux = 0, dflux1 = 0, dflux2 = 0;

    spec->type = 3;
    spec->n = 0;

    /* Check file existence */
    if ((stream = fopen(filename, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", SM_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, SM_ERROR_FOF, "%s", errtxt);
    }

    /* Find number of header and data lines */
    while (fgets(str, SM_LENLINE+2, stream) != NULL) {
        if (str[0] == '#') {
            nhead++;
            if (nrec != 0) {
                fclose(stream);
                sprintf(errtxt, "%s: %s (comment line in data part)",
                        SM_ERROR_UFS_TXT, filename);
                return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s",
                                             errtxt);
            }
        } else if (isdigit(str[0]) || isspace(str[0]) || str[0] == '-') {
            if (nrec == 0) {
                /* Number of values per line */
                ncol0 = sscanf(str, "%le %le %le %le",
                               &lam, &flux, &dflux1, &dflux2);
                if (ncol0 == 0) {
                    fclose(stream);
                    sprintf(errtxt, "%s: %s (empty line)",
                        SM_ERROR_UFS_TXT, filename);
                    return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s",
                                                 errtxt);
                }
            }
            nrec++;
        } else {
            fclose(stream);
            sprintf(errtxt, "%s: %s (unexpected first character at line)",
                    SM_ERROR_UFS_TXT, filename);
            return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s",
                                         errtxt);
        }
    }
    rewind(stream);

    /* No data points */
    if (nrec == 0) {
        fclose(stream);
        sprintf(errtxt, "%s: %s (no data)",
                SM_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s", errtxt);
    }

    /* Allocate memory (set spec->type = 3) */
    if ((int) sm_spec_malloc(spec, nrec) == (int) SM_ERROR_ISM) {
        fclose(stream);
        return SM_ERROR_ISM;
    }

    /* Existence of flux and error data? */
    spec->type = ncol0;

    /* Skip header */
    for (i = 0; i < nhead; i++) {
        if (fgets(str, SM_LENLINE+2, stream)) {};
    }

    /* Read spectral data */

    for (i = 0; i < nrec; i++) {

        if (spec->type == 1) {
            ncol = fscanf(stream, "%le", &lam);
            flux = 0;
            dflux1 = 0;
            dflux2 = 0;
        } else if (spec->type == 2) {
            ncol = fscanf(stream, "%le %le", &lam, &flux);
            dflux1 = 0;
            dflux2 = 0;
        } else if (spec->type == 3) {
            ncol = fscanf(stream, "%le %le %le", &lam, &flux, &dflux1);
            dflux2 = dflux1;
        } else {
            ncol = fscanf(stream, "%le %le %le %le",
                          &lam, &flux, &dflux1, &dflux2);
        }

        if (ncol != ncol0) {
            sprintf(errtxt, "%s: %s (unexpected number of values at line)",
                    SM_ERROR_UFS_TXT, filename);
            return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s",
                                         errtxt);

        }

        spec->dat[i].lam = lam;
        spec->dat[i].flux = flux;

        if (dflux1 < 0) {
            spec->dat[i].dflux1 = HUGE_VAL;
            isoutrange = T;
        } else {
            spec->dat[i].dflux1 = dflux1;
        }

        if (dflux2 < 0) {
            spec->dat[i].dflux2 = HUGE_VAL;
            isoutrange = T;
        } else {
            spec->dat[i].dflux2 = dflux2;
        }

    }

    fclose(stream);

    if (isoutrange == T) {
        sprintf(errtxt, "%s (negative error(s)) [warning]", filename);
        return cpl_error_set_message(cpl_func, CPL_ERROR_ACCESS_OUT_OF_RANGE,
                                     "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_readrange(smspec *spec, const char *filename,
                                 const double limlam[2], const int step)
{
    /*!
     * Reading of a spectrum from file in a given wavelength range.
     * The wavelength range of the extracted spectrum tends to be wider than
     * the interval between the provided limits and depends on the given step
     * size.
     * The presence of one or two optional error columns is detected
     * automatically.
     * Header lines are allowed if they are marked by #.
     *
     * \b INPUT:
     * \param filename  name of file with wavelength and flux (and error)
     *                  columns
     * \param limlam    wavelength limits for extraction
     * \param step      step in file lines for search of given wavelength
     *                  range (optimum: square root of total number of lines)
     *
     * \b OUTPUT:
     * \param spec      spectrum of structure ::smspec
     *
     * ERRORS:
     * - IIP: Invalid input parameter(s)
     * - FOF: File opening failed
     * - UFS: Unexpected file structure
     * - CPL: Access out of range [warning]
     */

    FILE *stream;
    char errtxt[SM_MAXLEN+1], str[SM_LENLINE+2];
    smbool isoutrange = F;
    int l = 0, ipos[2] = {-1, -1}, ospos = 0, ncol0 = 3, spos = 0, spos0 = 0;
    int i = 0, ncol = 0;
    double lam = 0, flux = 0, dflux1 = 0, dflux2 = 0;

    spec->type = 3;
    spec->n = 0;

    /* Reasonable wavelength limits and step size? */

    if (limlam[1] <= limlam[0] || step < 1) {
        sprintf(errtxt, "%s: limlam[1] <= limlam[0] (wavelength limits) || "
                "step < 1 (line jump for searching)" , SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Check file existence */

    if ((stream = fopen(filename, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", SM_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, SM_ERROR_FOF, "%s", errtxt);
    }

    /* Skip header (marked by "#") */

    while (fgets(str, SM_LENLINE+2, stream) != NULL) {
        if (str[0] == '#') {
            ospos = ftell(stream);
        } else if (isdigit(str[0]) || isspace(str[0]) || str[0] == '-') {
            fseek(stream, ospos, SEEK_SET);
            /* Number of values per line */
            ncol0 = sscanf(str, "%le %le %le %le",
                           &lam, &flux, &dflux1, &dflux2);
            if (ncol0 == 0) {
                fclose(stream);
                sprintf(errtxt, "%s: %s (empty line)", SM_ERROR_UFS_TXT,
                        filename);
                return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s",
                                             errtxt);
            }
            break;
        } else {
            fclose(stream);
            sprintf(errtxt, "%s: %s (unexpected first character at line)",
                    SM_ERROR_UFS_TXT, filename);
            return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s",
                                         errtxt);
        }
    }

    /* Read file line by line and check wavelength every [step] lines.
       Identify file positions of limiting wavelengths. */

    while (fgets(str,SM_LENLINE+2,stream) != NULL) {

        if (l % step == 0) {

            sscanf(str, "%le", &lam);

            if (ipos[0] < 0 && lam > limlam[0]) {
                if (l == 0) {
                    ipos[0] = l;
                } else {
                    ipos[0] = l - step;
                }
                spos0 = ospos;
            }

            if (ipos[0] >= 0 && lam >= limlam[1]) {
                if (l == 0) {
                    ipos[0] = -1;
                } else {
                    ipos[1] = l;
                }
                break;
            }

        }

        if ((l + 1) % step == 0) {
            ospos = spos;
            spos = ftell(stream);
        }

        l++;

    }

    if (ipos[0] >= 0 && ipos[1] < 0) {
        ipos[1] = l - 1;
    }

    /* Next steps for overlapping wavelength ranges only */

    if (ipos[0] >= 0 && ipos[1] >= 0) {

        /* Derivation of number of wavelengths in output spectrum and
           allocation of memory (set spec->type = 3) */

        if ((int) sm_spec_malloc(spec, ipos[1] - ipos[0] + 1) ==
                (int) SM_ERROR_ISM) {
            fclose(stream);
            return SM_ERROR_ISM;
        }

        /* Existence of flux and error data? */
        spec->type = ncol0;

        /* Jump to position of first relevant data point in file and read
           data */

        fseek(stream, spos0, SEEK_SET);

        for (i = 0; i < spec->n; i++) {

            if (spec->type == 1) {
                ncol = fscanf(stream, "%le", &lam);
                flux = 0;
                dflux1 = 0;
                dflux2 = 0;
            } else if (spec->type == 2) {
                ncol = fscanf(stream, "%le %le", &lam, &flux);
                dflux1 = 0;
                dflux2 = 0;
            } else if (spec->type == 3) {
                ncol = fscanf(stream, "%le %le %le", &lam, &flux, &dflux1);
                dflux2 = dflux1;
            } else {
                ncol = fscanf(stream, "%le %le %le %le",
                              &lam, &flux, &dflux1, &dflux2);
            }

            if (ncol != ncol0) {
                sprintf(errtxt,
                        "%s: %s (unexpected number of values at line)",
                        SM_ERROR_UFS_TXT, filename);
                return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s",
                                             errtxt);
            }

            spec->dat[i].lam = lam;
            spec->dat[i].flux = flux;

            if (dflux1 < 0) {
                spec->dat[i].dflux1 = HUGE_VAL;
                isoutrange = T;
            } else {
                spec->dat[i].dflux1 = dflux1;
            }

            if (dflux2 < 0) {
                spec->dat[i].dflux2 = HUGE_VAL;
                isoutrange = T;
            } else {
                spec->dat[i].dflux2 = dflux2;
            }

        }

    }

    fclose(stream);

    if (isoutrange == T) {
        sprintf(errtxt, "%s (negative error(s)) [warning]", filename);
        return cpl_error_set_message(cpl_func, CPL_ERROR_ACCESS_OUT_OF_RANGE,
                                     "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_readcpl(smspec *spec, const cpl_table *cpltab)
{
    /*!
     * Reads ::smspec structure from CPL table with one to four columns.
     * The column IDs "lam", "flux", "dflux" (for symmetric errors), "dflux1"
     * (lower error), and "dflux2" (upper error) are mandatory.
     *
     * \b INPUT:
     * \param cpltab  CPL table with one to four columns named
     *                "lam", "flux", "dflux"/"dflux1", or "dflux2"
     *
     * \b OUTPUT:
     * \param spec    spectrum of structure ::smspec
     *
     * \b ERRORS:
     * - NDA: No data
     * - IOS: Invalid object structure
     * - CPL: Access out of range [warning]
     */

    char errtxt[SM_MAXLEN+1];
    char colid[][7] = {"lam", "flux", "dflux", "dflux1", "dflux2"};
    char colid0[7];
    smbool isoutrange = F;
    int nrow, ncol, j, i;
    const double *lam, *flux, *dflux1, *dflux2;

    spec->n = 0;
    spec->type = 3;

    /* Get number of table rows */
    nrow = cpl_table_get_nrow(cpltab);
    if (nrow <= 0) {
        sprintf(errtxt, "%s: cpl_table *cpltab", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    }

    /* Check number of table columns */
    ncol = cpl_table_get_ncol(cpltab);
    if (ncol < 1 || ncol > 4) {
        sprintf(errtxt,
                "%s: cpl_table *cpltab (number of columns < 1 or > 4)",
                SM_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IOS, "%s", errtxt);
    }

    /* Check column names */
    for (j = 0; j < 5; j++) {
        if (cpl_table_has_column(cpltab, colid[j]) != 1) {
            if ((j < ncol && ncol <= 3) ||
                (j != 2 && ncol == 4)) {
                sprintf(errtxt, "%s: cpl_table *cpltab (invalid column name)",
                        SM_ERROR_IOS_TXT);
                return cpl_error_set_message(cpl_func, SM_ERROR_IOS, "%s",
                                             errtxt);
            }
        }
    }

    /* Allocate memory (set spec->type = 3) */
    if ((int) sm_spec_malloc(spec, nrow) == (int) SM_ERROR_ISM) {
        return SM_ERROR_ISM;
    }

    /* Existence of flux and error data? */
    spec->type = ncol;

    /* Transfer lam column to smspec structure */
    lam = cpl_table_get_data_double_const(cpltab, colid[0]);
    for (i = 0; i < spec->n; i++) {
        spec->dat[i].lam = lam[i];
    }

    /* Transfer flux column to smspec structure if it exists */
    if (spec->type >= 2) {
        flux = cpl_table_get_data_double_const(cpltab, colid[1]);
        for (i = 0; i < spec->n; i++) {
            spec->dat[i].flux = flux[i];
        }
    }

    /* Transfer dflux/dflux1 column to smspec structure if it exists */
    if (spec->type >= 3) {
        if (spec->type == 3) {
            strcpy(colid0, colid[2]);
        } else {
            strcpy(colid0, colid[3]);
        }
        dflux1 = cpl_table_get_data_double_const(cpltab, colid0);
        for (i = 0; i < spec->n; i++) {
            spec->dat[i].dflux1 = dflux1[i];
            if (spec->dat[i].dflux1 < 0) {
                spec->dat[i].dflux1 = HUGE_VAL;
                isoutrange = T;
            }
        }
    }

    /* Transfer dflux2 column to smspec structure if it exists */
    if (spec->type == 4) {
        dflux2 = cpl_table_get_data_double_const(cpltab, colid[4]);
        for (i = 0; i < spec->n; i++) {
            spec->dat[i].dflux2 = dflux2[i];
            if (spec->dat[i].dflux2 < 0) {
                spec->dat[i].dflux2 = HUGE_VAL;
                isoutrange = T;
            }
        }
    }

    /* Copy dflux1 to dflux2 if only dflux1 exists */
    if (spec->type == 3) {
        for (i = 0; i < spec->n; i++) {
            spec->dat[i].dflux2 = spec->dat[i].dflux1;
        }
    }

    if (isoutrange == T) {
        sprintf(errtxt, "cpl_table *cpltab (negative error(s)) [warning]");
        return cpl_error_set_message(cpl_func, CPL_ERROR_ACCESS_OUT_OF_RANGE,
                                     "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_readfits(smspec *spec, const char *filename)
{
    /*!
     * Reads ::smspec structure from FITS file by means of CPL
     *
     * \b INPUT:
     * \param filename  name of FITS file (extensions "fits" or "mt")
     *
     * \b OUTPUT:
     * \param spec      spectrum of structure ::smspec
     *
     * \b ERRORS:
     * - IFE: Invalid file name extension
     * - FOF: File opening failed
     * - UFS: Unexpected file structure
     */

    cpl_table *cpltab;
    cpl_array *colnames;
    char errtxt[SM_MAXLEN+1];
    const char *col;
    char colid[][7] = {"lam", "flux", "dflux", "dflux1", "dflux2"};
    char colid0[7];
    int ncol = 0;

    spec->n = 0;

    /* Check file name extension */
    if (strstr(filename, ".fits") == NULL &&
        strstr(filename, ".mt") == NULL) {
        sprintf(errtxt, "%s: %s (neither 'fits' nor 'mt')", SM_ERROR_IFE_TXT,
                filename);
        return cpl_error_set_message(cpl_func, SM_ERROR_IFE, "%s", errtxt);
    }

    /* Load FITS file into CPL table */
    cpltab = cpl_table_load(filename, 1, 0);
    if (cpltab == NULL) {
        sprintf(errtxt, "%s: %s", SM_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, SM_ERROR_FOF, "%s", errtxt);
    }

    /* Check number of table columns */
    ncol = cpl_table_get_ncol(cpltab);
    if (ncol < 1 || ncol > 4) {
        sprintf(errtxt, "%s: %s (number of columns < 1 or > 4)",
                SM_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s", errtxt);
    }

    /* Make sure that the CPL table columns have the correct names */

    colnames = cpl_table_get_column_names(cpltab);
    if (ncol >= 1 && ncol <= 4) {
        col = cpl_array_get_string(colnames, 0);
        if (strncmp(col, colid[0], strlen(colid[0])) != 0) {
            cpl_table_name_column(cpltab, cpl_array_get_string(colnames, 0),
                                  colid[0]);
        }
    }
    if (ncol >= 2 && ncol <= 4) {
        col = cpl_array_get_string(colnames, 1);
        if (strncmp(col, colid[1], strlen(colid[1])) != 0) {
            cpl_table_name_column(cpltab, cpl_array_get_string(colnames, 1),
                                  colid[1]);
        }
    }
    if (ncol >= 3 && ncol <= 4) {
        col = cpl_array_get_string(colnames, 2);
        if (ncol == 3) {
            strcpy(colid0, colid[2]);
        } else {
            strcpy(colid0, colid[3]);
        }
        if (strncmp(col, colid0, strlen(colid0)) != 0) {
            cpl_table_name_column(cpltab, cpl_array_get_string(colnames, 2),
                                  colid0);
        }
    }
    if (ncol == 4) {
        col = cpl_array_get_string(colnames, 3);
        if (strncmp(col, colid[4], strlen(colid[4])) != 0) {
            cpl_table_name_column(cpltab, cpl_array_get_string(colnames, 3),
                                  colid[4]);
        }
    }
    cpl_array_delete(colnames);

    /* Transfer table elements to smspec structure */
    sm_spec_readcpl(spec, cpltab);

    /* Delete temporary CPL table */
    cpl_table_delete(cpltab);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_readfitsrange(smspec *spec, const char *filename,
                                     const double limlam[2], const int step)
{
    /*!
     * Reading of a spectrum from FITS file in a given wavelength range by
     * means of CPL.
     * The wavelength range of the extracted spectrum tends to be wider than
     * the interval between the provided limits and depends on the given step
     * size.
     *
     * \b INPUT:
     * \param filename  name of FITS file (extensions "fits or "mt")
     * \param limlam    wavelength limits for extraction
     * \param step      step in FITS table rows for search of given wavelength
     *                  range (optimum: square root of total number of rows)
     *
     * \b OUTPUT:
     * \param spec      spectrum of structure ::smspec
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     * - IFE: Invalid file name extension
     * - FOF: File opening failed
     * - UFS: Unexpected file structure
     */

    cpl_table *temptab, *cpltab;
    cpl_array *colnames;
    cpl_propertylist *header;
    cpl_property *prop;
    char errtxt[SM_MAXLEN+1];
    const char *col;
    char colid[][7] = {"lam", "flux", "dflux", "dflux1", "dflux2"};
    char colid0[7];
    int nrow = 0, i = 0, ipos[2] = {-1, -1}, npos = 0, ncol = 0;

    spec->n = 0;

    /* Reasonable wavelength limits and step size? */
    if (limlam[1] <= limlam[0] || step < 1) {
        sprintf(errtxt, "%s: limlam[1] <= limlam[0] (wavelength limits) || "
                "step < 1 (line jump for searching)", SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Check file name extension */
    if (strstr(filename, ".fits") == NULL &&
        strstr(filename, ".mt") == NULL) {
        sprintf(errtxt, "%s: %s (neither 'fits' nor 'mt')", SM_ERROR_IFE_TXT,
                filename);
        return cpl_error_set_message(cpl_func, SM_ERROR_IFE, "%s", errtxt);
    }

    /* Read FITS table header */
    header = cpl_propertylist_load(filename, 1);
    if (header == NULL) {
        sprintf(errtxt, "%s: %s", SM_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, SM_ERROR_FOF, "%s", errtxt);
    }

    /* Read FITS keyword NAXIS */
    prop = cpl_propertylist_get_property(header, "NAXIS2");
    if (prop == NULL) {
        cpl_propertylist_delete(header);
        sprintf(errtxt, "%s: %s (keyword NAXIS not found)",
                SM_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s", errtxt);
    }
    nrow = cpl_property_get_int(prop);

    /* Delete temporary CPL property list */
    cpl_propertylist_delete(header);

    /* Start caching of FITS file */
    cpl_fits_set_mode(CPL_FITS_START_CACHING);

    /* Check wavelength every [step] file rows to identify positions of
       limiting wavelengths */

    while ((temptab = cpl_table_load_window(filename, 1, 0, NULL, i, 1))
           != NULL || i == 0) {

        if (ipos[0] < 0 &&
            cpl_table_get_double(temptab, "lam", 0, NULL) > limlam[0]) {
            if (i == 0) {
                ipos[0] = i;
            } else {
                ipos[0] = i - step;
            }
        }

        if (ipos[0] >= 0 &&
            cpl_table_get_double(temptab, "lam", 0, NULL) >= limlam[1]) {
            if (i == 0) {
                ipos[0] = -1;
            } else {
                ipos[1] = i;
            }
            break;
        }

        if (i == nrow - 1) {
            break;
        }

        cpl_table_delete(temptab);

        i += step;
        if (i >= nrow) {
            i = nrow - 1;
        }

    }

    if (ipos[0] >= 0 && ipos[1] < 0) {
        ipos[1] = i;
    }

    if (ipos[0] < 0 || ipos[1] < 0) {
        /* No overlap of wavelength ranges */
        cpl_fits_set_mode(CPL_FITS_STOP_CACHING);
        return CPL_ERROR_NONE;
    }

    /* Load FITS file into CPL table */
    npos = ipos[1] - ipos[0] + 1;
    cpltab = cpl_table_load_window(filename, 1, 0, NULL, ipos[0], npos);

    /* Stop caching of FITS file */
    cpl_fits_set_mode(CPL_FITS_STOP_CACHING);

    /* Check number of table columns */
    ncol = cpl_table_get_ncol(cpltab);
    if (ncol < 1 || ncol > 4) {
        sprintf(errtxt, "%s: %s (number of columns < 1 or > 4)",
                SM_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s", errtxt);
    }

    /* Make sure that the CPL table columns have the correct names */

    colnames = cpl_table_get_column_names(cpltab);
    if (ncol >= 1 && ncol <= 4) {
        col = cpl_array_get_string(colnames, 0);
        if (strncmp(col, colid[0], strlen(colid[0])) != 0) {
            cpl_table_name_column(cpltab, cpl_array_get_string(colnames, 0),
                                  colid[0]);
        }
    }
    if (ncol >= 2 && ncol <= 4) {
        col = cpl_array_get_string(colnames, 1);
        if (strncmp(col, colid[1], strlen(colid[1])) != 0) {
            cpl_table_name_column(cpltab, cpl_array_get_string(colnames, 1),
                                  colid[1]);
        }
    }
    if (ncol >= 3 && ncol <= 4) {
        col = cpl_array_get_string(colnames, 2);
        if (ncol == 3) {
            strcpy(colid0, colid[2]);
        } else {
            strcpy(colid0, colid[3]);
        }
        if (strncmp(col, colid0, strlen(colid0)) != 0) {
            cpl_table_name_column(cpltab, cpl_array_get_string(colnames, 2),
                                  colid0);
        }
    }
    if (ncol == 4) {
        col = cpl_array_get_string(colnames, 3);
        if (strncmp(col, colid[4], strlen(colid[4])) != 0) {
            cpl_table_name_column(cpltab, cpl_array_get_string(colnames, 3),
                                  colid[4]);
        }
    }
    cpl_array_delete(colnames);

    /* Transfer table elements to smspec structure */
    sm_spec_readcpl(spec, cpltab);

    /* Delete temporary CPL tables */
    cpl_table_delete(temptab);
    cpl_table_delete(cpltab);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_copy(smspec *outspec, const smspec *inspec)
{
    /*!
     * Copies ::smspec structure
     *
     * \b INPUT:
     * \param inspec   input spectrum
     *
     * \b OUTPUT:
     * \param outspec  output spectrum
     *
     * \b ERRORS:
     * - none
     */

    int i;

    /* Get number of data points and allocate memory */
    if ((int) sm_spec_malloc(outspec, inspec->n) == (int) SM_ERROR_ISM) {
        return SM_ERROR_ISM;
    }

    /* Copy type number of spectrum */
    outspec->type = inspec->type;

    /* Copy wavelengths and fluxes */

    if (outspec->n > 0) {
        for (i = 0; i < outspec->n; i++) {
            outspec->dat[i].lam = inspec->dat[i].lam;
            outspec->dat[i].flux = inspec->dat[i].flux;
            if (outspec->type >= 3) {
                outspec->dat[i].dflux1 = inspec->dat[i].dflux1;
                outspec->dat[i].dflux2 = inspec->dat[i].dflux2;
            }
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_compgrids(const smspec *spec1, const smspec *spec2)
{
    /*!
     * Tests agreement of wavelength grids of two smspec structures.
     * Returns \f$\ne 0\f$ for disagreement.
     *
     * \b INPUT:
     * \param spec1  first spectrum
     * \param spec2  second spectrum
     *
     * \b ERRORS:
     * - IDG: Inconsistent data grids
     */

    char errtxt[SM_MAXLEN+1];
    int i;
    double llam, ulam;

    if ((spec1 != NULL && spec2 == NULL) ||
        (spec1 == NULL && spec2 != NULL) ||
        spec1->n != spec2->n) {
        sprintf(errtxt, "%s: smspec *spec1 != smspec *spec2",
                SM_ERROR_IDG_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IDG, "%s", errtxt);
    }

    if (spec1->n > 0) {
        for (i = 0; i < spec1->n; i++) {
            llam = spec1->dat[i].lam * (1 - SM_TOL);
            ulam = spec1->dat[i].lam * (1 + SM_TOL);
            if (spec2->dat[i].lam < llam || spec2->dat[i].lam > ulam) {
                sprintf(errtxt, "%s: smspec *spec1 != smspec *spec2",
                        SM_ERROR_IDG_TXT);
                return cpl_error_set_message(cpl_func, SM_ERROR_IDG, "%s",
                                             errtxt);
            }
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_join(smspec *spec, const smspec *errfunc,
                            const int errflag)
{
    /*!
     * Uses the flux of a second ::smspec structure to fill the error vector
     * of the first spectrum. If two error columns are used, lower or upper
     * error has to be specified and the function has to be called for either
     * case (starting with the lower error).
     *
     * \note Pre-existing error data in spec will be overwritten.
     *
     * \b INPUT:
     * \param spec     spectrum without errors (2 columns)
     * \param errfunc  error function (2 columns)
     * \param errflag  error flag (0 = symmetric error, 1 = lower error,
     *                 2 = upper error)
     *
     * \b OUTPUT:
     * \param spec     spectrum with errors (3 columns)
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     */

    char errtxt[SM_MAXLEN+1];
    int i;

    /* Test agreement of wavelength grids of input spectra */
    if ((int) sm_spec_compgrids(spec, errfunc) == (int) SM_ERROR_IDG) {
        return SM_ERROR_IDG;
    }

    /* No data -> no copying */
    if (spec->n <= 0) {
        return CPL_ERROR_NONE;
    }

    /* Valid error flag? */
    if (errflag < 0 || errflag > 2) {
        sprintf(errtxt, "%s: errflag (flag for flux error type < 0 or > 2)",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Adapt spectrum type */
    if (errflag == 0) {
        spec->type = 3;
    } else {
        spec->type = 4;
    }

    /* Copy error function to spec */
    if (errflag == 2) {
        for (i = 0; i < spec->n; i++) {
            spec->dat[i].dflux2 = errfunc->dat[i].flux;
        }
    } else {
        for (i = 0; i < spec->n; i++) {
            spec->dat[i].dflux1 = errfunc->dat[i].flux;
            if (errflag == 0) {
                spec->dat[i].dflux2 = spec->dat[i].dflux1;
            }
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_split(smspec *spec, smspec *errfunc, const int errflag)
{
    /*!
     * Writes error function of an ::smspec structure to another spectrum.
     * The error vector of the former spectrum is filled with zero.
     * If two error columns are used, lower or upper error has to be
     * specified and the function has to be called for either case.
     *
     * \b INPUT:
     * \param spec     spectrum with errors (3 columns)
     *
     * \b OUTPUT:
     * \param spec     spectrum without errors (2 columns)
     * \param errfunc  error function (2 columns)
     * \param errflag  error flag (0 = symmetric error, 1 = lower error,
     *                 2 = upper error)
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     */

    char errtxt[SM_MAXLEN+1];
    int i;

    /* Set error status = false */
    spec->type = 2;

    /* Create smspec structure for error function */
    sm_spec_copy(errfunc, spec);

    /* No data -> no copying */
    if (spec->n <= 0) {
        return CPL_ERROR_NONE;
    }

    /* Valid error flag? */
    if (errflag < 0 || errflag > 2) {
        sprintf(errtxt, "%s: errflag (flag for flux error type < 0 or > 2)",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Move error vector of spec to flux vector of errfunc */
    if (errflag == 2) {
        for (i = 0; i < spec->n; i++) {
            errfunc->dat[i].flux = spec->dat[i].dflux2;
            spec->dat[i].dflux2 = 0;
        }
    } else {
        for (i = 0; i < spec->n; i++) {
            errfunc->dat[i].flux = spec->dat[i].dflux1;
            spec->dat[i].dflux1 = 0;
            if (errflag == 0) {
                spec->dat[i].dflux2 = 0;
            }
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_changetype(smspec *spec, const int type)
{
    /*!
     * Changes spectrum type (1 = no flux, 2 = no errors, 3 = symmetric error,
     * 4 = lower and upper error).
     * For changes to 1 or 2 the obsolete columns are set to zero.
     * For a change from 4 to 3 the lower and upper errors are made symmetric
     * by simple averaging.
     * For a change from 1-3 to 4 or 3 it is ensured that the lower and upper
     * errors are equal (using the lower errors as reference).
     * Only for a change from 4 to 4 nothing is done.
     *
     * \b INPUT:
     * \param spec  spectrum of arbitrary type
     * \param type  spectrum type (see ::smspec)
     *
     * \b OUTPUT:
     * \param spec  spectrum of given \e type
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     */

    char errtxt[SM_MAXLEN+1];
    int i = 0;

    /* Valid type? */
    if (type < 1 || type > 4) {
        sprintf(errtxt, "%s: type (spectrum type < 1 or > 4)",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Treat all cases */

    if (spec->type < 4 && type >= 3) {
        for (i = 0; i < spec->n; i++) {
            spec->dat[i].dflux2 = spec->dat[i].dflux1;
        }
    } else if (spec->type == 4 && type == 3) {
        for (i = 0; i < spec->n; i++) {
            spec->dat[i].dflux1 = (spec->dat[i].dflux1 +
                                   spec->dat[i].dflux2) / 2;
            spec->dat[i].dflux2 = spec->dat[i].dflux1;
        }
    } else if (type <= 2) {
        for (i = 0; i < spec->n; i++) {
            spec->dat[i].dflux1 = 0;
            spec->dat[i].dflux2 = 0;
            if (type == 1) {
                spec->dat[i].flux = 0;
            }
        }
    }

    /* Change type */
    spec->type = type;

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_getval(double *flux, const smspec *spec,
                              const double lam)
{
    /*!
     * Provides the interpolated flux for a given wavelength. If the input
     * wavelength is outside the range of the input spectrum, the flux of the
     * first or last spectrum wavelength is returned.
     *
     * \b INPUT:
     * \param spec  input spectrum
     * \param lam   input wavelength
     *
     * \b OUTPUT:
     * \param flux  flux for given wavelength
     *
     * \b ERRORS:
     * - NDA: No data
     */

    char errtxt[SM_MAXLEN+1];
    int i = 0, i0 = -1;
    double frac = 0;

    /* Default value */
    *flux = 0;

    /* No data */
    if (spec->n <= 0) {
        sprintf(errtxt, "%s: smspec *spec", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    }

    /* Find index for first larger wavelength */

    for (i = 0; i < spec->n; i++) {
        if (spec->dat[i].lam > lam) {
            i0 = i;
            break;
        }
    }

    if (i0 == -1) {
        i0 = spec->n;
    }

    /* Get pixel fraction */
    if (i0 != 0 && i0 != spec->n) {
        frac = (lam - spec->dat[i0].lam) /
               (spec->dat[i0-1].lam - spec->dat[i0].lam);
    }

    /* Calculate output flux */
    if (i0 == 0) {
        *flux = spec->dat[i0].flux;
    } else if (i0 == spec->n) {
        *flux = spec->dat[i0-1].flux;
    } else {
        *flux = frac * (spec->dat[i0-1].flux - spec->dat[i0].flux) +
                spec->dat[i0].flux;
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_scalerange(smspec *spec, const double limlam[2],
                                  const char op, const double c)
{
    /*!
     * Modifies the flux of an ::smspec structure in a given wavelength range
     * by adding, subtracting, multiplying, or dividing a constant. Moreover,
     * the flux can be raised to a constant power or the selected wavelength
     * range can be given a constant flux value.
     *
     * The two subtraction operators '-' and '_' differ by the order of
     * spectrum and constant, i.e., '-' and '_' imply \f$x - c\f$ and
     * \f$c - x\f$, respectively.
     * The two division operators '/' and '|' differ by the order of
     * spectrum and constant, i.e., '/' and '|' imply \f$x / c\f$ and
     * \f$c / x\f$, respectively.
     *
     * \b INPUT:
     * \param spec    input ::smspec structure
     * \param limlam  lower and upper wavelength limit for operation
     *                (Take HUGE_VAL for unconstrained limits.)
     * \param op      operation: '+', '-', '_', '*', '/', '|', '^', or '='
     * \param c       constant (double precision float)
     *
     * \b OUTPUT:
     * \param spec    modified spectrum
     *                (no modification in the case of errors)
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     * - CPL: Division by zero [warning]
     * - CPL: Access out of range [warning]
     */

    cpl_error_code err = CPL_ERROR_NONE;
    smbool isoutrange = F, isdiv0 = F;
    char errtxt[SM_MAXLEN+1];
    int i = 0;
    double dflux1 = 0., dflux2 = 0.;

    /* Check wavelength limits */

    if (limlam[1] <= limlam[0]) {
        /* Return unmodified spectrum */
        sprintf(errtxt, "%s: limlam[1] <= limlam[0] (wavelength limits)",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* No data -> no calculation */

    if (spec->n <= 0) {
        return CPL_ERROR_NONE;
    }

    /* Modify errors */

    if (spec->type >= 3) {

        for (i = 0; i < spec->n; i++) {

            /* Errors have to be positive */
            if (spec->dat[i].dflux1 < 0 || spec->dat[i].dflux2 < 0) {
                spec->dat[i].dflux1 = HUGE_VAL;
                spec->dat[i].dflux2 = HUGE_VAL;
                isoutrange = T;
                continue;
            }

            if (spec->dat[i].lam >= limlam[0] &&
                spec->dat[i].lam <= limlam[1]) {

                if (op == '+') {
                    /* No change */
                } else if (op == '-') {
                    /* No change */
                } else if (op == '_') {
                    /* No change */
                } else if (op == '*') {
                    if (c < 0.) {
                        /* Avoid undefined expressions */
                        spec->dat[i].dflux1 = HUGE_VAL;
                        spec->dat[i].dflux2 = HUGE_VAL;
                        isoutrange = T;
                    } else {
                        spec->dat[i].dflux1 *= c;
                        spec->dat[i].dflux2 *= c;
                    }
                } else if (op == '/') {
                    if (c == 0.) {
                        /* Avoid division by zero */
                        spec->dat[i].dflux1 = HUGE_VAL;
                        spec->dat[i].dflux2 = HUGE_VAL;
                        isdiv0 = T;
                    } else if (c < 0.) {
                        /* Avoid undefined expressions */
                        spec->dat[i].dflux1 = HUGE_VAL;
                        spec->dat[i].dflux2 = HUGE_VAL;
                        isoutrange = T;
                    } else {
                        spec->dat[i].dflux1 /= c;
                        spec->dat[i].dflux2 /= c;
                    }
                } else if (op == '=') {
                    spec->dat[i].dflux1 = 0.;
                    spec->dat[i].dflux2 = 0.;
                }

                if (spec->type == 3) {

                    /* Case: symmetric errors */

                    if (op == '|') {
                        if (spec->dat[i].flux == 0.) {
                            /* Avoid division by zero */
                            spec->dat[i].dflux1 = HUGE_VAL;
                            isdiv0 = T;
                        } else if (c < 0.) {
                            /* Avoid undefined expressions */
                            spec->dat[i].dflux1 = HUGE_VAL;
                            isoutrange = T;
                        } else {
                            spec->dat[i].dflux1 *= c /
                              pow(spec->dat[i].flux, 2);
                        }
                    } else if (op == '^') {
                        if (spec->dat[i].flux == 0. && c < 0.) {
                            /* Avoid division by zero */
                            spec->dat[i].dflux1 = HUGE_VAL;
                            isdiv0 = T;
                        } else {
                            spec->dat[i].dflux1 *=
                              c * pow(spec->dat[i].flux, c - 1);
                        }
                    }

                    spec->dat[i].dflux2 = spec->dat[i].dflux1;

                } else if (spec->type == 4) {

                    /* Case: asymmetric errors */

                    if (op == '|') {
                        if (spec->dat[i].flux == 0. ||
                            spec->dat[i].dflux1 == spec->dat[i].flux) {
                            /* Avoid division by zero */
                            dflux2 = HUGE_VAL;
                            isdiv0 = T;
                        } else if (spec->dat[i].dflux1 > spec->dat[i].flux ||
                                   c < 0.) {
                            /* Avoid undefined expressions */
                            dflux2 = HUGE_VAL;
                            isoutrange = T;
                        } else {
                            dflux2 = c * spec->dat[i].dflux1 /
                              (spec->dat[i].flux *
                               (spec->dat[i].flux - spec->dat[i].dflux1));
                        }
                        if (spec->dat[i].flux == 0.) {
                            /* Avoid division by zero */
                            dflux1 = HUGE_VAL;
                            isdiv0 = T;
                        } else if (c < 0.) {
                            /* Avoid undefined expressions */
                            dflux1 = HUGE_VAL;
                            isoutrange = T;
                        } else {
                            dflux1 = c * spec->dat[i].dflux2 /
                              (spec->dat[i].flux *
                               (spec->dat[i].flux + spec->dat[i].dflux2));
                        }
                        spec->dat[i].dflux1 = dflux1;
                        spec->dat[i].dflux2 = dflux2;
                    } else if (op == '^') {
                        if ((spec->dat[i].flux == 0. ||
                             spec->dat[i].dflux1 == spec->dat[i].flux) &&
                            c < 0.) {
                            /* Avoid division by zero */
                            dflux1 = -HUGE_VAL;
                            isdiv0 = T;
                        } else if (spec->dat[i].dflux1 > spec->dat[i].flux) {
                            /* Avoid undefined expressions */
                            dflux1 = HUGE_VAL;
                            isoutrange = T;
                        } else {
                            dflux1 =
                              pow(spec->dat[i].flux, c) -
                              pow(spec->dat[i].flux - spec->dat[i].dflux1, c);
                        }
                        if (spec->dat[i].flux == 0. && c < 0.) {
                            /* Avoid division by zero */
                            dflux2 = -HUGE_VAL;
                            isdiv0 = T;
                        } else {
                            dflux2 =
                              pow(spec->dat[i].flux + spec->dat[i].dflux2, c)-
                              pow(spec->dat[i].flux, c);
                        }
                        if (dflux1 < 0 && dflux2 < 0) {
                            spec->dat[i].dflux1 = -dflux2;
                            spec->dat[i].dflux2 = -dflux1;
                        } else if (dflux1 >= 0 && dflux2 >= 0) {
                            spec->dat[i].dflux1 = dflux1;
                            spec->dat[i].dflux2 = dflux2;
                        } else {
                            /* Avoid undefined expressions */
                            spec->dat[i].dflux1 = HUGE_VAL;
                            spec->dat[i].dflux2 = HUGE_VAL;
                            isoutrange = T;
                        }
                    }

                }

            }

        }

    }

    /* Perform operation */

    for (i = 0; i < spec->n; i++) {
        if (spec->dat[i].lam >= limlam[0] && spec->dat[i].lam <= limlam[1]) {
            if (op == '+') {
                spec->dat[i].flux += c;
            } else if (op == '-') {
                spec->dat[i].flux -= c;
            } else if (op == '_') {
                spec->dat[i].flux = c - spec->dat[i].flux;
            } else if (op == '*') {
                spec->dat[i].flux *= c;
            } else if (op == '/') {
                if (c == 0.) {
                    /* Avoid division by zero */
                    spec->dat[i].flux = 0.;
                    isdiv0 = T;
                } else {
                    spec->dat[i].flux /= c;
                }
            } else if (op == '|') {
                if (spec->dat[i].flux == 0.) {
                    /* Avoid division by zero */
                    spec->dat[i].flux = 0.;
                    isdiv0 = T;
                } else {
                    spec->dat[i].flux = c / spec->dat[i].flux;
                }
            } else if (op == '^') {
                if (spec->dat[i].flux == 0. && c < 0.) {
                    /* Avoid division by zero */
                    spec->dat[i].flux = 0.;
                    isdiv0 = T;
                } else if (spec->dat[i].flux < 0. &&
                           (c != floor(c) && c != ceil(c))) {
                    /* Avoid undefined expressions */
                    spec->dat[i].flux = 0.;
                    isoutrange = T;
                } else {
                    spec->dat[i].flux = pow(spec->dat[i].flux, c);
                }
            } else if (op == '=') {
                spec->dat[i].flux = c;
            } else {
                sprintf(errtxt, "%s: %c (operators '+', '-', '_', '*', '/', "
                        "'|', '^', or '=' only)", SM_ERROR_IIP_TXT, op);
                return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s",
                                             errtxt);
            }
        }
    }

    if (isdiv0 == T) {
        sprintf(errtxt, "smspec *spec, operator %c, constant %g "
                "[warning]", op, c);
        err = cpl_error_set_message(cpl_func, CPL_ERROR_DIVISION_BY_ZERO,
                                    "%s", errtxt);
    }

    if (isoutrange == T) {
        sprintf(errtxt, "smspec *spec, operator %c, constant %g "
                "(undefined expression(s)) [warning]", op, c);
        err = cpl_error_set_message(cpl_func, CPL_ERROR_ACCESS_OUT_OF_RANGE,
                                    "%s", errtxt);
    }

    return err;
}


cpl_error_code sm_spec_scale(smspec *spec, const char op, const double c)
{
    /*!
     * Modifies the flux of an ::smspec structure by adding, subtracting,
     * multiplying, or dividing a constant. Moreover, the flux can be raised
     * to a constant power or the full spectrum can be given a constant flux
     * value.
     *
     * The two subtraction operators '-' and '_' differ by the order of
     * spectrum and constant, i.e., '-' and '_' imply \f$x - c\f$ and
     * \f$c - x\f$, respectively.
     * The two division operators '/' and '|' differ by the order of
     * spectrum and constant, i.e., '/' and '|' imply \f$x / c\f$ and
     * \f$c / x\f$, respectively.
     *
     * \b INPUT:
     * \param spec  input ::smspec structure
     * \param op    operation: '+', '-', '*', '/', '^', '=', '_', or '|'
     * \param c     constant (double precision float)
     *
     * \b OUTPUT:
     * \param spec  modified spectrum
     *              (no modification in the case of errors)
     *
     * \b ERRORS:
     * - see ::sm_spec_scalerange
     */

    double limlam[2];

    limlam[0] = -HUGE_VAL;
    limlam[1] = HUGE_VAL;

    sm_spec_scalerange(spec, limlam, op, c);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_modval(smspec *spec, const double lam, const char op,
                              const double c)
{
    /*!
     * Modifies the flux of a given wavelength belonging to an ::smspec
     * structure by adding, subtracting, multiplying, or dividing a constant.
     * Moreover, the flux can be raised to a constant power or the flux of the
     * selected wavelength can be substituted.
     *
     * The two subtraction operators '-' and '_' differ by the order of
     * spectrum and constant, i.e., '-' and '_' imply \f$x - c\f$ and
     * \f$c - x\f$, respectively.
     * The two division operators '/' and '|' differ by the order of
     * spectrum and constant, i.e., '/' and '|' imply \f$x / c\f$ and
     * \f$c / x\f$, respectively.
     *
     * \b INPUT:
     * \param spec  input ::smspec structure
     * \param lam   wavelength (accuracy ruled by ::SM_TOL)
     * \param op    operation: '+', '-', '*', '/', '^', '=', '_', or '|'
     * \param c     constant (double precision float)
     *
     * \b OUTPUT:
     * \param spec  modified spectrum
     *              (no modification in the case of errors)
     *
     * \b ERRORS:
     * - see ::sm_spec_scalerange
     */

    double limlam[2];

    limlam[0] = lam * (1 - SM_TOL);
    limlam[1] = lam * (1 + SM_TOL);

    sm_spec_scalerange(spec, limlam, op, c);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_calc(smspec *spec, const char op, const smspec *opspec)
{
    /*!
     * Performs addition, subtraction, multiplication, division, or
     * equalisation of two spectra with error functions. Moreover, a spectrum
     * can be used as a power of another spectrum.
     * Error propagation (adding of squared errors) is considered.
     * For asymmetric errors positive fluxes are required for multiplicative
     * operations.
     *
     * The two subtraction operators '-' and '_' differ by the order of
     * the spectra, i.e., '-' and '_' imply \f$x_1 - x_2\f$ and
     * \f$x_2 - x_1\f$, respectively.
     * The two division operators '/' and '|' differ by the order of
     * the spectra, i.e., '/' and '|' imply \f$x_1 / x_2\f$ and
     * \f$x_2 / x_1\f$, respectively.
     *
     * \b INPUT:
     * \param spec    input ::smspec structure
     * \param op      operation: '+', '-', '_', '*', '/', '|', '^', or '='
     * \param opspec  ::smspec structure to modify spec
     *                (identical wavelengths needed!)
     *
     * \b OUTPUT:
     * \param spec    modified spectrum
     *                (no modification in the case of errors)
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     * - CPL: Division by zero [warning]
     * - CPL: Access out of range [warning]
     */

    cpl_error_code err = CPL_ERROR_NONE;
    smbool isdiv0 = F, isoutrange = F;
    char errtxt[SM_MAXLEN+1];
    int i = 0;
    double flux = 0., dflux1 = 0., dflux2 = 0., dopflux = 0., term = 0.;

    /* Test agreement of wavelength grids of input spectra */

    if ((int) sm_spec_compgrids(spec, opspec) == (int) SM_ERROR_IDG) {
        return SM_ERROR_IDG;
    }

    /* No data -> no calculation */

    if (spec->n <= 0) {
        return CPL_ERROR_NONE;
    }

    /* Perform operation */

    for (i = 0; i < spec->n; i++) {

        /* Calculate output spectrum */

        if (op == '+') {
            flux = spec->dat[i].flux + opspec->dat[i].flux;
        } else if (op == '-') {
            flux = spec->dat[i].flux - opspec->dat[i].flux;
        } else if (op == '_') {
            flux = opspec->dat[i].flux - spec->dat[i].flux;
        } else if (op == '*') {
            flux = spec->dat[i].flux * opspec->dat[i].flux;
        } else if (op == '/') {
            if (opspec->dat[i].flux == 0.) {
                /* Avoid division by zero */
                flux = 0.;
                isdiv0 = T;
            } else {
                flux = spec->dat[i].flux / opspec->dat[i].flux;
            }
        } else if (op == '|') {
            if (spec->dat[i].flux == 0.) {
                /* Avoid division by zero */
                flux = 0.;
                isdiv0 = T;
            } else {
                flux = opspec->dat[i].flux / spec->dat[i].flux;
            }
        } else if (op == '^') {
            if (spec->dat[i].flux == 0. && opspec->dat[i].flux < 0.) {
                /* Avoid division by zero */
                flux = 0.;
                isdiv0 = T;
            } else {
                flux = pow(spec->dat[i].flux, opspec->dat[i].flux);
            }
        } else if (op == '=') {
            flux = opspec->dat[i].flux;
        } else {
            sprintf(errtxt, "%s: %c (operators '+', '-', '_', '*', '/', "
                    "'|', '^', or '=' only)", SM_ERROR_IIP_TXT, op);
            return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s",
                                         errtxt);
        }

        /* Errors have to be positive */
        if (spec->type >= 3 &&
            (spec->dat[i].dflux1 < 0 || spec->dat[i].dflux2 < 0 ||
             opspec->dat[i].dflux1 < 0 || opspec->dat[i].dflux2 < 0)) {
            dflux1 = HUGE_VAL;
            dflux2 = HUGE_VAL;
            isoutrange = T;
            continue;
        }

        /* Calculate output error functions */

        if (spec->type == 3) {

            /* Case: symmetric errors for target input spectrum */

            if (opspec->type == 4) {
                /* Make errors of 2nd spectrum symmetric if necessary */
                dopflux = (opspec->dat[i].dflux1 + opspec->dat[i].dflux2) / 2;
            } else {
                dopflux = opspec->dat[i].dflux1;
            }

            if (op == '+' || op == '-' || op == '_') {
                dflux1 = sqrt(pow(spec->dat[i].dflux1, 2) +
                              pow(dopflux, 2));
            } else if (op == '*' ) {
                dflux1 = sqrt(pow(opspec->dat[i].flux *
                                  spec->dat[i].dflux1, 2) +
                              pow(spec->dat[i].flux *
                                  dopflux, 2));
            } else if (op == '/') {
                if (opspec->dat[i].flux == 0.) {
                    /* Avoid division by zero */
                    dflux1 = HUGE_VAL;
                    isdiv0 = T;
                } else {
                    dflux1 = sqrt(pow(spec->dat[i].dflux1, 2) +
                                  pow(flux * dopflux, 2)) /
                             fabs(opspec->dat[i].flux);
                }
            } else if (op == '|') {
                if (spec->dat[i].flux == 0.) {
                    /* Avoid division by zero */
                    dflux1 = HUGE_VAL;
                    isdiv0 = T;
                } else {
                    dflux1 = sqrt(pow(dopflux, 2) +
                                  pow(flux * spec->dat[i].dflux1, 2)) /
                             fabs(spec->dat[i].flux);
                }
            } else if (op == '^') {
                if (spec->dat[i].flux <= 0.) {
                    /* Avoid undefined expressions */
                    dflux1 = HUGE_VAL;
                    isoutrange = T;
                } else {
                    dflux1 = flux *
                      sqrt(pow((opspec->dat[i].flux / spec->dat[i].flux) *
                               spec->dat[i].dflux1, 2) +
                           pow(log(spec->dat[i].flux) * dopflux, 2));
                }
            } else if (op == '=') {
                dflux1 = dopflux;
            }

            dflux2 = dflux1;

        }

        if (spec->type == 4) {

            /* Case: asymmetric errors for target input spectrum */

            if (op == '+') {
                dflux1 = sqrt(pow(spec->dat[i].dflux1, 2) +
                              pow(opspec->dat[i].dflux1, 2));
                dflux2 = sqrt(pow(spec->dat[i].dflux2, 2) +
                              pow(opspec->dat[i].dflux2, 2));
            } else if (op == '-' ) {
                dflux1 = sqrt(pow(spec->dat[i].dflux1, 2) +
                              pow(opspec->dat[i].dflux2, 2));
                dflux2 = sqrt(pow(spec->dat[i].dflux2, 2) +
                              pow(opspec->dat[i].dflux1, 2));
            } else if (op == '_' ) {
                dflux1 = sqrt(pow(spec->dat[i].dflux2, 2) +
                              pow(opspec->dat[i].dflux1, 2));
                dflux2 = sqrt(pow(spec->dat[i].dflux1, 2) +
                              pow(opspec->dat[i].dflux2, 2));
            } else if (op == '*' ) {
                dflux1 = sqrt(pow(opspec->dat[i].flux *
                                  spec->dat[i].dflux1, 2) +
                              pow(spec->dat[i].flux *
                                  opspec->dat[i].dflux1, 2));
                dflux2 = sqrt(pow(opspec->dat[i].flux *
                                  spec->dat[i].dflux2, 2) +
                              pow(spec->dat[i].flux *
                                  opspec->dat[i].dflux2, 2));
            } else if (op == '/') {
                if (opspec->dat[i].flux == 0.) {
                    /* Avoid division by zero */
                    dflux1 = HUGE_VAL;
                    isdiv0 = T;
                } else if (spec->dat[i].flux < 0 || opspec->dat[i].flux < 0) {
                    /* Avoid undefined expressions */
                    dflux1 = HUGE_VAL;
                    isoutrange = T;
                } else {
                    dflux1 = sqrt(pow(spec->dat[i].dflux1, 2) +
                                  pow(spec->dat[i].flux *
                                      opspec->dat[i].dflux2 /
                                      (opspec->dat[i].flux +
                                       opspec->dat[i].dflux2), 2)) /
                             opspec->dat[i].flux;
                }
                if (opspec->dat[i].flux == 0. ||
                    opspec->dat[i].flux == opspec->dat[i].dflux1) {
                    /* Avoid division by zero */
                    dflux2 = HUGE_VAL;
                    isdiv0 = T;
                } else if (spec->dat[i].flux < 0 ||
                           opspec->dat[i].flux < opspec->dat[i].dflux1) {
                    /* Avoid undefined expressions */
                    dflux2 = HUGE_VAL;
                    isoutrange = T;
                } else {
                    dflux2 = sqrt(pow(spec->dat[i].dflux2, 2) +
                                  pow(spec->dat[i].flux *
                                      opspec->dat[i].dflux1 /
                                      (opspec->dat[i].flux -
                                       opspec->dat[i].dflux1), 2)) /
                             opspec->dat[i].flux;
                }
            } else if (op == '|') {
                if (spec->dat[i].flux == 0.) {
                    /* Avoid division by zero */
                    dflux1 = HUGE_VAL;
                    isdiv0 = T;
                } else if (spec->dat[i].flux < 0 || opspec->dat[i].flux < 0) {
                    /* Avoid undefined expressions */
                    dflux1 = HUGE_VAL;
                    isoutrange = T;
                } else {
                    dflux1 = sqrt(pow(opspec->dat[i].dflux1, 2) +
                                  pow(opspec->dat[i].flux *
                                      spec->dat[i].dflux2 /
                                      (spec->dat[i].flux +
                                       spec->dat[i].dflux2), 2)) /
                             spec->dat[i].flux;
                }
                if (spec->dat[i].flux == 0. ||
                    spec->dat[i].flux == spec->dat[i].dflux1) {
                    /* Avoid division by zero */
                    dflux2 = HUGE_VAL;
                    isdiv0 = T;
                } else if (spec->dat[i].flux < spec->dat[i].dflux1 ||
                           opspec->dat[i].flux < 0) {
                    /* Avoid undefined expressions */
                    dflux2 = HUGE_VAL;
                    isoutrange = T;
                } else {
                    dflux2 = sqrt(pow(opspec->dat[i].dflux2, 2) +
                                  pow(opspec->dat[i].flux *
                                      spec->dat[i].dflux1 /
                                      (spec->dat[i].flux -
                                       spec->dat[i].dflux1), 2)) /
                             spec->dat[i].flux;
                }
            } else if (op == '^') {
                if (spec->dat[i].flux == 0. &&
                    opspec->dat[i].flux < 0) {
                    /* Avoid division by zero */
                    dflux1 = HUGE_VAL;
                    isdiv0 = T;
                } else if (spec->dat[i].flux < spec->dat[i].dflux1 &&
                           opspec->dat[i].flux > 0) {
                    /* Avoid undefined expressions */
                    dflux1 = HUGE_VAL;
                    isoutrange = T;
                } else {
                    term = pow(flux -
                               pow(spec->dat[i].flux,
                                   opspec->dat[i].flux -
                                   opspec->dat[i].dflux1), 2);
                    if (opspec->dat[i].flux < 0) {
                        dflux1 = sqrt(pow(pow(spec->dat[i].flux +
                                              spec->dat[i].dflux2,
                                              opspec->dat[i].flux) - flux, 2)+
                                      term);
                    } else {
                        dflux1 = sqrt(pow(flux -
                                          pow(spec->dat[i].flux -
                                              spec->dat[i].dflux1,
                                              opspec->dat[i].flux), 2) +
                                      term);
                    }
                }
                if (spec->dat[i].flux == 0. &&
                    opspec->dat[i].flux < 0) {
                    /* Avoid division by zero */
                    dflux2 = HUGE_VAL;
                    isdiv0 = T;
                } else if (spec->dat[i].flux < spec->dat[i].dflux1 &&
                           opspec->dat[i].flux < 0) {
                    /* Avoid undefined expressions */
                    dflux2 = HUGE_VAL;
                    isoutrange = T;
                } else {
                    term = pow(pow(spec->dat[i].flux,
                                   opspec->dat[i].flux +
                                   opspec->dat[i].dflux2) - flux, 2);
                    if (opspec->dat[i].flux < 0) {
                        dflux2 = sqrt(pow(flux -
                                          pow(spec->dat[i].flux -
                                              spec->dat[i].dflux1,
                                              opspec->dat[i].flux), 2) +
                                      term);
                    } else {
                        dflux2 = sqrt(pow(pow(spec->dat[i].flux +
                                              spec->dat[i].dflux2,
                                              opspec->dat[i].flux) - flux, 2)+
                                      term);
                    }
                }
           } else if (op == '=' ) {
                dflux1 = opspec->dat[i].dflux1;
                dflux2 = opspec->dat[i].dflux2;
            }

        }

        spec->dat[i].flux = flux;
        if (spec->type >= 3) {
            spec->dat[i].dflux1 = dflux1;
            spec->dat[i].dflux2 = dflux2;
        }

    }

    if (isdiv0 == T) {
        sprintf(errtxt, "smspec *spec, operator %c, smspec *opspec "
                "[warning]", op);
        err = cpl_error_set_message(cpl_func, CPL_ERROR_DIVISION_BY_ZERO,
                                    "%s", errtxt);
    }

    if (isoutrange == T) {
        sprintf(errtxt, "smspec *spec, operator %c, smspec *opspec "
                "(undefined expression(s)) [warning]", op);
        err = cpl_error_set_message(cpl_func, CPL_ERROR_ACCESS_OUT_OF_RANGE,
                                    "%s", errtxt);
    }

    return err;
}


cpl_error_code sm_spec_funct(smspec *spec, const char *funct,
                             const char baselab)
{
    /*!
     * Applies exponential or logarithmic function. Names of functions equal
     * those of <math.h> excepting "exp10", "expmag", and "logmag" which
     * correspond to \f$10^x\f$, \f$10^{-0.4\,x}\f$, and
     * \f$-2.5\,\log_{10}(x)\f$, respectively.
     * Correction of resulting error asymmetry by means of the recipe
     * \f$\frac{1}{2} (f(x + {\rm d}x) - f(x - {\rm d}x))\f$ for spectrum
     * type 3 (see ::smspec).
     *
     * \b INPUT:
     * \param spec     input ::smspec structure
     * \param funct    "exp" or "log"
     * \param baselab  'e', 'd' (10), or 'm' (\f$10^{-0.4}\f$)
     *
     * \b OUTPUT:
     * \param spec     modified spectrum if function is valid
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     * - CPL: Access out of range [warning]
     */

    smbool isoutrange = F;
    char errtxt[SM_MAXLEN+1];
    int i;
    double base = 0., dflux1 = 0., dflux2 = 0., relerr = 0.;

    /* No data -> no calculation */

    if (spec->n <= 0) {
        return CPL_ERROR_NONE;
    }

    /* Derive base */

    if (baselab == 'e') {
        base = exp(1);
    } else if (baselab == 'd') {
        base = 10;
    } else if (baselab == 'm') {
        base = pow(10, -0.4);
    } else {
        sprintf(errtxt, "%s: %c ('e', 'd', or 'm' for base label only)",
                SM_ERROR_IIP_TXT, baselab);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Apply function */

    for (i = 0; i < spec->n; i++) {

        /* Errors have to be positive */
        if (spec->type >= 3 &&
            (spec->dat[i].dflux1 < 0 || spec->dat[i].dflux2 < 0)) {
            spec->dat[i].dflux1 = HUGE_VAL;
            spec->dat[i].dflux2 = HUGE_VAL;
            isoutrange = T;
            continue;
        }

        if (strncmp(funct, "exp", strlen(funct)) == 0) {
            spec->dat[i].flux = pow(base, spec->dat[i].flux);
            if (spec->type == 3) {
                /* Average error (correction of error asymmetry) */
                spec->dat[i].dflux1 = sinh(fabs(log(base)) *
                                           spec->dat[i].dflux1) *
                                      spec->dat[i].flux;
                spec->dat[i].dflux2 = spec->dat[i].dflux1;
            } else if (spec->type == 4) {
                /* Exact lower and upper errors */
                if (base < 1.) {
                    dflux1 = (1 - pow(base, spec->dat[i].dflux2)) *
                             spec->dat[i].flux;
                    dflux2 = (pow(base, - spec->dat[i].dflux1) - 1) *
                             spec->dat[i].flux;
                } else {
                    dflux1 = (1 - pow(base, - spec->dat[i].dflux1)) *
                             spec->dat[i].flux;
                    dflux2 = (pow(base, spec->dat[i].dflux2) - 1) *
                             spec->dat[i].flux;
                }
                spec->dat[i].dflux1 = dflux1;
                spec->dat[i].dflux2 = dflux2;
            }
        } else if (strncmp(funct, "log", strlen(funct)) == 0) {
            if (spec->dat[i].flux <= 0.) {
                spec->dat[i].flux = 0.;
                if (spec->type == 3) {
                    spec->dat[i].dflux1 = HUGE_VAL;
                    spec->dat[i].dflux2 = HUGE_VAL;
                }
                isoutrange = T;
            } else {
                if (spec->type == 3) {
                    /* Average error (correction of error asymmetry) */
                    relerr = spec->dat[i].dflux1 / spec->dat[i].flux;
                    spec->dat[i].dflux1 = log(relerr +
                                             sqrt(pow(relerr, 2) + 1)) /
                                             fabs(log(base));
                    spec->dat[i].dflux2 = spec->dat[i].dflux1;
                } else if (spec->type == 4) {
                    /* Exact lower and upper errors */
                    if (base < 1.) {
                        dflux1 = (log(spec->dat[i].flux) -
                                  log(spec->dat[i].flux +
                                      spec->dat[i].dflux2)) / log(base);
                        dflux2 = (log(spec->dat[i].flux -
                                      spec->dat[i].dflux1) -
                                  log(spec->dat[i].flux)) / log(base);
                    } else {
                        dflux1 = (log(spec->dat[i].flux) -
                                  log(spec->dat[i].flux -
                                      spec->dat[i].dflux1)) / log(base);
                        dflux2 = (log(spec->dat[i].flux +
                                      spec->dat[i].dflux2) -
                                  log(spec->dat[i].flux)) / log(base);
                    }
                    spec->dat[i].dflux1 = dflux1;
                    spec->dat[i].dflux2 = dflux2;
                }
                spec->dat[i].flux = log(spec->dat[i].flux) / log(base);
            }
        } else {
            sprintf(errtxt, "%s: %s (function 'exp' or 'log' only)",
                    SM_ERROR_IIP_TXT, funct);
            return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s",
                                         errtxt);
        }

    }

    if (isoutrange == T) {
        sprintf(errtxt, "smspec *spec, function %s, base label %c"
                "(undefined expression(s)) [warning]", funct, baselab);
        return cpl_error_set_message(cpl_func, CPL_ERROR_ACCESS_OUT_OF_RANGE,
                                     "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_functnoerr(smspec *spec, const char *funct)
{
    /*!
     * Applies mathematical function to flux. A possible error function is not
     * changed.
     * Names of functions equal those of <math.h> excepting "exp10", "expmag",
     * and "logmag" which correspond to \f$10^x\f$, \f$10^{-0.4\,x}\f$, and
     * \f$-2.5\,\log_{10}(x)\f$, respectively.
     *
     * \b INPUT:
     * \param spec   input ::smspec structure
     * \param funct  "acos", "asin", "atan", "cos", "sin", "tan", "cosh",
     *               "sinh", "tanh", "exp", "exp10", "expmag", "log", "log10",
     *               "logmag", "sqrt", "ceil", "fabs", "floor"
     *
     * \b OUTPUT:
     * \param spec   modified spectrum if function is valid
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     * - CPL: Access out of range [warning]
     */

    smbool isoutrange = F;
    char errtxt[SM_MAXLEN+1];
    int i;

    /* No data -> no calculation */

    if (spec->n <= 0) {
        return CPL_ERROR_NONE;
    }

    /* Apply function */

    for (i = 0; i < spec->n; i++) {
        if (strncmp(funct, "acos", strlen(funct)) == 0) {
            if (spec->dat[i].flux < -1 || spec->dat[i].flux > 1) {
                spec->dat[i].flux = 0.;
                isoutrange = T;
            } else {
                spec->dat[i].flux = acos(spec->dat[i].flux);
            }
        } else if (strncmp(funct, "asin", strlen(funct)) == 0) {
            if (spec->dat[i].flux < -1 || spec->dat[i].flux > 1) {
                spec->dat[i].flux = 0.;
                isoutrange = T;
            } else {
                spec->dat[i].flux = asin(spec->dat[i].flux);
            }
        } else if (strncmp(funct, "atan", strlen(funct)) == 0) {
            spec->dat[i].flux = atan(spec->dat[i].flux);
        } else if (strncmp(funct, "cos", strlen(funct)) == 0) {
            spec->dat[i].flux = cos(spec->dat[i].flux);
        } else if (strncmp(funct, "sin", strlen(funct)) == 0) {
            spec->dat[i].flux = sin(spec->dat[i].flux);
        } else if (strncmp(funct, "tan", strlen(funct)) == 0) {
            spec->dat[i].flux = tan(spec->dat[i].flux);
        } else if (strncmp(funct, "cosh", strlen(funct)) == 0) {
            spec->dat[i].flux = cosh(spec->dat[i].flux);
        } else if (strncmp(funct, "sinh", strlen(funct)) == 0) {
            spec->dat[i].flux = sinh(spec->dat[i].flux);
        } else if (strncmp(funct, "tanh", strlen(funct)) == 0) {
            spec->dat[i].flux = tanh(spec->dat[i].flux);
        } else if (strncmp(funct, "exp", strlen(funct)) == 0) {
            spec->dat[i].flux = exp(spec->dat[i].flux);
        } else if (strncmp(funct, "exp10", strlen(funct)) == 0) {
            spec->dat[i].flux = pow(10., spec->dat[i].flux);
        } else if (strncmp(funct, "expmag", strlen(funct)) == 0) {
            spec->dat[i].flux = pow(10., -0.4 * spec->dat[i].flux);
        } else if (strncmp(funct, "log", strlen(funct)) == 0) {
            if (spec->dat[i].flux <= 0.) {
                spec->dat[i].flux = 0.;
                isoutrange = T;
            } else {
                spec->dat[i].flux = log(spec->dat[i].flux);
            }
        } else if (strncmp(funct, "log10", strlen(funct)) == 0) {
            if (spec->dat[i].flux <= 0.) {
                spec->dat[i].flux = 0.;
                isoutrange = T;
            } else {
                spec->dat[i].flux = log10(spec->dat[i].flux);
            }
        } else if (strncmp(funct, "logmag", strlen(funct)) == 0) {
            if (spec->dat[i].flux <= 0.) {
                spec->dat[i].flux = 0.;
                isoutrange = T;
            } else {
                spec->dat[i].flux = -2.5 * log10(spec->dat[i].flux);
            }
        } else if (strncmp(funct, "sqrt", strlen(funct)) == 0) {
            if (spec->dat[i].flux < 0.) {
                spec->dat[i].flux = 0.;
                isoutrange = T;
            } else {
                spec->dat[i].flux = sqrt(spec->dat[i].flux);
            }
        } else if (strncmp(funct, "ceil", strlen(funct)) == 0) {
            spec->dat[i].flux = ceil(spec->dat[i].flux);
        } else if (strncmp(funct, "fabs", strlen(funct)) == 0) {
            spec->dat[i].flux = fabs(spec->dat[i].flux);
        } else if (strncmp(funct, "floor", strlen(funct)) == 0) {
            spec->dat[i].flux = floor(spec->dat[i].flux);
        } else {
            sprintf(errtxt, "%s: %s (math.h functions, 'exp10', 'expmag', "
                    "or 'logmag' only)", SM_ERROR_IIP_TXT, funct);
            return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s",
                                         errtxt);
        }
    }

    if (isoutrange == T) {
        sprintf(errtxt, "smspec *spec, function %s "
                "(undefined expression(s)) [warning]", funct);
        return cpl_error_set_message(cpl_func, CPL_ERROR_ACCESS_OUT_OF_RANGE,
                                     "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_convunits(smspec *spec, const double factor,
                                 const int lamexp)
{
    /*!
     * Modifies the flux of an ::smspec structure by multiplying a \e factor
     * times \f$\lambda^{lamexp}\f$.
     *
     * \b INPUT:
     * \param spec    input ::smspec structure
     * \param factor  constant factor
     * \param lamexp  power for wavelengths
     *
     * \b OUTPUT:
     * \param spec    modified spectrum
     *
     * \b ERRORS:
     * - IDR: Invalid data range
     */

    char errtxt[SM_MAXLEN+1];
    int i;

    /* No data -> no calculation */

    if (spec->n <= 0) {
        return CPL_ERROR_NONE;
    }

    /* Convert units */

    for (i = 0; i < spec->n; i++) {
        if (spec->dat[i].lam <= 0.) {
            spec->dat[i].flux = 0;
            if (spec->type >= 3) {
                spec->dat[i].dflux1 = 0;
                spec->dat[i].dflux2 = 0;
            }
            sprintf(errtxt, "%s: smspec *spec (wavelengths <= 0)",
                    SM_ERROR_IDR_TXT);
            return cpl_error_set_message(cpl_func, SM_ERROR_IDR, "%s",
                                         errtxt);
        } else {
            spec->dat[i].flux *= factor * pow(spec->dat[i].lam, lamexp);
            if (spec->type >= 3) {
                spec->dat[i].dflux1 *= factor * pow(spec->dat[i].lam, lamexp);
                if (spec->type == 3) {
                    spec->dat[i].dflux2 = spec->dat[i].dflux1;
                } else {
                    spec->dat[i].dflux2 *= factor *
                                           pow(spec->dat[i].lam, lamexp);
                }
            }
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_changegrid(smspec *spec, const double factor,
                                  const char *scale)
{
    /*!
     * Modifies the wavelengths of an ::smspec structure by multiplying a
     * factor and/or performing \f$\log_{10}(x)\f$ or \f$10^x\f$.
     *
     * \b INPUT:
     * \param spec    input ::smspec structure
     * \param factor  wavelengths are multiplied by this value (\f$c\f$).
     * \param scale   wavelength scale, options:
     *                - "log": logarithmic (\f$\log_{10}(c\,x)\f$)
     *                - "exp": exponential (\f$c\,10^x\f$)
     *                         \f$\to\f$ reverses "log" option
     *                - "lin": linear (\f$c\,x\f$)
     *
     * \b OUTPUT:
     * \param spec    spectrum with modified wavelength scale
     *                (no modification in the case of errors)
     *
     * \b ERRORS:
     * - IDR: Invalid data range
     * - IIP: Invalid input parameter(s)
     */

    char errtxt[SM_MAXLEN+1];
    int i;

    /* No data -> no calculation */

    if (spec->n <= 0) {
        return CPL_ERROR_NONE;
    }

    /* Change wavelength grid */

    for (i = 0; i < spec->n; i++) {
        if (strncmp(scale, "log", strlen(scale)) == 0) {
            if (spec->dat[i].lam <= 0.) {
                sprintf(errtxt, "%s: smspec *spec (wavelengths <= 0)",
                        SM_ERROR_IDR_TXT);
                return cpl_error_set_message(cpl_func, SM_ERROR_IDR, "%s",
                                             errtxt);
            } else {
                spec->dat[i].lam = log10(factor * spec->dat[i].lam);
            }
        } else if (strncmp(scale, "exp", strlen(scale)) == 0) {
            spec->dat[i].lam = factor * pow(10., spec->dat[i].lam);
        } else if (strncmp(scale, "lin", strlen(scale)) == 0) {
            spec->dat[i].lam *= factor;
        } else {
            sprintf(errtxt, "%s: %s ('log', 'exp', or 'lin' for scale only)",
                    SM_ERROR_IIP_TXT, scale);
            return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s",
                                         errtxt);
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_minmax(double limflux[2], const smspec *spec,
                              const double limlam[2])
{
    /*!
     * Calculates minimum and maximum flux of a spectrum in the given
     * wavelength range
     *
     * \b INPUT:
     * \param spec     input spectrum
     * \param limlam   lower and upper wavelength limit for operation
     *                 (Take HUGE_VAL for unconstrained limits.)
     *
     * \b OUTPUT:
     * \param limflux  minimum and maximum flux
     *
     * \b ERRORS:
     * - NDA: No data
     * - IIP: Invalid input parameter(s)
     */

    char errtxt[SM_MAXLEN+1];
    int i;

    /* Default values */
    limflux[0] = HUGE_VAL;
    limflux[1] = -HUGE_VAL;

    /* No data */
    if (spec->n <= 0) {
        sprintf(errtxt, "%s: smspec *spec", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    }

    /* Check wavelength limits */
    if (limlam[1] <= limlam[0]) {
        sprintf(errtxt, "%s: limlam[1] <= limlam[0] (wavelength limits)",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Check wavelength ranges */
    if (spec->dat[0].lam > limlam[1] ||
        spec->dat[(spec->n)-1].lam < limlam[0]) {
        sprintf(errtxt, "%s: limlam[2] (no overlap with smspec *spec)",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Compute minimum and maximum flux for given wavelength range */
    for (i = 0; i < spec->n; i++) {
        if (spec->dat[i].lam >= limlam[0] && spec->dat[i].lam <= limlam[1]) {
            if (spec->dat[i].flux < limflux[0]) {
                limflux[0] = spec->dat[i].flux;
            }
            if (spec->dat[i].flux > limflux[1]) {
                limflux[1] = spec->dat[i].flux;
            }
       }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_average(double mean[3], const smspec *spec,
                               const double limlam[2])
{
    /*!
     * Calculates mean of a spectrum (and its error function if given) in
     * the given wavelength range
     *
     * \b INPUT:
     * \param spec    input spectrum
     * \param limlam  lower and upper wavelength limit for operation
     *                (Take HUGE_VAL for unconstrained limits.)
     *
     * \b OUTPUT:
     * \param mean    mean values for spectrum (and lower and upper error
     *                function if given)
     *
     * \b ERRORS:
     * - NDA: No data
     * - IIP: Invalid input parameter(s)
     */

    char errtxt[SM_MAXLEN+1];
    int i, n = 0;

    /* Default mean = 0 */
    mean[0] = 0;
    mean[1] = 0;
    mean[2] = 0;

    /* No data */
    if (spec->n <= 0) {
        sprintf(errtxt, "%s: smspec *spec", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    }

    /* Check wavelength limits */
    if (limlam[1] <= limlam[0]) {
        sprintf(errtxt, "%s: limlam[1] <= limlam[0] (wavelength limits)",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Check wavelength ranges */
    if (spec->dat[0].lam > limlam[1] ||
        spec->dat[(spec->n)-1].lam < limlam[0]) {
        sprintf(errtxt, "%s: limlam[2] (no overlap with smspec *spec)",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Compute average for given wavelength range */

    for (i = 0; i < spec->n; i++) {
        if (spec->dat[i].lam >= limlam[0] && spec->dat[i].lam <= limlam[1]) {
            n++;
            mean[0] = mean[0] + spec->dat[i].flux;
            if (spec->type >= 3) {
                mean[1] = mean[1] + spec->dat[i].dflux1;
            }
            if (spec->type == 4) {
                mean[2] = mean[2] + spec->dat[i].dflux2;
            }
        }
    }

    mean[0] = mean[0] / (double) n;
    if (spec->type >= 3) {
        mean[1] = mean[1] / (double) n;
        if (spec->type == 3) {
            /* Same average lower and upper error for spectrum type 3 */
            mean[2] = mean[1];
        } else {
            mean[2] = mean[2] / (double) n;
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_write(const smspec *spec, const char *filename)
{
    /*!
     * Writes ::smspec structure to ASCII file or on stdout
     * (indicated by "stdout" as filename)
     *
     * \b INPUT:
     * \param spec      input ::smspec structure
     * \param filename  name of output file or "stdout"
     *
     * \b ERRORS:
     * - FOF: File opening failed
     * - NDA: No data
     * - IOV: Invalid object value(s)
     */

    FILE *stream;
    char errtxt[SM_MAXLEN+1];
    int i;

    /* Output: file or stdout? */

    if (strncmp(filename, "stdout", 6) == 0) {
        /* Output on stdout */
        stream = stdout;
    } else {
        /* Check file existence */
        if ((stream = fopen(filename, "w")) == NULL) {
            sprintf(errtxt, "%s: %s", SM_ERROR_FOF_TXT, filename);
            return cpl_error_set_message(cpl_func, SM_ERROR_FOF, "%s",
                                         errtxt);
        }
    }

    /* Write header and spectral data */

    if (spec->n == 0) {
        sprintf(errtxt, "%s: smspec *spec", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    } else {
        switch (spec->type) {
        case 1:
            fprintf(stream, "# lam\n");
            for(i = 0; i < spec->n; i++) {
                fprintf(stream, "%e\n", spec->dat[i].lam);
            }
            break;
        case 2:
            fprintf(stream, "# lam        flux\n");
            for(i = 0; i < spec->n; i++) {
                fprintf(stream, "%e %e\n", spec->dat[i].lam,
                        spec->dat[i].flux);
            }
            break;
        case 3:
            fprintf(stream, "# lam        flux         dflux\n");
            for(i = 0; i < spec->n; i++) {
                fprintf(stream, "%e %e %e\n", spec->dat[i].lam,
                        spec->dat[i].flux, spec->dat[i].dflux1);
            }
            break;
        case 4:
            fprintf(stream,
                    "# lam        flux         dflux1       dflux2\n");
            for(i = 0; i < spec->n; i++) {
                fprintf(stream, "%e %e %e %e\n", spec->dat[i].lam,
                        spec->dat[i].flux, spec->dat[i].dflux1,
                        spec->dat[i].dflux2);
            }
            break;
        }
    }

    if (strncmp(filename, "stdout", 6) != 0) {
        fclose(stream);
    }

    if (spec->type < 1 || spec->type > 4) {
        sprintf(errtxt, "%s: smspec *spec (spectrum type < 1 or > 4)",
                SM_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IOV, "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_print(const smspec *spec)
{
    /*!
     * Prints ::smspec structure on stdout
     *
     * \b INPUT:
     * \param spec  input ::smspec structure
     *
     * \b ERRORS:
     * - none
     */

    sm_spec_write(spec, "stdout");

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_writecpl(cpl_table *cpltab, const smspec *spec)
{
    /*!
     * Writes ::smspec structure to CPL table
     *
     * \b INPUT:
     * \param spec    input ::smspec structure
     *
     * \b OUTPUT:
     * \param cpltab  output CPL table
     *
     * \b ERRORS:
     * - NDA: No data
     * - IOV: Invalid object value(s)
     */

    char errtxt[SM_MAXLEN+1];
    char colid[][7] = {"lam", "flux", "dflux", "dflux1", "dflux2"};
    char colid0[7];
    int i;
    double *lam, *flux, *dflux1, *dflux2;

    cpl_table_set_size(cpltab, spec->n);

    if (spec->n == 0) {
        sprintf(errtxt, "%s: cpl_table *cpltab", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    } else {
        if (spec->type >= 1 && spec->type <= 4) {
            cpl_table_new_column(cpltab, colid[0], CPL_TYPE_DOUBLE);
            cpl_table_fill_column_window_double(cpltab, colid[0], 0, spec->n,
                                                0.);
            lam = cpl_table_get_data_double(cpltab, colid[0]);
            for (i = 0; i < spec->n; i++) {
                lam[i] = spec->dat[i].lam;
            }
        }
        if (spec->type >= 2 && spec->type <= 4) {
            cpl_table_new_column(cpltab, colid[1], CPL_TYPE_DOUBLE);
            flux = cpl_table_get_data_double(cpltab, colid[1]);
            cpl_table_fill_column_window_double(cpltab, colid[1], 0, spec->n,
                                                0.);
            for (i = 0; i < spec->n; i++) {
                flux[i] = spec->dat[i].flux;
            }
        }
        if (spec->type >= 3 && spec->type <= 4) {
            if (spec->type == 3) {
                strcpy(colid0, colid[2]);
            } else {
                strcpy(colid0, colid[3]);
            }
            cpl_table_new_column(cpltab, colid0, CPL_TYPE_DOUBLE);
            cpl_table_fill_column_window_double(cpltab, colid0, 0, spec->n,
                                                0.);
            dflux1 = cpl_table_get_data_double(cpltab, colid0);
            for (i = 0; i < spec->n; i++) {
                dflux1[i] = spec->dat[i].dflux1;
            }
        }
        if (spec->type == 4) {
            cpl_table_new_column(cpltab, colid[4], CPL_TYPE_DOUBLE);
            cpl_table_fill_column_window_double(cpltab, colid[4], 0, spec->n,
                                                0.);
            dflux2 = cpl_table_get_data_double(cpltab, colid[4]);
            for (i = 0; i < spec->n; i++) {
                dflux2[i] = spec->dat[i].dflux2;
            }
        }
    }

    if (spec->type < 1 || spec->type > 4) {
        sprintf(errtxt, "%s: smspec *spec (spectrum type < 1 or > 4)",
                SM_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IOV, "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_writecplcolumn(cpl_table *cpltab, const char *colname,
                                      const smspec *spec, const int datatype)
{
    /*!
     * Writes data type (lam, flux, dflux1, or dflux2) of ::smspec structure
     * to column of CPL table. The column of the given name is created if it
     * does not exist. The CPL table size is increased if the number of data
     * points is higher in the ::smspec structure.
     *
     * \b INPUT:
     * \param cpltab    CPL table to be changed
     * \param colname   name of CPL table column for new data
     * \param spec      input ::smspec structure
     * \param datatype  type of data to be copied
     *                  (1 = lam, 2 = flux, 3 = dflux1, 4 = dflux2)
     *
     * \b OUTPUT:
     * \param cpltab    output CPL table with additional spectral data
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     */

    char errtxt[SM_MAXLEN+1];
    int n = 0, i = 0;
    double *col;

    /* Check data type */
    if (datatype < 1 || datatype > 4) {
        sprintf(errtxt, "%s: datatype < 1 or > 4", SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Return in the case of no data */
    if (spec->n == 0) {
        return CPL_ERROR_NONE;
    }

    /* Increase table size if input smspec structure has more rows */
    n = cpl_table_get_nrow(cpltab);
    if (spec->n > n) {
        cpl_table_set_size(cpltab, spec->n);
    }

    /* Create column in CPL table if necessary */
    if (cpl_table_has_column(cpltab, colname) != 1) {
        cpl_table_new_column(cpltab, colname, CPL_TYPE_DOUBLE);
    }

    /* Fill zero in targeted table elements */
    cpl_table_fill_column_window_double(cpltab, colname, 0, spec->n, 0.);

    /* Get pointer to CPL table column */
    col = cpl_table_get_data_double(cpltab, colname);

    /* Fill CPL table column with spectral data depending on data type */
    for (i = 0; i < spec->n; i++) {
        if (datatype == 1) {
            col[i] = spec->dat[i].lam;
        } else if (datatype == 2) {
            col[i] = spec->dat[i].flux;
        } else if (datatype == 3) {
            col[i] = spec->dat[i].dflux1;
        } else if (datatype == 4) {
            col[i] = spec->dat[i].dflux2;
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_writefits(const smspec *spec, const char *filename)
{
    /*!
     * Writes ::smspec structure to FITS file by means of CPL
     *
     * \b INPUT:
     * \param spec      input ::smspec structure
     * \param filename  name of output file
     *
     * \b ERRORS:
     * - none
     */

    cpl_table *cpltab;

    cpltab = cpl_table_new(0);
    sm_spec_writecpl(cpltab, spec);
    cpl_table_save(cpltab, NULL, NULL, filename, CPL_IO_CREATE);
    cpl_table_delete(cpltab);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_free(smspec *spec)
{
    /*!
     * Frees memory occupied by an ::smspec structure
     *
     * \b INPUT:
     * \param spec  ::smspec structure with \e n data points
     *
     * \b OUTPUT:
     * \param spec  ::smspec structure with zero data points
     *
     * \b ERRORS:
     * - none
     */

    if (spec->n > 0) {
        free(spec->dat);
        spec->dat = NULL;
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_extract(smspec *outspec, const smspec *inspec,
                               const double limlam[2])
{
    /*!
     * Extracts subspectrum in a given wavelength range of the input spectrum
     * (:: smspec structure)
     *
     * \b INPUT:
     * \param inspec   input spectrum
     * \param limlam   wavelength limits for extraction
     *
     * \b OUTPUT:
     * \param outspec  extracted spectrum
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     */

    char errtxt[SM_MAXLEN+1];
    int i = 0, pix1 = 0, count = 0;

    /* Copy type flag */
    outspec->type = inspec->type;

    /* No data points -> no extraction */
    if (inspec->n <= 0) {
        outspec->n = 0;
        return CPL_ERROR_NONE;
    }

    /* Test input parameters */
    if (limlam[1] <= limlam[0]) {
        /* Return spectrum with zero data points */
        outspec->n = 0;
        sprintf(errtxt, "%s: limlam[1] <= limlam[0] (wavelength limits)",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Get number of data points in given wavelength range */
    for (i = 0; i < inspec->n; i++) {
        if (inspec->dat[i].lam >= limlam[0]) {
            if (inspec->dat[i].lam > limlam[1]) {
                break;
            }
            if (count == 0) {
                pix1 = i;
            }
            count++;
        }
    }

    /* Allocate memory (set outspec->type = 3) */
    if ((int) sm_spec_malloc(outspec, count) == (int) SM_ERROR_ISM) {
        return SM_ERROR_ISM;
    }

    /* Copy type flag again */
    outspec->type = inspec->type;

    /* Copy wavelengths and fluxes */
    if (outspec->n > 0) {
        for (i = 0; i < outspec->n; i++) {
            outspec->dat[i].lam = inspec->dat[i+pix1].lam;
            outspec->dat[i].flux = inspec->dat[i+pix1].flux;
            if (outspec->type >= 3) {
                outspec->dat[i].dflux1 = inspec->dat[i+pix1].dflux1;
                outspec->dat[i].dflux2 = inspec->dat[i+pix1].dflux2;
            }
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_rebin(smspec *outspec, const smspec *inspec)
{
    /*!
     * Rebins ::smspec structure by using bins with constant fluxes.
     * Conserves integral of fluxes.
     * Errors are rebinned in the same way as fluxes
     * \f$\to\f$ not suitable for statistical noise.
     *
     * \b INPUT:
     * \param outspec  desired wavelength grid
     * \param inspec   wavelengths and fluxes (per wavelength unit) of input
     *                 spectrum
     *
     * \b OUTPUT:
     * \param outspec  rebinned input spectrum
     *
     * \b ERRORS:
     * - none
     */

    int i = 0, jo = -1, j = 0;
    double olmin = 0., olmax = 0., dol = 0., ilmin = 0., ilmax = 0., dil = 0.;
    double rdl = 0.;

    /* Copy type flag */
    outspec->type = inspec->type;

    /* No data points -> no rebinning */
    if (inspec->n <= 0 || outspec->n <= 0) {
        return CPL_ERROR_NONE;
    }

    /* One input data point only */
    if (inspec->n == 1) {
        for (j = 0; j < outspec->n; j++) {
            outspec->dat[j].flux = inspec->dat[0].flux;
            if (outspec->type >= 3) {
                outspec->dat[j].dflux1 = inspec->dat[0].dflux1;
                outspec->dat[j].dflux2 = inspec->dat[0].dflux2;
            }
        }
        return CPL_ERROR_NONE;
    }

    for (i = 0; i < outspec->n; i++) {

        /* Limits of wavelength bin in output spectrum */

        if (outspec->n == 1) {

          /* Full range of input spectrum for one output data point */
          olmin = 1.5 * inspec->dat[0].lam
                - 0.5 * inspec->dat[1].lam;
          olmax = 1.5 * inspec->dat[(inspec->n)-1].lam
                - 0.5 * inspec->dat[(inspec->n)-2].lam;

        } else {

            if (i == 0) {
                olmin = 1.5 * outspec->dat[i].lam
                      - 0.5 * outspec->dat[i+1].lam;
            } else {
                olmin = olmax;
            }

            if (i == outspec->n - 1) {
                olmax = 1.5 * outspec->dat[i].lam
                      - 0.5 * outspec->dat[i-1].lam;
            } else {
                olmax = 0.5 * (outspec->dat[i].lam + outspec->dat[i+1].lam);
            }

        }

        dol = olmax - olmin;

        outspec->dat[i].flux = 0.;
        if (outspec->type >= 3) {
            outspec->dat[i].dflux1 = 0.;
            outspec->dat[i].dflux2 = 0.;
        }

        do {

            /* Limits of wavelength bin in input spectrum */

            if (j != jo) {

                if (j == 0) {
                    ilmin = 1.5 * inspec->dat[j].lam
                          - 0.5 * inspec->dat[j+1].lam;
                } else {
                    ilmin = ilmax;
                }

                if (j == inspec->n - 1) {
                    ilmax = 1.5 * inspec->dat[j].lam
                          - 0.5 * inspec->dat[j-1].lam;
                } else {
                    ilmax = 0.5 * (inspec->dat[j].lam + inspec->dat[j+1].lam);
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
                outspec->dat[i].flux += inspec->dat[j].flux * rdl;
                if (outspec->type >= 3) {
                    /* Not suitable for statistical noise */
                    outspec->dat[i].dflux1 += inspec->dat[j].dflux1 * rdl;
                    outspec->dat[i].dflux2 += inspec->dat[j].dflux2 * rdl;
                }
            }

            j++;

        } while (ilmax <= olmax && j < inspec->n);

        j--;
        jo = j;

    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_interpol(smspec *outspec, const smspec *inspec)
{
    /*!
     * Linear interpolation of ::smspec structure for conversion to finer
     * wavelength grid.
     * Optimal for sparse, unbinned data.
     *
     * \b INPUT:
     * \param outspec  desired wavelength grid
     * \param inspec   wavelengths and fluxes (per wavelength unit) of input
     *                 spectrum
     *
     * \b OUTPUT:
     * \param outspec  interpolated input spectrum
     *
     * \b ERRORS:
     * - none
     */

    int i, j = 0;
    double dlam = 0., m = 0., me1 = 0., me2 = 0., dx = 0.;

    /* Copy type flag */
    outspec->type = inspec->type;

    /* No data points -> no interpolation */
    if (inspec->n <= 0 || outspec->n <= 0) {
        return CPL_ERROR_NONE;
    }

    /* One input data point only */
    if (inspec->n == 1) {
        for (j = 0; j < outspec->n; j++) {
            outspec->dat[j].flux = inspec->dat[0].flux;
            if (outspec->type >= 3) {
                outspec->dat[j].dflux1 = inspec->dat[0].dflux1;
                outspec->dat[j].dflux2 = inspec->dat[0].dflux2;
            }
        }
        return CPL_ERROR_NONE;
    }

    /* Constant extrapolation below lower limit */
    while (j < outspec->n && outspec->dat[j].lam < inspec->dat[0].lam) {
        outspec->dat[j].flux = inspec->dat[0].flux;
        if (outspec->type >= 3) {
            outspec->dat[j].dflux1 = inspec->dat[0].dflux1;
            outspec->dat[j].dflux2 = inspec->dat[0].dflux2;
        }
        j++;
    }

    /* Interpolation */

    for (i = 1; i < inspec->n; i++) {

        /* Slope for flux */
        dlam = inspec->dat[i].lam - inspec->dat[i-1].lam;
        m = (inspec->dat[i].flux - inspec->dat[i-1].flux) / dlam;

        /* Slope for error */
        if (outspec->type == 3) {
            me1 = (inspec->dat[i].dflux1 - inspec->dat[i-1].dflux1) / dlam;
        } else if (outspec->type == 4) {
            me2 = (inspec->dat[i].dflux2 - inspec->dat[i-1].dflux2) / dlam;
        }

        /* Linear interpolation */
        while (j < outspec->n && outspec->dat[j].lam < inspec->dat[i].lam) {
            dx = outspec->dat[j].lam - inspec->dat[i-1].lam;
            outspec->dat[j].flux = m * dx + inspec->dat[i-1].flux;
            if (outspec->type >= 3) {
                outspec->dat[j].dflux1 = me1 * dx + inspec->dat[i-1].dflux1;
            }
            if (outspec->type == 3) {
                outspec->dat[j].dflux2 = outspec->dat[j].dflux1;
            } else if (outspec->type == 4) {
                outspec->dat[j].dflux2 = me2 * dx + inspec->dat[i-1].dflux2;
            }
            j++;
        }

    }

    /* Constant extrapolation above upper limit */
    while (j < outspec->n &&
           outspec->dat[j].lam >= inspec->dat[(inspec->n)-1].lam) {
        outspec->dat[j].flux = inspec->dat[(inspec->n)-1].flux;
        if (outspec->type >= 3) {
            outspec->dat[j].dflux1 = inspec->dat[(inspec->n)-1].dflux1;
            outspec->dat[j].dflux2 = inspec->dat[(inspec->n)-1].dflux2;
        }
        j++;
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_spec_convolve(smspec *spec, const int nkpix,
                                const double *kernel)
{
    /*!
     * Convolution of ::smspec structure with given kernel.
     * Errors are convolved in the same way as fluxes
     * \f$\to\f$ not suitable for statistical noise.
     *
     * \note The centre of the convolution function is shifted by -0.5 pixels
     *       for an even number of kernel pixels.
     *
     * \b INPUT:
     * \param spec    input spectrum (::smspec structure)
     * \param nkpix   number of kernel pixels
     * \param kernel  vector of kernel values
     *                (required: sum of all values = 1)
     *
     * \b OUTPUT:
     * \param spec    convolved input spectrum
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     */

    smspec tempspec;
    char errtxt[SM_MAXLEN+1];
    int k, kmin, kmax, i, j;
    double sum = 0.;

    /* No data points -> no convolution */

    if (spec->n <= 0) {
        return CPL_ERROR_NONE;
    }

    /* No kernel data -> no convolution */

    if (nkpix <= 0) {
        return CPL_ERROR_NONE;
    }

    /* Check kernel */

    if (nkpix > 0 && kernel == NULL) {
        sprintf(errtxt, "%s: nkpix (N_kernelpixel > 0 but kernel == NULL)",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    for (k = 0; k < nkpix; k++) {
        if (kernel[k] < 0 || kernel[k] > 1) {
            sprintf(errtxt, "%s: *kernel (kernel element(s) < 0 or > 1)",
                    SM_ERROR_IIP_TXT);
            return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s",
                                         errtxt);
        }
        sum += kernel[k];
    }

    if (sum < 1 - SM_TOL || sum > 1 + SM_TOL) {
        sprintf(errtxt, "%s: *kernel (sum of kernel elements != 1)",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Skip convolution if number of kernel pixels is one */

    if (nkpix == 1) {
        return CPL_ERROR_NONE;
    }

    /* Create temporary spectrum */

    sm_spec_copy(&tempspec, spec);
    sm_spec_scale(spec, '=', 0.);

    /* Kernel with even or odd pixel number?
       Note: centre of kernel at -0.5 pixels for even pixel number */

    if (nkpix % 2 == 0) {
        kmin = - nkpix / 2;
    } else {
        kmin = - (nkpix - 1) / 2;
    }
    kmax = kmin + nkpix - 1;

    /* Convolve spectrum with kernel */

    for (i = 0; i < spec->n; i++) {

        for (k = kmin; k <= kmax; k++) {

            j = i + k;

            /* Flux of first or last valid pixel for invalid pixels */
            if (j < 0) {
                j = 0;
            } else if (j >= spec->n) {
                j = spec->n - 1;
            }

            spec->dat[i].flux += tempspec.dat[j].flux * kernel[k - kmin];
            if (spec->type >= 3) {
                /* Not suitable for statistical noise */
                spec->dat[i].dflux1 +=
                    tempspec.dat[j].dflux1 * kernel[k - kmin];
            }
            if (spec->type == 4) {
                spec->dat[i].dflux2 +=
                    tempspec.dat[j].dflux2 * kernel[k - kmin];
            }

        }

        if (spec->type == 3) {
            spec->dat[i].dflux2 = spec->dat[i].dflux1;
        }

    }

    /* Delete temporary spectrum */

    sm_spec_free(&tempspec);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_param_read(FILE *stream, smparam par[])
{
    /*!
     * Reads a line from a file and returns the space-, tab-, or '='-separated
     * individual strings in an array of ::SM_MAXPAR ::smparam structures. A
     * ::smparam structure contains each string converted into integer and
     * double precision floating point numbers. Comment lines marked by # and
     * empty lines are not considered. In this case, the next line with data
     * is read.
     *
     * \b INPUT:
     * \param stream  name of stream (to be defined by fopen command)
     * \param par     array of ::SM_MAXPAR ::smparam structures
     *                (to transfer parameter values)
     *
     * \b OUTPUT:
     * \param par     ::smparam array that contains the read value(s) as
     *                character string ("c"), integer ("i"), or double
     *                precision floating point number ("d").
     *                The different data types can be accessed by adding the
     *                corresponding suffixes "c", "i", or "d" to the name of
     *                the ::smparam structure.
     *                An additional suffix "n" element indicates the number
     *                of read values minus the array index, i.e. \e par[0].n
     *                provides the number of values in the read line.
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     * - UFS: Unexpected file structure
     * - ISM: Insufficient memory
     *
     * In the case of problems, the default values "", 0, and 0.0 are used
     * for the corresponding data types.
     */

    char errtxt[SM_MAXLEN+1], line[SM_MAXLEN+2], **str = NULL;
    cpl_boolean fl_data = CPL_FALSE;
    int npar0 = 0, i = 0, nparx = 0;

    /* Set maximum number of parameters */
    npar0 = SM_MAXPAR;

    /* Valid number of parameters? */

    if (npar0 < 1) {
        sprintf(errtxt, "%s: npar < 1 (parameters at line)",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Fill smparam structure with 0 */

    for (i = 0; i < npar0; i++) {
        par[i].c[0] = '\0';
        par[i].i = 0;
        par[i].d = 0.;
        par[i].n = 0;
    }

    /* Skip comments and empty lines and check existence of data */

    while (fgets(line, SM_MAXLEN+2, stream) != NULL) {
        if (line[0] != '#' && strspn(line, "\n\t =") != strlen(line)) {
            fl_data = CPL_TRUE;
            break;
        }
    }

    if (fl_data == CPL_FALSE) {
        sprintf(errtxt, "%s: FILE *stream (no parameters)",
                SM_ERROR_UFS_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s", errtxt);
    }

    /* Allocate memory to string array */

    str = (char **) calloc(npar0+1, sizeof(char *));
    if (str == NULL) {
        sprintf(errtxt, "%s: char **str", SM_ERROR_ISM_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_ISM, "%s", errtxt);
    }

    /* Read parameter values in strings and check number of values */

    str[0] = strtok(line, "\n\t =");
    if (npar0 > 1) {
        for (i = 1; i < npar0; i++) {
            str[i] = strtok(NULL, "\n\t =");
            if (str[i] == NULL) {
                nparx = i;
                break;
            }
        }
    }

    if (nparx == npar0) {
        str[npar0] = strtok(NULL, "\n\t =");
        if (str[npar0] != NULL) {
            sprintf(errtxt, "%s: FILE *stream (more parameters than "
                    "expected)", SM_ERROR_UFS_TXT);
            return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s",
                                         errtxt);
        }
    }

    /* Convert strings in integer and double numbers */

    for (i = 0; i < nparx; i++) {
        strcpy(par[i].c, str[i]);
        strtok(par[i].c, "\r");
        par[i].i = strtol(str[i], NULL, 10);
        if (sm_basic_isnumber(str[i]) == CPL_TRUE) {
            par[i].d = strtod(str[i], NULL);
        } else {
            par[i].d = 0.;
        }
        par[i].n = nparx - i;
    }

    /* Free memory space of string array */
    free(str);
    str = NULL;

    return CPL_ERROR_NONE;
}


cpl_error_code sm_param_readcheck(FILE *stream, smparam par[],
                                  const char *parname, const int npar)
{
    /*!
     * Reads and checks a line from a file and returns the space-, tab-, or
     * '='-separated individual strings in an array of ::SM_MAXPAR ::smparam
     * structures. A ::smparam structure contains each string converted into
     * integer and double precision floating point numbers. Comment lines
     * marked by # and empty lines are not considered. In this case, the next
     * line with data is read.
     *
     * The first string of the read line is expected to be the parameter name.
     * It is checked whether it matches the input name. In the case of a
     * deviation, an error is returned. The routine also checks the number of
     * parameter values (not counting the parameter name), which is also an
     * input parameter of the routine. An error is returned if there is a
     * mismatch.
     *
     * \b INPUT:
     * \param stream   name of stream (to be defined by fopen command)
     * \param par      array of ::SM_MAXPAR ::smparam structures
     *                 (to transfer parameter values)
     * \param parname  name of parameter
     * \param npar     number of expected parameter values
     *
     * \b OUTPUT:
     * \param par      ::smparam array that contains the read value(s) as
     *                 character string ("c"), integer ("i"), or double
     *                 precision floating point number ("d").
     *                 The different data types can be accessed by adding the
     *                 corresponding suffixes "c", "i", or "d" to the name of
     *                 the ::smparam structure.
     *                 An additional suffix "n" element indicates the number
     *                 of read values minus the array index, i.e. \e par[0].n
     *                 provides the number of values in the read line.
     *
     * \b ERRORS:
     * - EIS: Error in subroutine (see ::sm_param_read)
     * - IOV: Invalid object value(s)
     *
     * In the case of problems, the default values "", 0, and 0.0 are used
     * for the corresponding data types.
     */

    cpl_error_code status = CPL_ERROR_NONE;
    char errtxt[SM_MAXLEN+1];
    int i = 0;

    /* Read data from line */
    status = sm_param_read(stream, par);

    /* Check for errors returned by sm_param_read */
    if (status != CPL_ERROR_NONE) {
        sprintf(errtxt, "%s: FILE *stream (reading of parameter %s failed)",
                SM_ERROR_EIS_TXT, parname);
        return cpl_error_set_message(cpl_func, SM_ERROR_EIS, "%s", errtxt);
    }

    /* Check the number of parameter values */
    if (par[0].n != npar + 1) {
        sprintf(errtxt, "%s: smparam par[] (unexpected number of values for "
                "parameter %s)", SM_ERROR_IOV_TXT, parname);
        for (i = 0; i < par[0].n; i++) {
            par[i].c[0] = '\0';
            par[i].i = 0;
            par[i].d = 0.;
            par[i].n = 0;
        }
        return cpl_error_set_message(cpl_func, SM_ERROR_IOV, "%s", errtxt);
    }

    /* Check the parameter name */
    if (strcmp(par[0].c, parname) != 0) {
        sprintf(errtxt, "%s: smparam par[] (parameter %s instead of %s)",
                SM_ERROR_IOV_TXT, par[0].c, parname);
        for (i = 0; i < par[0].n; i++) {
            par[i].c[0] = '\0';
            par[i].i = 0;
            par[i].d = 0.;
            par[i].n = 0;
        }
        return cpl_error_set_message(cpl_func, SM_ERROR_IOV, "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_param_check(smparam par[], const char *parname,
                              const int npar)
{
    /*!
     * Checks a parameter definition held by an array of ::smparam structures.
     * The first string of the read line is expected to be the parameter name.
     * It is checked whether it matches the input name. In the case of a
     * deviation, an error is returned. The routine also checks the number of
     * parameter values (not counting the parameter name), which is also an
     * input parameter of the routine. An error is returned if there is a
     * mismatch. In the case of errors, all entries of the ::smparam array
     * are set to the default values "", 0, and 0.0.
     *
     * \b INPUT:
     * \param par      array of ::SM_MAXPAR ::smparam structures
     * \param parname  name of parameter
     * \param npar     number of expected parameter values
     *
     * \b OUTPUT:
     * \param par      ::smparam array (error: set to the data type default
     *                 values, otherwise untouched)
     *
     * \b ERRORS:
     * - EIS: Error in subroutine (see ::sm_param_read)
     * - IOV: Invalid object value(s)
     */

    char errtxt[SM_MAXLEN+1];
    int i = 0;

    /* Check for errors returned by sm_param_read */
    if (par[0].c[0] == '\0') {
        sprintf(errtxt, "%s: FILE *stream (reading of parameter %s failed)",
                SM_ERROR_EIS_TXT, parname);
        return cpl_error_set_message(cpl_func, SM_ERROR_EIS, "%s", errtxt);
    }

    /* Check the parameter name */
    if (strcmp(par[0].c, parname) != 0) {
        sprintf(errtxt, "%s: smparam par[] (parameter %s instead of %s)",
                SM_ERROR_IOV_TXT, par[0].c, parname);
        for (i = 0; i < npar+1; i++) {
            par[i].c[0] = '\0';
            par[i].i = 0;
            par[i].d = 0.;
            par[i].n = 0;
        }
        return cpl_error_set_message(cpl_func, SM_ERROR_IOV, "%s", errtxt);
    }

    /* Check the number of parameter values */
    if (par[0].n != npar + 1) {
        sprintf(errtxt, "%s: smparam par[] (unexpected number of values for "
                "parameter %s)", SM_ERROR_IOV_TXT, parname);
        for (i = 0; i < npar+1; i++) {
            par[i].c[0] = '\0';
            par[i].i = 0;
            par[i].d = 0.;
            par[i].n = 0;
        }
        return cpl_error_set_message(cpl_func, SM_ERROR_IOV, "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_grid_malloc(smgrid *xy, const int nx, const int ny)
{
    /*!
     * Allocates memory for an ::smgrid structure of size
     * \f$n_x \times n_y\f$. All values are zero.
     *
     * \b INPUT:
     * \param xy  ::smgrid structure of size zero
     * \param nx  number of \e x coordinates
     * \param ny  number of \e y coordinates
     *
     * \b OUTPUT:
     * \param xy  ::smgrid structure of size \f$n_x \times n_y\f$ if there is
     *            sufficient memory (otherwise size remains zero)
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     * - ISM: Insufficient memory
     */

    cpl_error_code err = CPL_ERROR_NONE;
    smbool fl_size = T, fl_mem = T;
    char errtxt[SM_MAXLEN+1];
    int i;

    /* Get number of data points */

    if (nx <= 0 || ny <= 0) {
        /* Return spectrum with zero data points */
        xy->nx = 0;
        xy->ny = 0;
        fl_size = F;
    } else {
        xy->nx = nx;
        xy->ny = ny;
    }

    /* Allocate memory */

    if (fl_size == T && fl_mem == T) {
        for (i = 0; i < 2; i++) {
            xy->xpos = (double *) calloc(xy->nx, sizeof(double));
            xy->ypos = (double *) calloc(xy->ny, sizeof(double));
            xy->val = (double **) calloc(xy->nx, sizeof(double));
            if (xy->xpos == NULL || xy->ypos == NULL || xy->val == NULL) {
                xy->nx = 0;
                xy->ny = 0;
                fl_mem = F;
            } else {
                break;
            }
        }
    }

    if (fl_size == T && fl_mem == T) {
        for (i = 0; i < xy->nx; i++) {
            xy->val[i] = (double *) calloc(xy->ny, sizeof(double));
            if (xy->val[i] == NULL && xy->ny > 0) {
                xy->ny = 0;
                fl_mem = F;
                i = 0;
            }
        }
        if (xy->ny == 0) {
            xy->nx = 0;
        }
    }

    /* Manage errors */

    if (fl_mem == F) {
        sprintf(errtxt, "%s: smgrid *xy", SM_ERROR_ISM_TXT);
        err = cpl_error_set_message(cpl_func, SM_ERROR_ISM, "%s", errtxt);
    }

    if (fl_size == F) {
        sprintf(errtxt, "%s: nx <= 0 || ny <= 0", SM_ERROR_IIP_TXT);
        err = cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    return err;
}


cpl_error_code sm_grid_read(smgrid *xy, const char *filename)
{
    /*!
     * Reads ::smgrid structure from ASCII file.
     *
     * Required file structure:
     * - Comment lines have to be marked by #.
     * - Data lines:
     *   - 1: \e nx \e ny
     *   - 2: \e xpos (\e nx values)
     *   - 3: \e ypos (\e ny values)
     *   - 4-(\e nx+3): \e val (\e ny values each line)
     *
     * \b INPUT:
     * \param filename  name of input file
     *
     * \b OUTPUT:
     * \param xy        output ::smgrid structure
     *
     * \b ERRORS:
     * - FOF: File opening failed
     * - UFS: Unexpected file structure
     * - ISM: Insufficient memory
     */

    FILE *stream;
    char errtxt[SM_MAXLEN+1], line[SM_MAXLEN+2], **str = NULL;
    int datline = 0, nval = 0, nx = 0, ny = 0, i = 0;

    xy->nx = 0;
    xy->ny = 0;

    /* Check file existence */

    if ((stream = fopen(filename, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", SM_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, SM_ERROR_FOF, "%s", errtxt);
    }

    /* Read file line by line */

    while (fgets(line, SM_MAXLEN+2, stream) != NULL) {

        if (line[0] == '#') {
            continue;
        }

        if (datline == 3 + nx) {
            fclose(stream);
            sprintf(errtxt, "%s: %s (more lines than expected)",
                    SM_ERROR_UFS_TXT, filename);
            return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s",
                                         errtxt);
        }

        /* Number of expected values */

        datline++;
        if (datline == 1) {
            nval = 2;
        } else if (datline == 2) {
            nval = nx;
        } else if (datline == 3) {
            nval = ny;
        }

        /* Allocate memory to string array */

        if (datline < 4) {
            str = (char **) calloc(nval, sizeof(char *));
            if (str == NULL) {
                fclose(stream);
                sprintf(errtxt, "%s: char **str", SM_ERROR_ISM_TXT);
                return cpl_error_set_message(cpl_func, SM_ERROR_ISM, "%s",
                                             errtxt);
            }
        }

        /* Read values in strings */

        str[0] = strtok(line, "\n\t ");
        if (nval > 1) {
            for (i = 1; i < nval; i++) {
                str[i] = strtok(NULL, "\n\t ");
                if (str[i] == NULL) {
                    fclose(stream);
                    sprintf(errtxt,
                            "%s: %s (less parameters at line than expected)",
                            SM_ERROR_UFS_TXT, filename);
                    return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s",
                                                 errtxt);
                }
            }
        }

        /* Convert strings in integer or double numbers */

        if (datline == 1) {
            nx = strtol(str[0], NULL, 10);
            ny = strtol(str[1], NULL, 10);
            /* Create smgrid structure */
            if ((int) sm_grid_malloc(xy, nx, ny) == (int) SM_ERROR_ISM) {
                fclose(stream);
                return SM_ERROR_ISM;
            }
        } else if (datline == 2) {
            for (i = 0; i < nval; i++) {
                xy->xpos[i] = strtod(str[i], NULL);
            }
        } else if (datline == 3) {
            for (i = 0; i < nval; i++) {
                xy->ypos[i] = strtod(str[i], NULL);
            }
        } else {
            for (i = 0; i < nval; i++) {
                xy->val[datline-4][i] = strtod(str[i], NULL);
            }
        }

        /* Free memory space of string array */

        if (datline < 3 || datline == 3 + nx) {
            free(str);
            str = NULL;
        }

    }

    if (datline < 3 + nx) {
        sprintf(errtxt, "%s: %s (less lines than expected)",
                SM_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s", errtxt);
    }

    fclose(stream);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_grid_extract(double *outval, const smgrid *xy,
                               const double x0, const double y0)
{
    /*!
     * Extracts value of an ::smgrid structure at position
     * (\f$x_0\f$, \f$y_0\f$) by means of bilinear interpolation
     *
     * \b INPUT:
     * \param xy      input ::smgrid structure
     * \param x0      \e x coordinate position for extraction
     * \param y0      \e y coordinate position for extraction
     *
     * \b OUTPUT:
     * \param outval  value of ::smgrid structure at position
     *                (\f$x_0\f$, \f$y_0\f$)
     *
     * \b ERRORS:
     * - NDA: No data
     * - ISD: Insufficient data points
     * - IOD: Invalid order of data points
     * - IIP: Invalid input parameter(s)
     */

    char errtxt[SM_MAXLEN+1];
    int i = 0, ip0 = 0, ip1 = 0, j = 0, jp0 = 0, jp1 = 0;
    double xsig = 0, ysig = 0, dx = 0, dy = 0, fx = 0, fy = 0;

    /* Default output value */

    *outval = 0;

    /* Basic tests of input data */

    if (xy->nx < 1 || xy->ny < 1) {
        sprintf(errtxt, "%s: smgrid *xy", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    }
    if (xy->nx == 1 || xy->ny == 1) {
        sprintf(errtxt, "%s: smgrid *xy (xy->nx == 1 || xy->ny == 1)",
                SM_ERROR_ISD_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_ISD, "%s", errtxt);
    }

    /* Ascending or descending coordinates? -> Derive signum */

    xsig = xy->xpos[(xy->nx)-1] - xy->xpos[0];
    xsig /= fabs(xsig);
    ysig = xy->ypos[(xy->ny)-1] - xy->ypos[0];
    ysig /= fabs(ysig);

    /* Monotonic increase/decrease of coordinates? */

    for (i = 1; i < xy->nx; i++) {
        if ((xsig > 0 && xy->xpos[i-1] > xy->xpos[i]) ||
            (xsig < 0 && xy->xpos[i-1] < xy->xpos[i])) {
            sprintf(errtxt, "%s: smgrid *xy (x coordinate xy->xpos)",
                    SM_ERROR_IOD_TXT);
            return cpl_error_set_message(cpl_func, SM_ERROR_IOD, "%s",
                                         errtxt);
        }
    }

    for (j = 1; j < xy->ny; j++) {
        if ((ysig > 0 && xy->ypos[j-1] > xy->ypos[j]) ||
            (ysig < 0 && xy->ypos[j-1] < xy->ypos[j])) {
            sprintf(errtxt, "%s: smgrid *xy (y coordinate xy->ypos)",
                    SM_ERROR_IOD_TXT);
            return cpl_error_set_message(cpl_func, SM_ERROR_IOD, "%s",
                                         errtxt);
        }
    }

    /* Desired coordinates inside the grid? */

    if ((xsig > 0 && (x0 < xy->xpos[0] || x0 > xy->xpos[(xy->nx)-1])) ||
        (xsig < 0 && (x0 > xy->xpos[0] || x0 < xy->xpos[(xy->nx)-1])) ||
        (ysig > 0 && (y0 < xy->ypos[0] || y0 > xy->ypos[(xy->ny)-1])) ||
        (ysig < 0 && (y0 > xy->ypos[0] || y0 < xy->ypos[(xy->ny)-1]))) {
        sprintf(errtxt,
                "%s: x0, y0 (selected coordinates outside of grid smgrid *xy",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Determine the grid positions that encircle the desired coordinates */

    for (i = 1; i < xy->nx; i++) {
        dx = x0 - xy->xpos[i];
        if (xsig > 0) {
          if (dx <= 0) {
                ip0 = i - 1;
                ip1 = i;
                break;
            }
        } else {
            if (dx >= 0) {
                ip0 = i;
                ip1 = i - 1;
                break;
            }
        }
    }

    for (j = 1; j < xy->ny; j++) {
        dy = y0 - xy->ypos[j];
        if (ysig > 0) {
            if (dy <= 0) {
                jp0 = j - 1;
                jp1 = j;
                break;
            }
        } else {
            if (dy >= 0) {
                jp0 = j;
                jp1 = j - 1;
                break;
            }
        }
    }

    /* Bilinear interpolation (approximation) */

    fx = (x0 - xy->xpos[ip0]) / (xy->xpos[ip1] - xy->xpos[ip0]);
    fy = (y0 - xy->ypos[jp0]) / (xy->ypos[jp1] - xy->ypos[jp0]);

    *outval = xy->val[ip0][jp0] * (1 - fx) * (1 - fy)
            + xy->val[ip1][jp0] * fx * (1 - fy)
            + xy->val[ip0][jp1] * (1 - fx) * fy
            + xy->val[ip1][jp1] * fx * fy;

    return CPL_ERROR_NONE;
}


cpl_error_code sm_grid_write(const smgrid *xy, const char *filename)
{
    /*!
     * Writes ::smgrid structure to ASCII file or on stdout
     * (indicated by "stdout" as filename)
     *
     * \b INPUT:
     * \param xy        ::smgrid structure
     * \param filename  name of output file or "stdout"

     * \b ERRORS:
     * - FOF: File opening failed
     */

    FILE *stream;
    char errtxt[SM_MAXLEN+1];
    int i, j;
    double min = 0.001, max = 100000; /* float format within this range */

    if (strncmp(filename, "stdout", 6) == 0) {
        /* Output on stdout */
        stream = stdout;
    } else {
        /* Check file existence */
        if ((stream = fopen(filename, "w")) == NULL) {
            sprintf(errtxt, "%s: %s", SM_ERROR_FOF_TXT, filename);
            return cpl_error_set_message(cpl_func, SM_ERROR_FOF, "%s",
                                         errtxt);
        }
    }

    /* Write file */

    fprintf(stream, "# nx ny:\n");
    fprintf(stream, "%d %d\n", xy->nx, xy->ny);

    if (xy->nx > 0 && xy->ny > 0) {

        fprintf(stream, "# xpos[x]:\n");
        for (i = 0; i < xy->nx; i++) {
            if ((fabs(xy->xpos[i]) >= min && fabs(xy->xpos[i]) < max)
                || xy->xpos[i] == 0.) {
                fprintf(stream, "%g ", xy->xpos[i]);
            } else {
                fprintf(stream, "%e ", xy->xpos[i]);
            }
        }
        fprintf(stream, "\n");

        fprintf(stream, "# ypos[y]:\n");
        for (j = 0; j < xy->ny; j++) {
            if ((fabs(xy->ypos[j]) >= min && fabs(xy->ypos[j]) < max)
                || xy->ypos[j] == 0.) {
                fprintf(stream, "%g ", xy->ypos[j]);
            } else {
                fprintf(stream, "%e ", xy->ypos[j]);
            }
        }
        fprintf(stream, "\n");

        fprintf(stream, "# val[x][y]:\n");
        for (i = 0; i < xy->nx; i++) {
            for (j = 0; j < xy->ny; j++) {
                if ((fabs(xy->val[i][j]) >= min && fabs(xy->val[i][j]) < max)
                    || xy->val[i][j] == 0.) {
                    fprintf(stream, "%g ", xy->val[i][j]);
                } else {
                    fprintf(stream, "%e ", xy->val[i][j]);
                }
            }
            fprintf(stream, "\n");
        }

    }

    if (strncmp(filename, "stdout", 6) != 0) {
        fclose(stream);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_grid_print(const smgrid *xy)
{
    /*!
     * Prints ::smgrid structure on stdout
     *
     * \b INPUT:
     * \param xy  ::smgrid structure
     *
     * \b ERRORS:
     * - none
     */

    sm_grid_write(xy, "stdout");

    return CPL_ERROR_NONE;
}


cpl_error_code sm_grid_free(smgrid *xy)
{
    /*!
     * Frees memory occupied by an ::smgrid structure
     *
     * \b INPUT:
     * \param xy  ::smgrid structure of size \f$n_x \times n_y\f$
     *
     * \b OUTPUT:
     * \param xy  ::smgrid structure of size zero
     *
     * \b ERRORS:
     * - none
     */

    int i;

    if (xy->nx > 0 &&  xy->ny > 0) {

        free(xy->xpos);
        xy->xpos = NULL;
        free(xy->ypos);
        xy->ypos = NULL;

        for (i = 0; i < xy->nx; i++) {
            free(xy->val[i]);
            xy->val[i] = NULL;
        }
        free(xy->val);
        xy->val = NULL;

    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_basic_chdir(const char *dir)
{
    /*!
     * \brief
     *   Change working directory
     *
     * This function provides a wrapper for the intrinsic function chdir().
     * It reports back all chdir() errors.
     *
     * \b INPUT:
     * \param dir   new working directory
     *
     * \b OUTPUT:
     * \return      CPL_ERROR_NONE on success, or on failure:
     *              SM_ERROR_ACCES,
     *              SM_ERROR_LOOP,
     *              SM_ERROR_NAMETOOLONG,
     *              SM_ERROR_NOENT,
     *              SM_ERROR_NOTDIR,
     *              SM_ERROR_FAULT,
     *              SM_ERROR_IO,
     *              SM_ERROR_NOMEM,
     *              SM_ERROR_UNDEF
     *
     * \b ERRORS:
     * - see chdir manual for details
     */

    int err;

    errno = 0;

    err = chdir(dir);

    switch (errno) {
    case EACCES:
        return cpl_error_set_message(cpl_func, SM_ERROR_ACCES,
                                     "%s: %s", SM_ERROR_ACCES_TXT, dir);
        break;
    case ELOOP:
        return cpl_error_set_message(cpl_func, SM_ERROR_LOOP,
                                     "%s: %s", SM_ERROR_LOOP_TXT, dir);
        break;
    case ENAMETOOLONG:
        return cpl_error_set_message(cpl_func, SM_ERROR_NAMETOOLONG,
                                     "%s: %s", SM_ERROR_NAMETOOLONG_TXT, dir);
        break;
    case ENOENT:
        return cpl_error_set_message(cpl_func, SM_ERROR_NOENT,
                                     "%s: %s", SM_ERROR_NOENT_TXT, dir);
        break;
    case ENOTDIR:
        return cpl_error_set_message(cpl_func, SM_ERROR_NOTDIR,
                                     "%s: %s", SM_ERROR_NOTDIR_TXT, dir);
        break;
    case EFAULT:
        return cpl_error_set_message(cpl_func, SM_ERROR_FAULT,
                                     "%s: %s", SM_ERROR_FAULT_TXT, dir);
        break;
    case EIO:
        return cpl_error_set_message(cpl_func, SM_ERROR_IO,
                                     "%s: %s", SM_ERROR_IO_TXT, dir);
        break;
    case ENOMEM:
        return cpl_error_set_message(cpl_func, SM_ERROR_NOMEM,
                                     "%s: %s", SM_ERROR_NOMEM_TXT, dir);
        break;
    }

    if (err != 0) {
        return cpl_error_set_message(cpl_func, SM_ERROR_UNDEF,
                                     "%s: %s", SM_ERROR_UNDEF_TXT, dir);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_basic_access(const char *pathname, const int mode)
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
     *                  SM_ERROR_ROFS,
     *                  SM_ERROR_INVAL,
     *                  SM_ERROR_TXTBSY,
     *                  SM_ERROR_ACCES,
     *                  SM_ERROR_LOOP,
     *                  SM_ERROR_NAMETOOLONG,
     *                  SM_ERROR_NOENT,
     *                  SM_ERROR_NOTDIR,
     *                  SM_ERROR_FAULT,
     *                  SM_ERROR_IO,
     *                  SM_ERROR_NOMEM,
     *                  SM_ERROR_UNDEF
     *
     * \b ERRORS:
     * - see access manual for details
     */

    int err;

    errno = 0;

    err = access(pathname, mode);

    switch (errno) {
    case EROFS:
        return cpl_error_set_message(cpl_func, SM_ERROR_ROFS,
                                     "%s: %s", SM_ERROR_ROFS_TXT, pathname);
        break;
    case EINVAL:
        return cpl_error_set_message(cpl_func, SM_ERROR_INVAL,
                                     "%s: %s", SM_ERROR_INVAL_TXT, pathname);
        break;
    case ETXTBSY:
        return cpl_error_set_message(cpl_func, SM_ERROR_TXTBSY,
                                     "%s: %s", SM_ERROR_TXTBSY_TXT, pathname);
        break;
    case EACCES:
        return cpl_error_set_message(cpl_func, SM_ERROR_ACCES,
                                     "%s: %s", SM_ERROR_ACCES_TXT, pathname);
        break;
    case ELOOP:
        return cpl_error_set_message(cpl_func, SM_ERROR_LOOP,
                                     "%s: %s", SM_ERROR_LOOP_TXT, pathname);
        break;
    case ENAMETOOLONG:
        return cpl_error_set_message(cpl_func, SM_ERROR_NAMETOOLONG,
                                     "%s: %s", SM_ERROR_NAMETOOLONG_TXT,
                                     pathname);
        break;
    case ENOENT:
        return cpl_error_set_message(cpl_func, SM_ERROR_NOENT,
                                     "%s: %s", SM_ERROR_NOENT_TXT, pathname);
        break;
    case ENOTDIR:
        return cpl_error_set_message(cpl_func, SM_ERROR_NOTDIR,
                                     "%s: %s", SM_ERROR_NOTDIR_TXT, pathname);
        break;
    case EFAULT:
        return cpl_error_set_message(cpl_func, SM_ERROR_FAULT,
                                     "%s: %s", SM_ERROR_FAULT_TXT, pathname);
        break;
    case EIO:
        return cpl_error_set_message(cpl_func, SM_ERROR_IO,
                                     "%s: %s", SM_ERROR_IO_TXT, pathname);
        break;
    case ENOMEM:
        return cpl_error_set_message(cpl_func, SM_ERROR_NOMEM,
                                     "%s: %s", SM_ERROR_NOMEM_TXT, pathname);
        break;
    }

    if (err != 0) {
        return cpl_error_set_message(cpl_func, SM_ERROR_UNDEF,
                                     "%s: %s", SM_ERROR_UNDEF_TXT, pathname);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_basic_mkdir(const char *dir, const mode_t mode)
{
    /*!
     * \brief
     *   Create directory
     *
     * This function provides a wrapper for the intrinsic function mkdir().
     * It reports back all mkdir() errors.
     *
     * \b INPUT:
     * \param dir   name of directory
     * \param mode  access permissions
     *
     * \b OUTPUT:
     * \return      CPL_ERROR_NONE on success, or on failure:
     *              SM_ERROR_EXIST,
     *              SM_ERROR_NOSPC,
     *              SM_ERROR_PERM,
     *              SM_ERROR_ROFS,
     *              SM_ERROR_INVAL,
     *              SM_ERROR_TXTBSY,
     *              SM_ERROR_ACCES,
     *              SM_ERROR_LOOP,
     *              SM_ERROR_NAMETOOLONG,
     *              SM_ERROR_NOENT,
     *              SM_ERROR_NOTDIR,
     *              SM_ERROR_FAULT,
     *              SM_ERROR_IO,
     *              SM_ERROR_NOMEM,
     *              SM_ERROR_UNDEF
     *
     * \b ERRORS:
     * - see mkdir manual for details
     */

    int err;

    errno = 0;

    err = mkdir(dir, mode);

    switch (errno) {
    case EEXIST:
        return cpl_error_set_message(cpl_func, SM_ERROR_EXIST,
                                     "%s: %s", SM_ERROR_EXIST_TXT, dir);
        break;
    case ENOSPC:
        return cpl_error_set_message(cpl_func, SM_ERROR_NOSPC,
                                     "%s: %s", SM_ERROR_NOSPC_TXT, dir);
        break;
    case EPERM:
        return cpl_error_set_message(cpl_func, SM_ERROR_PERM,
                                     "%s: %s", SM_ERROR_PERM_TXT, dir);
        break;
    case EROFS:
        return cpl_error_set_message(cpl_func, SM_ERROR_ROFS,
                                     "%s: %s", SM_ERROR_ROFS_TXT, dir);
        break;
    case EINVAL:
        return cpl_error_set_message(cpl_func, SM_ERROR_INVAL,
                                     "%s: %s", SM_ERROR_INVAL_TXT, dir);
        break;
    case ETXTBSY:
        return cpl_error_set_message(cpl_func, SM_ERROR_TXTBSY,
                                     "%s: %s", SM_ERROR_TXTBSY_TXT, dir);
        break;
    case EACCES:
        return cpl_error_set_message(cpl_func, SM_ERROR_ACCES,
                                     "%s: %s", SM_ERROR_ACCES_TXT, dir);
        break;
    case ELOOP:
        return cpl_error_set_message(cpl_func, SM_ERROR_LOOP,
                                     "%s: %s", SM_ERROR_LOOP_TXT, dir);
        break;
    case ENAMETOOLONG:
        return cpl_error_set_message(cpl_func, SM_ERROR_NAMETOOLONG,
                                     "%s: %s", SM_ERROR_NAMETOOLONG_TXT, dir);
        break;
    case ENOENT:
        return cpl_error_set_message(cpl_func, SM_ERROR_NOENT,
                                     "%s: %s", SM_ERROR_NOENT_TXT, dir);
        break;
    case ENOTDIR:
        return cpl_error_set_message(cpl_func, SM_ERROR_NOTDIR,
                                     "%s: %s", SM_ERROR_NOTDIR_TXT, dir);
        break;
    case EFAULT:
        return cpl_error_set_message(cpl_func, SM_ERROR_FAULT,
                                     "%s: %s", SM_ERROR_FAULT_TXT, dir);
        break;
    case EIO:
        return cpl_error_set_message(cpl_func, SM_ERROR_IO,
                                     "%s: %s", SM_ERROR_IO_TXT, dir);
        break;
    case ENOMEM:
        return cpl_error_set_message(cpl_func, SM_ERROR_NOMEM,
                                     "%s: %s", SM_ERROR_NOMEM_TXT, dir);
        break;
    }

    if (err != 0) {
        return cpl_error_set_message(cpl_func, SM_ERROR_UNDEF,
                                     "%s: %s", SM_ERROR_UNDEF_TXT, dir);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_basic_createdir(const char *dir, const mode_t mode)
{
    /*!
     * \brief
     *   Create directory including error handling
     *
     * This function creates a new directory with the specified permissions.
     *
     * \b INPUT:
     * \param dir   name of directory
     * \param mode  access permissions
     *
     * \b OUTPUT:
     * \return      CPL_ERROR_NONE on success, SM_ERROR_SUBROUTINE else.
     *
     * \b ERRORS:
     * - see mkdir manual for details
     *
     * \b NOTE:
     * - routine uses ::sm_basic_mkdir
     */

    char str[FILENAME_MAX];

    cpl_error_code err_code;
    cpl_errorstate err_state;

    err_state = cpl_errorstate_get();

    sprintf(str, "%s", dir);
    err_code = sm_basic_access(str, R_OK | W_OK);
    if ((int) err_code == (int) SM_ERROR_NOENT) {
        // directory does not exist
        cpl_errorstate_set(err_state); // clean error state
        if (sm_basic_mkdir(str, mode) != CPL_ERROR_NONE) {
            // directory cannot be created
            return cpl_error_set_message(cpl_func, SM_ERROR_SUBROUTINE, "%s",
                                         SM_ERROR_SUBROUTINE_TXT);
        }
    } else if (err_code != 0) {
        // some other error
        return cpl_error_set_message(cpl_func, SM_ERROR_SUBROUTINE, "%s",
                                     SM_ERROR_SUBROUTINE_TXT);
    }

    return CPL_ERROR_NONE;
}


void sm_basic_initstring(char *str, const long n)
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


cpl_boolean sm_basic_isnumber(char *str)
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


char *sm_basic_replacestring(char *instring, char *oldsubstr, char *newsubstr)
{
    /*!
     * \brief Substitutes substring in string
     *
     * \param instring   input string
     * \param oldsubstr  substring to be replaced
     * \param newsubstr  new substring
     *
     * \return input string with substituted substring
     */

  static char buffer[FILENAME_MAX];
  char *ptr;

  if(!(ptr = strstr(instring, oldsubstr)))
  {
      return instring;
  }
  strncpy(buffer, instring, ptr-instring);
  buffer[ptr-instring] = '\0';  /* terminate */
  sprintf(buffer+(ptr-instring), "%s%s", newsubstr, ptr+strlen(oldsubstr));

  return buffer;
}


char *sm_basic_rmcntrl(char *str)
{
    /*!
     * \brief
     *   Remove control characters from string.
     *
     * This function removes all \c ASCII control characters from \em str
     * using \c iscntrl().
     *
     * \b INPUT:
     * \param str  string.
     *
     * \b OUTPUT:
     * \return     string.
     */

    if (str != NULL && *str != 0) {
        char *p, *q,       /* pointers for looping through the input string */
             *out = (char *)malloc(strlen(str) + 1);   /* new output string */

        /* set pointers to beginning of input string */
        q = str;
        p = out;

        /* loop through the input string */
        while (*q) {
            if (iscntrl(*q) && *q != '\0') {
                q++;             /* if control character found jump to next */
            }

            /* overwrite original string with non-control chars */
            *p = *q;
            p++;
            q++;
        }

        *p = '\0';                                      /* terminate string */

        return(out);
    } else {
        return(NULL);
    }
}


void sm_basic_rmcntrl_inplace(char *str)
{
    /*!
     * \brief
     *   Remove control characters from string.
     *
     * This function removes all \c ASCII control characters from \em str
     * using \c iscntrl(). In contrast to \c sm_basic_rmcntrl(), the action is
     * performed in place, i.e. the input gets overwritten.
     *
     * \b INPUT & OUTPUT:
     * \param str  string
     *
     * \b OUTPUT:
     * \return     nothing
     */

    if (str != NULL && *str != 0) {
        char *p, *q; /* pointers for looping through the input string */

        /* set pointers to beginning of input string */
        p = str;
        q = str;

        /* loop through the input string */
        while (*q) {
            if (iscntrl(*q) && *q != '\0') {
                q++;             /* if control character found jump to next */
            }

            /* overwrite original string with non-control chars */
            *p = *q;
            p++;
            q++;
        }

        *p = '\0';                                      /* terminate string */
    }
}


char *sm_basic_strtrim(char *str)
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


void sm_basic_strtrim_inplace(char *str)
{
    /*!
     * \brief
     *   Remove leading and trailing blanks from string.
     *
     * This function removes all leading and trailing " " characters from
     * \em str using \c isspace(). In contrast to \c sm_basic_strtrim(), the
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


void sm_basic_terminatestring(char *str)
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


cpl_error_code sm_basic_interpollin(const double *x_out, double *y_out,
                                    const long n_out, const double *x_ref,
                                    const double *y_ref, const long n_ref)
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
     *   Points outside the range of the reference vectors will be
     *   extrapolated based on the first / last values in the reference
     *   vectors for the low / high end, respectively.
     *
     * \b INPUT:
     * \param x_out  desired output x-spacing
     * \param n_out  length of \em x_out / \em y_out
     * \param x_ref  reference x-spacing
     * \param y_ref  reference y-values at \em x_ref
     * \param n_ref  length of \em x_ref / \em y_ref
     *
     * \b OUTPUT:
     * \param y_out  requested y-values at \em x_out
     *
     * \return       CPL_ERROR_NONE on success,
     *               CPL_ERROR_DIVISION_BY_ZERO else
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
            xr = x_ref + n_ref - 2;
            yr = y_ref + n_ref - 2;
        }

        /* calculate linear interpolation */
        if (*(xr + 1) - *xr == 0) {
            *yo = *yr;
            divide_by_zero = 1;
        } else {
            *yo = (*xo - *xr) * (*(yr + 1) - *yr) / (*(xr + 1) - *xr) + *yr;
        }
    }

    if (divide_by_zero == 1) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_DIVISION_BY_ZERO,
                                     "Duplicate x-values");
    } else {
        return CPL_ERROR_NONE;
    }
}

/**@}*/
