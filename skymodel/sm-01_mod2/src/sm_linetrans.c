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
 * \ingroup sm_module1a
 */

/**@{*/

/*!
 * \file sm_linetrans.c
 *
 * Routines for preparing a library of airglow line transmission files
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  06 Dec 2011
 * \date   06 Oct 2015
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <sm_linetrans.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

cpl_error_code sm_linetrans(const char *path)
{
    /*!
     * \callgraph
     * Calculates libraries of transmission lists for airglow lines. A time-
     * and a PWV-dependent library created by a radiative transfer code is
     * considered as input. The library structure is conserved. The file names
     * only differ in one letter indicating the file type ('L' instead of 'T'
     * or 'R'). The PWV-dependent library is skipped if the names of the two
     * library structure files in "sm_filenames.dat" are the same.
     *
     * \b INPUT:
     * \param path  path to "sm_filenames.dat"
     *
     * \b OUTPUT:
     * - none
     *
     * \b ERRORS:
     * - see subroutines
     */

    cpl_parameterlist *params;
    cpl_parameter *p;
    cpl_table *classdat, *filelist;
    smparmodel modelpar;
    smspec linetab;
    char linetabfile[FILENAME_MAX];
    int i = 0;

    /* Write input path to CPL parameter list */
    params = cpl_parameterlist_new();
    p = cpl_parameter_new_value("filepath", CPL_TYPE_STRING, "", "", path);
    cpl_parameterlist_append(params, p);

    /* Add required PWV parameter */
    p = cpl_parameter_new_value("pwv", CPL_TYPE_DOUBLE, "", "", -1);
    cpl_parameterlist_append(params, p);

    /* Fill smparmodel structure with crucial paths and file names */
    sm_etc_readfilenames(&modelpar, params);

    /* Read airglow variability file for mol masses and temperatures of
       variability classes */
    classdat = cpl_table_new(0);
    sm_linetrans_getclassdat(classdat, &modelpar);

    /* Read list of airglow emission lines */
    sprintf(linetabfile, "%s/%s", modelpar.datapath, modelpar.linetabname);
    sm_spec_read(&linetab, linetabfile);

    /* Handle the time- and PWV-dependent libraries */

    for (i = 0; i < 2; i++) {

        /* Skip second library if the names of both library structure files
           are the same */
        if (i == 1 && strcmp(modelpar.libstruct2, modelpar.libstruct1) == 0) {
            break;
        }

        /* Read library structure file to get file names of transmission
           spectra */
        filelist = cpl_table_new(0);
        sm_linetrans_getlibstruct(filelist, &modelpar, i);

        /* Compute lists of line transmission */
        sm_linetrans_calclib(filelist, classdat, &linetab);

        /* Free memory */
        cpl_table_delete(filelist);

    }

    /* Free memory */
    cpl_parameterlist_delete(params);
    cpl_table_delete(classdat);
    sm_spec_free(&linetab);

    /* Return error code of last error */
    return cpl_error_get_code();
}


cpl_error_code sm_linetrans_getclassdat(cpl_table *classdat,
                                        const smparmodel *modelpar)
{
    /*!
     * Writes a CPL table of effective mol masses and temperatures for the
     * different airglow variability classes. The required data are taken from
     * the airglow variability file whose path and name are provided by the
     * input ::smparmodel structure.
     *
     * \b INPUT:
     * \param classdat  empty CPL table
     * \param modelpar  sky emission parameters (see typedef of ::smparmodel)
     *
     * \b OUTPUT:
     * \param classdat  CPL table of mol masses and temperatures for the
     *                  different airglow variability classes
     *
     * \b ERRORS:
     * - FOF: File opening failed
     * - UFS: Unexpected file structure
     */

    FILE *stream;
    smparam x[SM_MAXPAR];
    char vardatfile[FILENAME_MAX], errtxt[SM_MAXLEN+1], line[SM_MAXLEN+2];
    int nfeat = 0, i = 0, repeat = 1;

    /* Check existence of airglow scaling parameter file */
    sprintf(vardatfile, "%s/%s", modelpar->datapath, modelpar->vardatname);
    if ((stream = fopen(vardatfile, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", SM_ERROR_FOF_TXT, vardatfile);
        return cpl_error_set_message(cpl_func, SM_ERROR_FOF, "%s", errtxt);
    }

    /* Get number of variability classes from file */
    sm_param_read(stream, x);
    nfeat = x[0].i;

    /* Return NULL vectors in the case of one or more zero for the data set
       size */
    if (nfeat == 0) {
        fclose(stream);
        sprintf(errtxt, "%s: %s (nfeat == 0)",
                SM_ERROR_UFS_TXT, vardatfile);
        return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s", errtxt);
    }

    /* Prepare output table */
    cpl_table_set_size(classdat, nfeat);
    cpl_table_new_column(classdat, "class", CPL_TYPE_INT);
    cpl_table_new_column(classdat, "molmass", CPL_TYPE_DOUBLE);
    cpl_table_new_column(classdat, "temp", CPL_TYPE_DOUBLE);

    /* Read airglow scaling parameter file and write molmass and temperature
       into output CPL table */

    do {

        /* Read next line (return in case of errors) */
        if (fgets(line, SM_MAXLEN+2, stream) == NULL) {
            fclose(stream);
            cpl_table_set_size(classdat, 0);
            sprintf(errtxt, "%s: %s (reading of line failed)",
                    SM_ERROR_UFS_TXT, vardatfile);
            return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s",
                                         errtxt);
        }

        /* Find required parameters of next feature */
        if (strncmp(line, "# molmass", 8) == 0) {
            /* Class number */
            cpl_table_set(classdat, "class", i, i + 1);
            /* Mol mass */
            sm_param_read(stream, x);
            cpl_table_set(classdat, "molmass", i, x[0].d);
            /* Temperature */
            sm_param_read(stream, x);
            cpl_table_set(classdat, "temp", i, x[0].d);
            /* Finish after reading of parameters of last feature */
            if (i == nfeat - 1) {
                repeat = 0;
            }
            /* Next feature */
            i++;
        }

    } while (repeat == 1);

    /* Close file */
    fclose(stream);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_linetrans_getlibstruct(cpl_table *filelist,
                                         const smparmodel *modelpar,
                                         const int libflag)
{
    /*!
     * Produces a CPL table which contains a list of paths and file names of
     * transmission spectra computed by radiative transfer codes and
     * corresponding lists of airglow line transmissions. The file names are
     * built from the information in the library structure file whose path and
     * name are provided by the input ::smparmodel structure. While libflag = 0
     * indicates a time-dependent library, libflag = 1 selects a PWV-dependent
     * one. As input transmission spectra are only those used that have the
     * highest available resolution.
     *
     * \b INPUT:
     * \param filelist  empty CPL table
     * \param modelpar  sky emission parameters (see typedef of ::smparmodel)
     * \param libflag   type of library
     *
     *
     * \b OUTPUT:
     * \param filelist  list of input and output file names as CPL table
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     * - FOF: File opening failed
     * - UFS: Unexpected file structure
     */

    FILE *stream;
    cpl_table *parnames, *parval;
    smparam nameform[SM_MAXPAR], pardef[SM_MAXPAR], val[SM_MAXPAR];
    char errtxt[SM_MAXLEN+1], libstructfile[FILENAME_MAX];
    char filename[FILENAME_MAX], transspec[FILENAME_MAX];
    char linetrans[FILENAME_MAX], label[SM_LENLINE+1], val0[SM_LENLINE+1];
    char infile[FILENAME_MAX], outfile[FILENAME_MAX];
    int i = 0, ascii0 = 48, j = 0, fpos[7] = {-1,-1,-1,-1,-1,-1,-1};
    int flen[7] = {0}, npar = 0, nfile = 0, j0 = 0, jres = 0, h = 0, k = 0;
    int l = 0, nval = 0;

    /* Build path and name of library structure file */
    if (libflag == 0) {
        sprintf(libstructfile, "%s/%s", modelpar->libpath,
                modelpar->libstruct1);
    } else if (libflag == 1) {
        sprintf(libstructfile, "%s/%s", modelpar->libpath,
                modelpar->libstruct2);
    } else {
        sprintf(errtxt, "%s: wrong libflag (only 0 and 1)", SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Check existence of library structure file */
    if ((stream = fopen(libstructfile, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", SM_ERROR_FOF_TXT, libstructfile);
        return cpl_error_set_message(cpl_func, SM_ERROR_FOF, "%s", errtxt);
    }

    /* Read structure of file names */
    sm_param_read(stream, nameform);
    strcpy(filename, nameform[0].c);

    /* Find positions and lengths of the parameter fields in the file name */

    while (filename[i] != '\0' && i < SM_LENLINE) {
        // printf("DEBUG i=%d j=%d filename[i]=%c filename-ascii0=%d\n", i, j, filename[i], filename[i]-ascii0);

        if (filename[i] - ascii0 == j + 1 && fpos[j] < 0) {
            fpos[j] = i;
        }
        if (filename[i] - ascii0 != j + 1 && fpos[j] >= 0) {
            flen[j] = i - fpos[j];
            j++;
            if (j >= 7) {
                break;
            }
            if (filename[i] - ascii0 == j + 1) {
                fpos[j] = i;
            }
        }
        i++;
    }

    /* Number of parameters found */
    npar = j;

    /* 1 to 7 library parameters are allowed */
    if (npar < 1 || npar > 7) {
        fclose(stream);
        sprintf(errtxt, "%s: %s (N_fields < 1 or > 7)", SM_ERROR_UFS_TXT,
                libstructfile);
        return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s", errtxt);
    }

    /* Copy format string to input and output filenames */
    strcpy(transspec, filename);
    strcpy(linetrans, filename);

    /* Create table of parameter names and numbers of values
       (neglect "spectype") */
    parnames = cpl_table_new(npar - 1);
    cpl_table_new_column(parnames, "ID", CPL_TYPE_INT);
    cpl_table_new_column(parnames, "label", CPL_TYPE_STRING);
    cpl_table_new_column(parnames, "N_val", CPL_TYPE_INT);

    /* Create counters for parameter values and set them to 0 */
    cpl_table_new_column(parnames, "pos", CPL_TYPE_INT);
    cpl_table_fill_column_window(parnames, "pos", 0, npar - 1, 0);

    /* Create table for parameter values */
    parval = cpl_table_new(0);

    /* Read parameter format and values and write fitting values in file name
       template */

    for (nfile = 1, j0 = -1, jres = -1, j = 0; j < npar; j++) {

        /* Read parameter format */
        sm_param_read(stream, pardef);

        /* Check parameter description and parameter ID */
        if (pardef[0].n != 3 || pardef[0].i != j + 1) {
            fclose(stream);
            cpl_table_delete(parnames);
            cpl_table_delete(parval);
            if (pardef[0].n != 3) {
                sprintf(errtxt, "%s: %s (# of parameter description values "
                        "!= 3)", SM_ERROR_UFS_TXT, libstructfile);
            } else {
                sprintf(errtxt, "%s: %s (ID != sequence number)",
                        SM_ERROR_UFS_TXT, libstructfile);
            }
            return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s",
                                         errtxt);
        }

        /* Read parameter values */
        sm_param_read(stream, val);

        /* Fill parameter table with parameter name and number of values
           (skip parameter "spectype") and create column for new parameter
           in table for parameter values (also add new rows if required) */
        if (strncmp(pardef[1].c, "spectype", 8) != 0) {
            cpl_table_set_int(parnames, "ID", h, pardef[0].i);
            cpl_table_set_string(parnames, "label", h, pardef[1].c);
            cpl_table_new_column(parval, pardef[1].c, CPL_TYPE_STRING);
            cpl_table_set_int(parnames, "N_val", h, val[0].n);
            if (cpl_table_get_nrow(parval) < val[0].n) {
                cpl_table_set_size(parval, val[0].n);
            }
            nfile *= val[0].n;
            h++;
        }

        if (strncmp(pardef[1].c, "spectype", 8) == 0) {

            /* Check existence of transmission spectra in library structure
               file and find "spectype" parameter number in file name */

            /* Transmission spectra available? */
            for (j0 = -1, k = 0; k < val[0].n; k++) {
                if (val[k].c[0] == 'T') {
                    j0 = j;
                }
            }


            /* Error if no "spectype" value is 'T' (transmission) */
            if (j0 == -1) {
                fclose(stream);
                cpl_table_delete(parnames);
                cpl_table_delete(parval);
                sprintf(errtxt, "%s: %s (no transmission spectra)",
                        SM_ERROR_UFS_TXT, libstructfile);
                return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s",
                                             errtxt);
            }

        } else {

            /* Identify number of resolution parameter */
            if (strncmp(pardef[1].c, "resol", 5) == 0) {
                jres = j;
            }

            /* Fill data table with parameter values */
            for (k = 0; k < val[0].n; k++) {
                cpl_table_set_string(parval, pardef[1].c, k, val[k].c);
            }

        }

        /* Copy first parameter value to file name template */
        for (i = 0; i < flen[j]; i++) {
            if (j == j0) {
                /* Values for parameter "spectype" depend on file type */
                transspec[i+fpos[j]] = 'T';
                linetrans[i+fpos[j]] = 'L';
            } else {
                if (j == jres) {
                    /* Take spectrum with highest resolution as input only */
                    transspec[i+fpos[j]] = val[(val[0].n)-1].c[i];
                } else {
                    transspec[i+fpos[j]] = val[0].c[i];
                }
                linetrans[i+fpos[j]] = val[0].c[i];
            }
        }

    }

    /* Close file stream */
    fclose(stream);

    /* Set columns and rows of file name table */
    cpl_table_new_column(filelist, "transspec", CPL_TYPE_STRING);
    cpl_table_new_column(filelist, "linetrans", CPL_TYPE_STRING);
    cpl_table_set_size(filelist, nfile);

    /* Copy parameter values to file name template */

    for (l = 0; l < nfile; l++) {

        /* Modify filenames */

        for (h = npar - 2; h >= 0; h--) {

            /* Get current position in the value list of the active
               parameter */
            k = cpl_table_get(parnames, "pos", h, NULL);

            /* Substitute required characters in the input and output file
               names */
            j = cpl_table_get_int(parnames, "ID", h, NULL) - 1;
            sprintf(label, "%s", cpl_table_get_string(parnames, "label", h));
            sprintf(val0, "%s", cpl_table_get_string(parval, label, k));
            for (i = 0; i < flen[j]; i++) {
                if (j != jres) {
                    /* Do not change resolution for input spectra
                       (take highest value instead) */
                    transspec[i+fpos[j]] = val0[i];
                }
                linetrans[i+fpos[j]] = val0[i];
            }

        }

        /* Compose name of input and output library files */
        sprintf(infile, "%s/%s", modelpar->libpath, transspec);
        sprintf(outfile, "%s/%s", modelpar->libpath, linetrans);

        /* Write file names to file list */
        cpl_table_set_string(filelist, "transspec", l, infile);
        cpl_table_set_string(filelist, "linetrans", l, outfile);

        /* Change those parameter values that are required for the next
           library spectrum */

        for (h = npar - 2; h >= 0; h--) {

            /* Get number of values and current position in the value list of
               the active parameter */
            nval = cpl_table_get(parnames, "N_val", h, NULL);
            k = cpl_table_get(parnames, "pos", h, NULL);

            /* Modify position in the value list or finish loop */
            if (k == nval - 1 && h > 0) {
                cpl_table_set(parnames, "pos", h, 0);
            } else {
                cpl_table_set(parnames, "pos", h, k + 1);
                break;
            }

        }

    }

    /* Free memory */
    cpl_table_delete(parnames);
    cpl_table_delete(parval);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_linetrans_calclib(const cpl_table *filelist,
                                    const cpl_table *classdat,
                                    const smspec *linetab)
{
    /*!
     * Calculates library of lists of airglow emission line transmissions. The
     * output files are related to corresponding transmission curves computed
     * by radiative transfer codes. The airglow line positions are taken from
     * the given line table. for the widths of the lines Doppler broadening is
     * assumed which depends on effective mol mass and temperature for the
     * different airglow variability classes.
     *
     * \b INPUT:
     * \param filelist  list of input and output file names as CPL table
     * \param classdat  CPL table of mol masses and temperatures for the
     *                  different airglow variability classes
     * \param linetab   table of airglow emission lines (provided as ::smspec
     *                  structure; dflux1 contains feature number)
     *
     * \b OUTPUT:
     * - none
     *
     * \b ERRORS:
     * - NDA: No data
     * - FOF: File opening failed
     * - IFE: Invalid file name extension
     */

    FILE *stream;
    cpl_table *outtab;
    smspec transspec, linetrans, linespec;
    char errtxt[SM_MAXLEN+1], infile[FILENAME_MAX] = "";
    char oinfile[FILENAME_MAX], outfile[FILENAME_MAX];
    int nfile = 0, nfeat = 0, i = 0, j = 0, feat = 0, k = 0, n = 0, l = 0;
    double limlam[2] = {0., 0.}, lam0 = 0., molmass = 0., temp = 0., sig = 0.;
    double llim[2] = {0., 0.}, dlam = 0., dpsig = 0., trans = 0., psum = 0.;
    double psig = 0., psiglim[2] = {0., 0.}, pxsum = 0., x = 0., px = 0.;

    double addlam = 0.01; // extra wavelength range in mum at limits of
                          // line list
    double cons = 9.618e-9; // constant for calculation of Doppler width
    double maxpsig = 6; // maximum integration range for Gaussian in sigma
    double dsig = 0.01; // bin width in sigma for integration of Gaussian

    /* Get number of table rows */
    nfile = cpl_table_get_nrow(filelist);
    nfeat = cpl_table_get_nrow(classdat);

    /* Check input data */
    if (nfile == 0) {
        sprintf(errtxt, "%s: cpl_table *filelist", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    }
    if (nfeat == 0) {
        sprintf(errtxt, "%s: cpl_table *classdat", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    }
    if (linetab->n <= 0) {
        sprintf(errtxt, "%s: smspec *linetab", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    }

    /* Get wavelengths limits of airglow line list */
    limlam[0] = linetab->dat[0].lam - addlam;
    limlam[1] = linetab->dat[(linetab->n)-1].lam + addlam;

    /* Create output CPL table */
    outtab = cpl_table_new(0);

    /* Create line transmission tables for each library transmission
       spectrum */

    for (i = 0; i < nfile; i++) {

        /* Save previous input filename */
        strcpy(oinfile, infile);

        /* Get path and name of library transmission spectrum */
        sprintf(infile, "%s", cpl_table_get_string(filelist, "transspec", i));

        /* Get path and name of output line transmission list */
        sprintf(outfile, "%s", cpl_table_get_string(filelist, "linetrans",
                                                    i));

        /* Write paths and names of input and output files to stdout */
        cpl_msg_info(cpl_func, "%s -> %s", infile, outfile);

        /* Take previous line transmission list if input file does not change
           (change output file name only) */
        if (strcmp(infile, oinfile) == 0) {
            cpl_table_save(outtab, NULL, NULL, outfile, CPL_IO_CREATE);
            continue;
        }

        /* Check file existence */
        if ((stream = fopen(infile, "r")) == NULL) {
            sprintf(errtxt, "%s: %s", SM_ERROR_FOF_TXT, infile);
            return cpl_error_set_message(cpl_func, SM_ERROR_FOF, "%s",
                                         errtxt);
        }
        fclose(stream);

        /* Check file name extension and read LBLRTM/RFM spectrum
           (each file has to contain a wavelength column and three flux
           columns) */
        if (strstr(infile, ".dat") != NULL ||
            strstr(infile, ".ascii") != NULL) {
            sm_spec_readrange(&transspec, infile, limlam, SM_RRSTEP);
        } else if (strstr(infile, ".fits") != NULL ||
                   strstr(infile, ".mt") != NULL) {
            sm_spec_readfitsrange(&transspec, infile, limlam, SM_RRSTEP);
        } else {
            fclose(stream);
            sprintf(errtxt, "%s: %s ('dat', 'ascii', 'fits', or 'mt' only)",
                    SM_ERROR_IFE_TXT, infile);
            return cpl_error_set_message(cpl_func, SM_ERROR_IFE, "%s",
                                         errtxt);
        }

        /* Create line transmission list */
        sm_spec_copy(&linetrans, linetab);
        sm_spec_changetype(&linetrans, 2);

        /* Get effective transmission for each line */

        for (j = 0; j < linetab->n; j++) {

            /* Get line wavelength and ID */
            lam0 = linetab->dat[j].lam;
            feat = (int) linetab->dat[j].dflux1;

            /* Get mol mass and temperature for given class */
            molmass = cpl_table_get(classdat, "molmass", feat - 1, NULL);
            temp = cpl_table_get(classdat, "temp", feat - 1, NULL);

            /* Compute Doppler width (sigma) of line */
            sig = cons * sqrt(temp / molmass) * lam0;

            /* Get relevant wavelength range for line */
            llim[0] = lam0 - maxpsig * sig;
            llim[1] = lam0 + maxpsig * sig;

            /* Extract relevant wavelength range from transmission spectrum */
            sm_spec_extract(&linespec, &transspec, llim);

            /* Get wavelength step */
            dlam = (linespec.dat[1].lam - linespec.dat[0].lam +
                    linespec.dat[(linespec.n)-1].lam -
                    linespec.dat[(linespec.n)-2].lam) / 2;

            /* Get step in sigma */
            dpsig = dlam / sig;

            /* Calculate effective transmission */

            for (trans = 0., psum = 0., k = 0; k < linespec.n; k++) {

                /* Range of Gaussian in sigma units covered by wavelength
                   bin */
                psig = (linespec.dat[k].lam - lam0) / sig;
                psiglim[0] = psig - dpsig / 2;
                psiglim[1] = psig + dpsig / 2;

                /* Number of integration points */
                n = floor((psiglim[1] - psiglim[0]) / dsig + 0.5);

                /* Compute probability sum for wavelength bin */
                for (pxsum = 0., l = 0; l < n; l++) {
                    x = psiglim[0] + ((double) l - 0.5) * dsig;
                    px = exp(-0.5 * pow(x, 2));
                    pxsum += px;
                }
                pxsum *= dsig / CPL_MATH_SQRT2PI;

                /* Compute transmission and probability sum for line */
                trans += pxsum * linespec.dat[k].flux;
                psum += pxsum;

            }

            /* Compute mean transmission for line */
            trans /= psum;

            /* Write result to line transmission list */
            linetrans.dat[j].flux = trans;

            /* Free memory */
            sm_spec_free(&linespec);

        }

        /* Write line transmission list to FITS file (via CPL table) */
        cpl_table_delete(outtab);
        outtab = cpl_table_new(0);
        sm_spec_writecpl(outtab, &linetrans);
        cpl_table_name_column(outtab, "flux", "trans");
        cpl_table_save(outtab, NULL, NULL, outfile, CPL_IO_CREATE);

        /* Free memory */
        sm_spec_free(&transspec);
        sm_spec_free(&linetrans);

    }

    /* Free memory */
    cpl_table_delete(outtab);

    return CPL_ERROR_NONE;
}

/**@}*/
