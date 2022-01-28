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
 * \file sc_par.c
 *
 * Routines for handling parameter files
 *
 * \author  Wolfgang Kausch, Stefan Noll, & ESO In-Kind Team Innsbruck
 * \since   02 Feb 2011
 * \date    28 Aug 2013
 *
 */


/*****************************************************************************
 *                                  INCLUDES                                 *
 ****************************************************************************/

#include <sc_par.h>


/*****************************************************************************
 *                                    CODE                                   *
 ****************************************************************************/

cpl_error_code sc_par_readfile(cpl_parameterlist *parlist,
                               const char *filename)
{
    /*!
     * Initialises general CPL parameter list (if not already done) and reads
     * parameters of this list from a provided driver file.
     *
     * \b INPUT:
     * \param parlist   empty CPL parameter list
     * \param filename  name of parameter file
     *
     * \b OUTPUT:
     * \param parlist   CPL parameter list with driver file parameters
     *
     * \b ERRORS:
     * - File not found
     * - Illegal input
     */

    FILE *parfile;              /* Stream to parameter file                 */
    int ncol=0;                 /* For counting numbers of columns of inspec*/
    int colstrlength=0;
    char tmpname[SC_MAXLEN] = "";    /* temporary name variable             */
    int loop=0;                 /* Used for loops                           */
    char line[SC_MAXLEN];       /* line to read in param file               */
    char err_msg[SC_MAXLEN];    /* error message to be returned             */
    char *entry;                /* string part in param line                */
    char cwd[SC_MAXLEN] = "";   /* path to working directory                */
    char basedir[SC_MAXLEN] = "";    /* installation directory              */
    cpl_parameter *p;

    /*----------------------------------------------------------------------*/

    /*------------ Initialise parameters and set default values ------------*/

    if (cpl_parameterlist_get_size(parlist) == 0) {

        /* Directories and files */
        p = cpl_parameter_new_value("inst_dir", CPL_TYPE_STRING, "", "",
                                    "./");
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("data_dir", CPL_TYPE_STRING, "", "",
                                    "sysdata/");
        cpl_parameterlist_append(parlist, p);  /* not to be changed by user */
        p = cpl_parameter_new_value("config_dir", CPL_TYPE_STRING, "", "",
                                    "config/");
        cpl_parameterlist_append(parlist, p);  /* not to be changed by user */
        p = cpl_parameter_new_value("spectype", CPL_TYPE_STRING, "", "",
                                    "SKY");
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("scispec", CPL_TYPE_STRING, "", "",
                                    "NONE");
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("skyspec", CPL_TYPE_STRING, "", "",
                                    "NONE");
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("output_dir", CPL_TYPE_STRING, "", "",
                                    "output/");
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("output_name", CPL_TYPE_STRING, "", "",
                                    "test");
        cpl_parameterlist_append(parlist, p);

        /* Input structure */
        p = cpl_parameter_new_value("col_lam", CPL_TYPE_STRING, "", "",
                                    "NONE");
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("col_flux", CPL_TYPE_STRING, "", "",
                                    "NONE");
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("col_dflux", CPL_TYPE_STRING, "", "",
                                    "NONE");
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("col_mask", CPL_TYPE_STRING, "", "",
                                    "NONE");
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("default_error", CPL_TYPE_DOUBLE, "", "",
                                    0.01);
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("wlgtomicron", CPL_TYPE_DOUBLE, "", "",
                                    0.);
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("vac_air", CPL_TYPE_STRING, "", "",
                                    "vac");
        cpl_parameterlist_append(parlist, p);

        /* FITS header keywords */
        p = cpl_parameter_new_value("date_key", CPL_TYPE_STRING, "", "",
                                    "MJD-OBS");
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("date_val", CPL_TYPE_DOUBLE, "", "", -1.);
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("mjd", CPL_TYPE_DOUBLE, "", "", -1.);
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("time_key", CPL_TYPE_STRING, "", "",
                                    "TM-START");
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("time_val", CPL_TYPE_DOUBLE, "", "", -1.);
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("telalt_key", CPL_TYPE_STRING, "", "",
                                    "ESO TEL ALT");
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("telalt_val", CPL_TYPE_DOUBLE, "", "",
                                    -1.);
        cpl_parameterlist_append(parlist, p);

        /* Required input data */
        p = cpl_parameter_new_value("linetabname", CPL_TYPE_STRING, "", "",
                                    "skylinegroups");
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("vardatname", CPL_TYPE_STRING, "", "",
                                    "sky_featvar");
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("soldaturl", CPL_TYPE_STRING, "", "",
                   "ftp.geolab.nrcan.gc.ca/data/solar_flux/monthly_averages");
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("soldatname", CPL_TYPE_STRING, "", "",
                                    "solflux_monthly_average.txt");
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("soldatsource", CPL_TYPE_STRING, "", "",
                                    "NONE");
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("solflux", CPL_TYPE_DOUBLE, "", "", -1.);
        cpl_parameterlist_append(parlist, p);

        /* Line identification */
        p = cpl_parameter_new_value("fwhm", CPL_TYPE_DOUBLE, "", "", 5.);
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_range("varfwhm", CPL_TYPE_INT, "", "", 0, 0, 1);
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("meanlam", CPL_TYPE_DOUBLE, "", "", 1.);
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("ltol", CPL_TYPE_DOUBLE, "", "", 0.01);
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("min_line_width_fac", CPL_TYPE_DOUBLE,
                                    "", "", 0.);
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("min_line_dist_fac", CPL_TYPE_DOUBLE,
                                    "", "", 2.5);
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("min_line_flux_fac", CPL_TYPE_DOUBLE,
                                    "", "", 0.);
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("fluxlim", CPL_TYPE_DOUBLE, "", "", -1.);
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("iteration",CPL_TYPE_INT, "", "", 0);
        cpl_parameterlist_append(parlist, p);

        /* Fitting of sky lines */
        p = cpl_parameter_new_value("ftol", CPL_TYPE_DOUBLE, "", "", 1e-3);
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("xtol", CPL_TYPE_DOUBLE, "", "", 1e-3);
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("wtol", CPL_TYPE_DOUBLE, "", "", 1e-3);
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("cheby_max", CPL_TYPE_INT, "", "", 7);
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("cheby_min", CPL_TYPE_INT, "", "", 3);
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("cheby_const", CPL_TYPE_DOUBLE, "", "",
                                    0.);
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("rebintype", CPL_TYPE_INT, "", "", 1);
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("weightlim", CPL_TYPE_DOUBLE, "", "",
                                    0.67);
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("siglim", CPL_TYPE_DOUBLE, "", "", 15.);
        cpl_parameterlist_append(parlist, p);
        p = cpl_parameter_new_value("fitlim", CPL_TYPE_DOUBLE, "", "", 0.);
        cpl_parameterlist_append(parlist, p);

        /* Plotting */
        p = cpl_parameter_new_value("plot_type", CPL_TYPE_STRING, "", "",
                                    "N");
        cpl_parameterlist_append(parlist, p);

    }

    /*----------------- Read parameters and replace values -----------------*/

    parfile=fopen(filename,"r");
    if (parfile == NULL) {
        sprintf(err_msg,"Cannot find parameter file '%s'. \n"
                        "Cannot continue. Emergency stop.",filename);
        return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_FOUND,
                                    "%s", err_msg);
    }

    sc_basic_initstring(line, SC_MAXLEN);

    while (fgets(line, SC_MAXLEN - 1, parfile) != NULL) {

        /* Skip comments and empty lines */
        if ((line[0] == '#') || (strlen(line) == 0) || (line[0] == '\n')) {
            continue;
        }
        entry=strtok(line,"=");
        sc_basic_strtrim_inplace(entry);

        /* +++++ Directories and files +++++ */

        /* Base directory --------------------------------------------------*/
        /* Base directory == installation directory                         */
        if (strcmp(entry,"INST_DIR")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (getcwd(cwd, sizeof(cwd)) == NULL) {
                return cpl_error_set_message(cpl_func, SC_ERROR_GETCWD,
                                             SC_ERROR_GETCWD_TXT);
            }
            sc_basic_abspath(basedir, entry, cwd);
            p = cpl_parameterlist_find(parlist, "inst_dir");
            cpl_parameter_set_string(p,basedir);
            continue;
        }

        /* Filename of input object spectrum -------------------------------*/
        if (strcmp(entry,"INPUT_OBJECT_SPECTRUM")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (strcmp(entry,"")==0)
            {
                fclose(parfile);
                sprintf(err_msg,"No object spectrum given in parameter file "
                                "'%s': \n"
                                "Cannot continue. Emergency stop.",filename);
                return cpl_error_set_message(cpl_func,CPL_ERROR_ILLEGAL_INPUT,
                                            "%s", err_msg);
            }
            else
            {
                p = cpl_parameterlist_find(parlist, "scispec");
                cpl_parameter_set_string(p,entry);
                continue;
            }
        }

        /* Filename of input sky spectrum ----------------------------------*/
        if (strcmp(entry,"INPUT_SKY_SPECTRUM")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (strcmp(entry,"")==0)
            {
                fclose(parfile);
                sprintf(err_msg,"No sky spectrum given in parameter file "
                                "'%s': \n"
                                "Cannot continue. Emergency stop.",filename);
                return cpl_error_set_message(cpl_func,CPL_ERROR_ILLEGAL_INPUT,
                                            "%s", err_msg);
            }
            else
            {
                p = cpl_parameterlist_find(parlist, "skyspec");
                cpl_parameter_set_string(p,entry);
                continue;
            }
        }

        /* Output directory ------------------------------------------------*/
        if (strcmp(entry,"OUTPUT_DIR")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            p = cpl_parameterlist_find(parlist, "output_dir");
            cpl_parameter_set_string(p,entry);
            continue;
        }

        /* Output name space -----------------------------------------------*/
        if (strcmp(entry,"OUTPUT_NAME")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            p = cpl_parameterlist_find(parlist, "output_name");
            cpl_parameter_set_string(p,entry);
            continue;
        }

        /* +++++ Input structure +++++ */

        /* Checking column names -------------------------------------------*/
        if (strcmp(entry,"COL_NAMES")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            /* Counting # of input columns */
            colstrlength=strlen(entry);
            /* Checking whether # columns != 0 */
            if (colstrlength != 0){
                ncol++;
            }
            for (loop=0;loop<=colstrlength;loop++)
            {
                if (entry[loop] == ' '){
                    ncol++;
                }
            }
            /* At least wavelength and fux column required */
            if (ncol != 4){
                fclose(parfile);
                sprintf(err_msg,"Wrong number of column names in parameter "
                                "COL_NAMES, must be = 4. \n"
                                "Cannot continue. Emergency stop.");
                return cpl_error_set_message(cpl_func,CPL_ERROR_ILLEGAL_INPUT,
                                            "%s", err_msg);
            }

            entry=strtok(entry," ");
            sprintf(tmpname,"col_lam");
            p = cpl_parameterlist_find(parlist, tmpname);
            cpl_parameter_set_string(p,entry);
            entry=strtok(NULL," ");
            sprintf(tmpname,"col_flux");
            p = cpl_parameterlist_find(parlist, tmpname);
            cpl_parameter_set_string(p,entry);
            entry=strtok(NULL," ");
            sprintf(tmpname,"col_dflux");
            p = cpl_parameterlist_find(parlist, tmpname);
            cpl_parameter_set_string(p,entry);
            entry=strtok(NULL," ");
            sprintf(tmpname,"col_mask");
            p = cpl_parameterlist_find(parlist, tmpname);
            cpl_parameter_set_string(p,entry);
            continue;
        }

        /* Relative error if no error column is provided -------------------*/
        if (strcmp(entry,"DEFAULT_ERROR")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (sc_basic_isnumber(entry) == CPL_TRUE)
            {
                p = cpl_parameterlist_find(parlist, "default_error");
                cpl_parameter_set_double(p,atof(entry));
            }
            else
            {
                cpl_msg_warning(cpl_func,"Invalid value for parameter "
                                "DEFAULT_ERROR.");
            }
            continue;
        }

        /* Conversion factor for given wavelength unit to micron -----------*/
        if (strcmp(entry,"WLG_TO_MICRON")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (sc_basic_isnumber(entry) == CPL_TRUE)
            {
                p = cpl_parameterlist_find(parlist, "wlgtomicron");
                cpl_parameter_set_double(p,atof(entry));
            }
            else
            {
                cpl_msg_warning(cpl_func,"Invalid value for parameter "
                                "WLG_TO_MICRON.");
            }
            continue;
        }

        /* Wavelengths in vacuum or air ------------------------------------*/
        if (strcmp(entry,"VAC_AIR")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (strcmp(entry,"vac")==0 || strcmp(entry,"air")==0)
            {
                p = cpl_parameterlist_find(parlist, "vac_air");
                cpl_parameter_set_string(p,entry);
            }
            else
            {
                cpl_msg_warning(cpl_func,"Wrong/Missing value for parameter "
                    "VAC_AIR, assume vacuum.");
            }
            continue;
        }

        /* +++++ FITS header keywords +++++ */

        /* FITS keyword for Modified Julian Day (MJD) or date in years -----*/
        if (strcmp(entry,"DATE_KEY")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            p = cpl_parameterlist_find(parlist, "date_key");
            cpl_parameter_set_string(p,entry);
            continue;
        }

        /* Modified Julian Day (MJD) or date in years ----------------------*/
        if (strcmp(entry,"DATE_VAL")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (sc_basic_isnumber(entry) == CPL_TRUE)
            {
                p = cpl_parameterlist_find(parlist, "date_val");
                cpl_parameter_set_double(p,atof(entry));
            }
            else
            {
                cpl_msg_warning(cpl_func,"Invalid value for parameter "
                                "DATE_VAL.");
            }
            continue;
        }

        /* FITS keyword for UTC time ---------------------------------------*/
        if (strcmp(entry,"TIME_KEY")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            p = cpl_parameterlist_find(parlist, "time_key");
            cpl_parameter_set_string(p,entry);
            continue;
        }

        /* UTC time in s ---------------------------------------------------*/
        if (strcmp(entry,"TIME_VAL")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (sc_basic_isnumber(entry) == CPL_TRUE)
            {
                p = cpl_parameterlist_find(parlist, "time_val");
                cpl_parameter_set_double(p,atof(entry));
            }
            else
            {
                cpl_msg_warning(cpl_func,"Invalid value for parameter "
                                "TIME_VAL.");
            }
            continue;
        }

        /* FITS keyword for telescope altitude angle -----------------------*/
        if (strcmp(entry,"TELALT_KEY")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            p = cpl_parameterlist_find(parlist, "telalt_key");
            cpl_parameter_set_string(p,entry);
            continue;
        }

        /* Telescope altitude angle ----------------------------------------*/
        if (strcmp(entry,"TELALT_VAL")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (sc_basic_isnumber(entry) == CPL_TRUE)
            {
                p = cpl_parameterlist_find(parlist, "telalt_val");
                cpl_parameter_set_double(p,atof(entry));
            }
            else
            {
                cpl_msg_warning(cpl_func,"Invalid value for parameter "
                                "TELALT_VAL.");
            }
            continue;
        }

        /* +++++ Required input data +++++ */

        /* Airglow line list -----------------------------------------------*/
        if (strcmp(entry,"LINETABNAME")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            p = cpl_parameterlist_find(parlist, "linetabname");
            cpl_parameter_set_string(p,entry);
            continue;
        }

        /* File for airglow scaling factors --------------------------------*/
        if (strcmp(entry,"VARDATNAME")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            p = cpl_parameterlist_find(parlist, "vardatname");
            cpl_parameter_set_string(p,entry);
            continue;
        }

        /* FTP address for folder with data on solar radio fluxes ----------*/
        if (strcmp(entry,"SOLDATURL")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            p = cpl_parameterlist_find(parlist, "soldaturl");
            cpl_parameter_set_string(p,entry);
            continue;
        }

        /* File with solar radio fluxes ------------------------------------*/
        if (strcmp(entry,"SOLDATNAME")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            p = cpl_parameterlist_find(parlist, "soldatname");
            cpl_parameter_set_string(p,entry);
            continue;
        }

        /* Solar radio flux ------------------------------------------------*/
        if (strcmp(entry,"SOLFLUX")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (sc_basic_isnumber(entry) == CPL_TRUE)
            {
                p = cpl_parameterlist_find(parlist, "solflux");
                cpl_parameter_set_double(p,atof(entry));
            }
            else
            {
                cpl_msg_warning(cpl_func,"Invalid value for parameter "
                                "SOLFLUX.");
            }
            continue;
        }

        /* +++++ Line identification +++++ */

        /* FWHM ------------------------------------------------------------*/
        if (strcmp(entry,"FWHM")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (sc_basic_isnumber(entry) == CPL_TRUE)
            {
                p = cpl_parameterlist_find(parlist, "fwhm");
                cpl_parameter_set_double(p,atof(entry));
            }
            else
            {
                cpl_msg_warning(cpl_func,"Invalid value for parameter "
                                "FWHM.");
            }
            continue;
        }

        /* Variable FWHM (linear increasing with wavelength)?---------------*/
        if (strcmp(entry,"VARFWHM")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (atoi(entry) == 0 || atoi(entry) == 1)
            {
                p = cpl_parameterlist_find(parlist, "varfwhm");
                cpl_parameter_set_int(p,atoi(entry));
            }
            else
            {
                cpl_msg_warning(cpl_func,"Invalid value for parameter "
                                "VARFWHM.");
            }
            continue;
        }

        /* Relative FWHM convergence criterion -----------------------------*/
        if (strcmp(entry,"LTOL")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (CPL_TRUE)
            {
                p = cpl_parameterlist_find(parlist, "ltol");
                cpl_parameter_set_double(p,atof(entry));
            }
            else
            {
                cpl_msg_warning(cpl_func,"Invalid value for parameter "
                                "LTOL.");
            }
            continue;
        }

        /* Minimum factor of distance to neighbouring line (minfac * FWHM) -*/
        if (strcmp(entry,"MIN_LINE_DIST")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (sc_basic_isnumber(entry) == CPL_TRUE)
            {
                p = cpl_parameterlist_find(parlist, "min_line_dist_fac");
                cpl_parameter_set_double(p,atof(entry));
            }
            else
            {
                cpl_msg_warning(cpl_func,"Invalid value for parameter "
                                "MIN_LINE_DIST.");
            }
            continue;
        }

        /* Relative lower flux limit for isolated lines --------------------*/
        if (strcmp(entry,"MIN_LINE_FLUX")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (sc_basic_isnumber(entry) == CPL_TRUE)
            {
                p = cpl_parameterlist_find(parlist, "min_line_flux_fac");
                cpl_parameter_set_double(p,atof(entry));
            }
            else
            {
                cpl_msg_warning(cpl_func,"Invalid value for parameter "
                                "MIN_LINE_FLUX.");
            }
            continue;
        }

        /* Relative lower flux limit for line identification ---------------*/
        if (strcmp(entry,"FLUXLIM")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (sc_basic_isnumber(entry) == CPL_TRUE)
            {
                p = cpl_parameterlist_find(parlist, "fluxlim");
                cpl_parameter_set_double(p,atof(entry));
            }
            else
            {
                cpl_msg_warning(cpl_func,"Invalid value for parameter "
                                "FLUXLIM.");
            }
            continue;
        }

        /* +++++ Fitting of sky lines +++++ */

        /* Relative chi^2 MPFIT convergence criterion ----------------------*/
        if (strcmp(entry,"FTOL")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (CPL_TRUE)
            {
                p = cpl_parameterlist_find(parlist, "ftol");
                cpl_parameter_set_double(p,atof(entry));
            }
            else
            {
                cpl_msg_warning(cpl_func,"Invalid value for parameter "
                                "FTOL.");
            }
            continue;
        }

        /* Relative parameter MPFIT convergence criterion ------------------*/
        if (strcmp(entry,"XTOL")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (CPL_TRUE)
            {
                p = cpl_parameterlist_find(parlist, "xtol");
                cpl_parameter_set_double(p,atof(entry));
            }
            else
            {
                cpl_msg_warning(cpl_func,"Invalid value for parameter "
                                "XTOL.");
            }
            continue;
        }

        /* Relative chi^2 convergence criterion for wavegrid polynomial ----*/
        if (strcmp(entry,"WTOL")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (CPL_TRUE)
            {
                p = cpl_parameterlist_find(parlist, "wtol");
                cpl_parameter_set_double(p,atof(entry));
            }
            else
            {
                cpl_msg_warning(cpl_func,"Invalid value for parameter "
                                "WTOL.");
            }
            continue;
        }

        /* Maximum degree of Chebyshev polynomial for wavegrid correction --*/
        if (strcmp(entry,"CHEBY_MAX")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (sc_basic_isnumber(entry) == CPL_TRUE)
            {
                p = cpl_parameterlist_find(parlist, "cheby_max");
                cpl_parameter_set_int(p,atoi(entry));
            }
            else
            {
                cpl_msg_warning(cpl_func,"Invalid value for parameter "
                                "CHEBY_MAX.");
            }
            continue;
        }

        /* Minimum degree of Chebyshev polynomial for wavegrid correction --*/
        if (strcmp(entry,"CHEBY_MIN")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (sc_basic_isnumber(entry) == CPL_TRUE)
            {
                p = cpl_parameterlist_find(parlist, "cheby_min");
                cpl_parameter_set_int(p,atoi(entry));
            }
            else
            {
                cpl_msg_warning(cpl_func,"Invalid value for parameter "
                                "CHEBY_MIN.");
            }
            continue;
        }

        /* Initial constant term of Chebyshev polynomial -------------------*/
        if (strcmp(entry,"CHEBY_CONST")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (sc_basic_isnumber(entry) == CPL_TRUE)
            {
                p = cpl_parameterlist_find(parlist, "cheby_const");
                cpl_parameter_set_double(p,atof(entry));
            }
            else
            {
                cpl_msg_warning(cpl_func,"Invalid value for parameter "
                                "CHEBY_CONST.");
            }
            continue;
        }

        /* Type of rebinning -----------------------------------------------*/
        if (strcmp(entry,"REBINTYPE")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (sc_basic_isnumber(entry) == CPL_TRUE)
            {
                p = cpl_parameterlist_find(parlist, "rebintype");
                cpl_parameter_set_int(p,atoi(entry));
            }
            else
            {
                cpl_msg_warning(cpl_func,"Invalid value for parameter "
                                "REBINTYPE.");
            }
            continue;
        }

        /* Lower pixel weight limit for line selection ---------------------*/
        if (strcmp(entry,"WEIGHTLIM")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (sc_basic_isnumber(entry) == CPL_TRUE)
            {
                p = cpl_parameterlist_find(parlist, "weightlim");
                cpl_parameter_set_double(p,atof(entry));
            }
            else
            {
                cpl_msg_warning(cpl_func,"Invalid value for parameter "
                                "WEIGHTLIM.");
            }
            continue;
        }

        /* Sigma-clipping limit --------------------------------------------*/
        if (strcmp(entry,"SIGLIM")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (sc_basic_isnumber(entry) == CPL_TRUE)
            {
                p = cpl_parameterlist_find(parlist, "siglim");
                cpl_parameter_set_double(p,atof(entry));
            }
            else
            {
                cpl_msg_warning(cpl_func,"Invalid value for parameter "
                                "SIGLIM.");
            }
            continue;
        }

        /* Selection criterion for line group fitting ----------------------*/
        if (strcmp(entry,"FITLIM")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            if (sc_basic_isnumber(entry) == CPL_TRUE)
            {
                p = cpl_parameterlist_find(parlist, "fitlim");
                cpl_parameter_set_double(p,atof(entry));
            }
            else
            {
                cpl_msg_warning(cpl_func,"Invalid value for parameter "
                                "FITLIM.");
            }
            continue;
        }

        /* +++++ Plotting +++++ */

        /* Checking gnuplot terminal ---------------------------------------*/
        if (strcmp(entry,"PLOT_TYPE")==0)
        {
            entry=strtok(NULL,"=");
            sc_basic_strtrim_inplace(entry);
            p = cpl_parameterlist_find(parlist, "plot_type");
            cpl_parameter_set_string(p,entry);
            continue;
        }

    }

    fclose(parfile);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_par_writefile(const cpl_parameterlist *parlist,
                                const char *filename)
{
    /*!
     * Writes parameters of the provided CPL parameter list to the given
     * driver file. Only CPL parameters with names matching parameter strings
     * in the file are considered.
     *
     * \b INPUT:
     * \param parlist   input CPL parameter list
     * \param filename  name of parameter file
     *
     * \b OUTPUT:
     * - none
     *
     * \b ERRORS:
     * - File opening failed
     */

    FILE *parfile;
    const cpl_parameter *p, *p2, *p3, *p4;
    cpl_array *text;
    char errtxt[SC_MAXLEN], line[SC_LENLINE+2], str[SC_LENLINE+2], *entry;
    int nline = 0, i = 0;

    /* Open parameter file for reading */
    if ((parfile = fopen(filename, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", SC_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, SC_ERROR_FOF, "%s", errtxt);
    }

    /* Get number of file lines */
    while (fgets(line, SC_LENLINE+2, parfile) != NULL) {
        nline++;
    }
    rewind(parfile);

    /* Create CPL string array for file content */
    text = cpl_array_new(nline, CPL_TYPE_STRING);

    /* Read parameter file */
    for (i = 0; i < nline; i++) {
        entry = fgets(line, SC_LENLINE+2, parfile);
        cpl_array_set_string(text, i, line);
    }

    /* Close parameter file */
    fclose(parfile);

    /* Open parameter file for writing */
    if ((parfile = fopen(filename, "w+")) == NULL) {
        cpl_array_delete(text);
        sprintf(errtxt, "%s: %s", SC_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, SC_ERROR_FOF, "%s", errtxt);
    }

    /* Rewrite driver file with updated parameter values */

    for (i = 0; i < nline; i++) {

        /* Get parameter label */
        strcpy(str, cpl_array_get_string(text, i));
        entry = strtok(str, "=");
        sc_basic_strtrim_inplace(entry);

        /* +++++ Directories and files +++++ */

        /* Base directory --------------------------------------------------*/
        if (strcmp(entry, "INST_DIR") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "inst_dir")) != NULL) {
            fprintf(parfile, "%s=%s\n", entry, cpl_parameter_get_string(p));
            continue;
        }

        /* Filename of input object spectrum -------------------------------*/
        if (strcmp(entry, "INPUT_OBJECT_SPECTRUM") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "scispec")) != NULL) {
            fprintf(parfile, "%s=%s\n", entry, cpl_parameter_get_string(p));
            continue;
        }

        /* Filename of input sky spectrum ----------------------------------*/
        if (strcmp(entry, "INPUT_SKY_SPECTRUM") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "skyspec")) != NULL) {
            fprintf(parfile, "%s=%s\n", entry, cpl_parameter_get_string(p));
            continue;
        }

        /* Output directory ------------------------------------------------*/
        if (strcmp(entry, "OUTPUT_DIR") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "output_dir"))
            != NULL) {
            fprintf(parfile, "%s=%s\n", entry, cpl_parameter_get_string(p));
            continue;
        }

        /* Output name space -----------------------------------------------*/
        if (strcmp(entry, "OUTPUT_NAME") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "output_name"))
            != NULL) {
            fprintf(parfile, "%s=%s\n", entry, cpl_parameter_get_string(p));
            continue;
        }

        /* +++++ Input structure +++++ */

        /* Checking column names -------------------------------------------*/
        if (strcmp(entry, "COL_NAMES") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "col_lam")) != NULL &&
            (p2 = cpl_parameterlist_find_const(parlist, "col_flux"))
            != NULL &&
            (p3 = cpl_parameterlist_find_const(parlist, "col_dflux"))
            != NULL &&
            (p4 = cpl_parameterlist_find_const(parlist, "col_mask"))
            != NULL) {
            fprintf(parfile, "%s=%s %s %s %s\n", entry,
                    cpl_parameter_get_string(p),
                    cpl_parameter_get_string(p2),
                    cpl_parameter_get_string(p3),
                    cpl_parameter_get_string(p4));
            continue;
        }

        /* Relative error if no error column is provided -------------------*/
        if (strcmp(entry, "DEFAULT_ERROR") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "default_error"))
            != NULL) {
            fprintf(parfile, "%s=%g\n", entry, cpl_parameter_get_double(p));
            continue;
        }

        /* Conversion factor for given wavelength unit to micron -----------*/
        if (strcmp(entry, "WLG_TO_MICRON") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "wlgtomicron"))
            != NULL) {
            fprintf(parfile, "%s=%g\n", entry, cpl_parameter_get_double(p));
            continue;
        }

        /* Wavelengths in vacuum or air ------------------------------------*/
        if (strcmp(entry, "VAC_AIR") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "vac_air")) != NULL) {
            fprintf(parfile, "%s=%s\n", entry, cpl_parameter_get_string(p));
            continue;
        }

        /* +++++ FITS header keywords +++++ */

        /* FITS keyword for Modified Julian Day (MJD) or date in years -----*/
        if (strcmp(entry, "DATE_KEY") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "date_key")) != NULL) {
            fprintf(parfile, "%s=%s\n", entry, cpl_parameter_get_string(p));
            continue;
        }

        /* Modified Julian Day (MJD) or date in years ----------------------*/
        if (strcmp(entry, "DATE_VAL") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "date_val")) != NULL) {
            fprintf(parfile, "%s=%g\n", entry, cpl_parameter_get_double(p));
            continue;
        }

        /* FITS keyword for UTC time ---------------------------------------*/
        if (strcmp(entry, "TIME_KEY") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "time_key")) != NULL) {
            fprintf(parfile, "%s=%s\n", entry, cpl_parameter_get_string(p));
            continue;
        }

        /* UTC time in s ---------------------------------------------------*/
        if (strcmp(entry, "TIME_VAL") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "time_val")) != NULL) {
            fprintf(parfile, "%s=%g\n", entry, cpl_parameter_get_double(p));
            continue;
        }

        /* FITS keyword for telescope altitude angle -----------------------*/
        if (strcmp(entry, "TELALT_KEY") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "telalt_key"))
            != NULL) {
            fprintf(parfile, "%s=%s\n", entry, cpl_parameter_get_string(p));
            continue;
        }

        /* Telescope altitude angle ----------------------------------------*/
        if (strcmp(entry, "TELALT_VAL") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "telalt_val"))
            != NULL) {
            fprintf(parfile, "%s=%g\n", entry, cpl_parameter_get_double(p));
            continue;
        }

        /* +++++ Required input data +++++ */

        /* Airglow line list -----------------------------------------------*/
        if (strcmp(entry, "LINETABNAME") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "linetabname"))
            != NULL) {
            fprintf(parfile, "%s=%s\n", entry, cpl_parameter_get_string(p));
            continue;
        }

        /* File for airglow scaling factors --------------------------------*/
        if (strcmp(entry, "VARDATNAME") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "vardatname"))
            != NULL) {
            fprintf(parfile, "%s=%s\n", entry, cpl_parameter_get_string(p));
            continue;
        }

        /* FTP address for folder with data on solar radio fluxes ----------*/
        if (strcmp(entry, "SOLDATURL") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "soldaturl"))
            != NULL) {
            fprintf(parfile, "%s=%s\n", entry, cpl_parameter_get_string(p));
            continue;
        }

        /* File with solar radio fluxes ------------------------------------*/
        if (strcmp(entry, "SOLDATNAME") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "soldatname"))
            != NULL) {
            fprintf(parfile, "%s=%s\n", entry, cpl_parameter_get_string(p));
            continue;
        }

        /* Solar radio flux ------------------------------------------------*/
        if (strcmp(entry, "SOLFLUX") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "solflux")) != NULL) {
            fprintf(parfile, "%s=%g\n", entry, cpl_parameter_get_double(p));
            continue;
        }

        /* +++++ Line identification +++++ */

        /* FWHM ------------------------------------------------------------*/
        if (strcmp(entry, "FWHM") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "fwhm")) != NULL) {
            fprintf(parfile, "%s=%g\n", entry, cpl_parameter_get_double(p));
            continue;
        }

        /* Variable FWHM (linear increasing with wavelength)?---------------*/
        if (strcmp(entry, "VARFWHM") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "varfwhm")) != NULL) {
            fprintf(parfile, "%s=%d\n", entry, cpl_parameter_get_int(p));
            continue;
        }

        /* Relative FWHM determination convergence criterion ---------------*/
        if (strcmp(entry, "LTOL") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "ltol")) != NULL) {
            fprintf(parfile, "%s=%g\n", entry, cpl_parameter_get_double(p));
            continue;
        }

        /* Minimum factor of distance to neighbouring line (minfac * FWHM) -*/
        if (strcmp(entry, "MIN_LINE_DIST") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "min_line_dist_fac"))
            != NULL) {
            fprintf(parfile, "%s=%g\n", entry, cpl_parameter_get_double(p));
            continue;
        }

        /* Relative lower flux limit for isolated lines --------------------*/
        if (strcmp(entry, "MIN_LINE_FLUX") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "min_line_flux_fac"))
            != NULL) {
            fprintf(parfile, "%s=%g\n", entry, cpl_parameter_get_double(p));
            continue;
        }

        /* Relative lower flux limit for line identification ---------------*/
        if (strcmp(entry, "FLUXLIM") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "fluxlim"))
            != NULL) {
            fprintf(parfile, "%s=%g\n", entry, cpl_parameter_get_double(p));
            continue;
        }

        /* +++++ Fitting of sky lines +++++ */

        /* Relative chi^2 MPFIT convergence criterion ----------------------*/
        if (strcmp(entry, "FTOL") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "ftol")) != NULL) {
            fprintf(parfile, "%s=%g\n", entry, cpl_parameter_get_double(p));
            continue;
        }

        /* Relative parameter MPFIT convergence criterion ------------------*/
        if (strcmp(entry, "XTOL") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "xtol")) != NULL) {
            fprintf(parfile, "%s=%g\n", entry, cpl_parameter_get_double(p));
            continue;
        }

        /* Relative chi^2 convergence criterion for wavegrid polynomial ----*/
        if (strcmp(entry, "WTOL") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "wtol")) != NULL) {
            fprintf(parfile, "%s=%g\n", entry, cpl_parameter_get_double(p));
            continue;
        }

        /* Maximum degree of Chebyshev polynomial for wavegrid correction --*/
        if (strcmp(entry, "CHEBY_MAX") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "cheby_max"))
            != NULL) {
            fprintf(parfile, "%s=%d\n", entry, cpl_parameter_get_int(p));
            continue;
        }

        /* Minimum degree of Chebyshev polynomial for wavegrid correction --*/
        if (strcmp(entry, "CHEBY_MIN") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "cheby_min"))
            != NULL) {
            fprintf(parfile, "%s=%d\n", entry, cpl_parameter_get_int(p));
            continue;
        }

        /* Initial constant term of Chebyshev polynomial -------------------*/
        if (strcmp(entry, "CHEBY_CONST") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "cheby_const"))
            != NULL) {
            fprintf(parfile, "%s=%g\n", entry, cpl_parameter_get_double(p));
            continue;
        }

        /* Type of rebinning -----------------------------------------------*/
        if (strcmp(entry, "REBINTYPE") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "rebintype"))
            != NULL) {
            fprintf(parfile, "%s=%d\n", entry, cpl_parameter_get_int(p));
            continue;
        }

        /* Lower pixel weight limit for line selection ---------------------*/
        if (strcmp(entry, "WEIGHTLIM") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "weightlim"))
            != NULL) {
            fprintf(parfile, "%s=%g\n", entry, cpl_parameter_get_double(p));
            continue;
        }

        /* Sigma-clipping limit --------------------------------------------*/
        if (strcmp(entry, "SIGLIM") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "siglim")) != NULL) {
            fprintf(parfile, "%s=%g\n", entry, cpl_parameter_get_double(p));
            continue;
        }

        /* Selection criterion for line group fitting ----------------------*/
        if (strcmp(entry, "FITLIM") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "fitlim")) != NULL) {
            fprintf(parfile, "%s=%g\n", entry, cpl_parameter_get_double(p));
            continue;
        }

        /* +++++ Plotting +++++ */

        /* Checking gnuplot terminal ---------------------------------------*/
        if (strcmp(entry, "PLOT_TYPE") == 0 &&
            (p = cpl_parameterlist_find_const(parlist, "plot_type"))
            != NULL) {
            fprintf(parfile, "%s=%s\n", entry, cpl_parameter_get_string(p));
            continue;
        }

        /* Rewrite lines without identified parameters */
        fprintf(parfile, "%s", cpl_array_get_string(text, i));

    }

    /* Close parameter file */
    fclose(parfile);

    /* Free memory space of string array */
    cpl_array_delete(text);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_par_addkeywords(cpl_propertylist *header,
                                  const cpl_parameterlist *parlist)
{
    /*!
     * Write parameters of the general CPL parameter list to a CPL property
     * list.
     *
     * \b INPUT:
     * \param header   input CPL property list
     * \param parlist  general CPL parameter list
     *
     * \b OUTPUT:
     * \param header   extended CPL property list
     *
     * \b ERRORS:
     * - none
     */

    const cpl_parameter *p;
    char skyfile[SC_MAXLEN] = "", *skyspec = NULL;

    /* Add name of reference sky spectrum */
    p = cpl_parameterlist_find_const(parlist, "skyspec");
    strncpy(skyfile, cpl_parameter_get_string(p), SC_MAXLEN);
    skyspec = strrchr(skyfile, '/');
    cpl_propertylist_append_string(header,"ESO SC REFSKY NAME", (skyspec+1));
    cpl_propertylist_set_comment(header, "ESO SC REFSKY NAME",
                                 "Name of reference sky");

    /* Add MJD of reference sky spectrum if given */
    p = cpl_parameterlist_find_const(parlist, "mjd");
    if (cpl_parameter_get_double(p) > 0.) {
        cpl_propertylist_append_double(header,"ESO SC REFSKY MJD",
                                       cpl_parameter_get_double(p));
        cpl_propertylist_set_comment(header, "ESO SC REFSKY MJD",
                                     "MJD [d] for reference sky");
    }

    /* Add date of reference sky spectrum in years */
    p = cpl_parameterlist_find_const(parlist, "date_val");
    cpl_propertylist_append_double(header,"ESO SC REFSKY DATE",
                                   cpl_parameter_get_double(p));
    cpl_propertylist_set_comment(header, "ESO SC REFSKY DATE",
                                 "Date [yr] for reference sky");

    /* Add time of reference sky spectrum in seconds */
    p = cpl_parameterlist_find_const(parlist, "time_val");
    cpl_propertylist_append_double(header,"ESO SC REFSKY TIME",
                                   cpl_parameter_get_double(p));
    cpl_propertylist_set_comment(header, "ESO SC REFSKY TIME",
                                 "Start time [s] for reference sky");

    /* Add solar radio flux */
    p = cpl_parameterlist_find_const(parlist, "solflux");
    cpl_propertylist_append_float(header,"ESO SC SOLRFLUX",
                                  cpl_parameter_get_double(p));
    cpl_propertylist_set_comment(header, "ESO SC SOLRFLUX",
                                 "Solar radio flux [s.f.u.]");

    /* Add estimated FWHM of sky lines */
    p = cpl_parameterlist_find_const(parlist, "fwhm");
    cpl_propertylist_append_float(header,"ESO SC FWHM",
                                  cpl_parameter_get_double(p));
    cpl_propertylist_set_comment(header, "ESO SC FWHM",
                                 "FWHM of sky lines [px]");

    return CPL_ERROR_NONE;
}

/**@}*/
