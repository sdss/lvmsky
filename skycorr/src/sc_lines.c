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
 * \file sc_lines.c
 *
 * Routines for handling the airglow line list
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  17 Jul 2013
 * \date   28 Aug 2013
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <sc_lines.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

cpl_error_code sc_lines(cpl_table *groups, cpl_parameterlist *parlist)
{
    /*!
     * \callgraph
     *
     * Reads airglow line list from an ASCII file and modifies the line
     * fluxes depending on variability class, period of the night, season,
     * solar activity, and zenith distance. Air wavelengths are provided if
     * required.
     *
     * \b INPUT:
     * \param parlist  input CPL parameter list
     *
     * \b OUTPUT:
     * \param groups   CPL table with airglow line information
     * \param parlist  parameter list complemented by parameters related to
     *                 airglow model
     *
     * \b ERRORS:
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_array *varpar;

    /* Read line groups from ASCII file */
    if ((status = sc_lines_readlist(groups, parlist)) != CPL_ERROR_NONE) {
        return status;
    }

    /* Get season and time bins for airglow model from time of observation */
    if ((status = sc_lines_getmodelbins(parlist)) != CPL_ERROR_NONE) {
        return status;
    }

    /* Read scaling factors of line variability classes from ASCII file */
    varpar = cpl_array_new(0, CPL_TYPE_DOUBLE);
    if ((status = sc_lines_readvarpar(varpar, parlist)) != CPL_ERROR_NONE) {
        cpl_array_delete(varpar);
        return status;
    }

    /* Modify line fluxes by multiplying factors depending on the variability
       class and the airmass */
    sc_lines_scalelines(groups, varpar, parlist);

    /* Convert vacuum to air wavelengths in line group list if required */
    if ((status = sc_lines_vactoair(groups, parlist)) != CPL_ERROR_NONE) {
        cpl_array_delete(varpar);
        return status;
    }

    /* Free allocated memory */
    cpl_array_delete(varpar);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_lines_readlist(cpl_table *groups,
                                 const cpl_parameterlist *parlist)
{
    /*!
     * Reads airglow line data from ASCII file in the sysdata/ folder and
     * writes them into a CPL table. Each line of the file must consist of
     * vacuum wavelength (1st), unextincted flux (2nd), atmospheric
     * transmission at zenith (3rd), variability class (4rd), roto-vibrational
     * system (5th), line group A number (6th), and line group B number (7th).
     * Header lines are allowed if they are marked by #.
     *
     * \b INPUT:
     * \param parlist  input CPL parameter list
     *
     * \b OUTPUT:
     * \param groups   CPL table with airglow line information
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     */

    FILE *stream;
    const cpl_parameter *p;
    char basedir[SC_MAXLEN], datadir[SC_MAXLEN];
    char linetabname[SC_MAXLEN], linetabfile[SC_MAXLEN];
    char errtxt[SC_MAXLEN], str[SC_LENLINE+2];
    cpl_boolean isoutrange = CPL_FALSE;
    int nhead = 0, nrec = 0, feat = 0, syst = 0, agroup = 0, bgroup = 0;
    int ncol0 = 0, ncolmin = 7, i = 0, ncol = 0;
    double lam = 0., flux = 0., trans = 0.;

    /* Create CPL table columns */
    cpl_table_new_column(groups, "lambda", CPL_TYPE_DOUBLE);
    cpl_table_new_column(groups, "flux", CPL_TYPE_DOUBLE);
    cpl_table_new_column(groups, "trans", CPL_TYPE_DOUBLE);
    cpl_table_new_column(groups, "feat", CPL_TYPE_INT);
    cpl_table_new_column(groups, "system", CPL_TYPE_INT);
    cpl_table_new_column(groups, "groupA", CPL_TYPE_INT);
    cpl_table_new_column(groups, "groupB", CPL_TYPE_INT);

    /* Open file with sky line data */
    p = cpl_parameterlist_find_const(parlist, "inst_dir");
    strncpy(basedir, cpl_parameter_get_string(p), SC_MAXLEN);
    p = cpl_parameterlist_find_const(parlist, "data_dir");
    sc_basic_abspath(datadir, cpl_parameter_get_string(p), basedir);
    p = cpl_parameterlist_find_const(parlist, "linetabname");
    strncpy(linetabname, cpl_parameter_get_string(p), SC_MAXLEN);
    sprintf(linetabfile, "%s%s", datadir, linetabname);
    if ((stream = fopen(linetabfile, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", SC_ERROR_FOF_TXT, linetabfile);
        return cpl_error_set_message(cpl_func, SC_ERROR_FOF, "%s", errtxt);
    }

    /* Write info message */
    cpl_msg_info(cpl_func, "Read line list %s", linetabfile);

    /* Find number of header and data lines */
    while (fgets(str, SC_LENLINE+2, stream) != NULL) {
        if (str[0] == '#') {
            nhead++;
            if (nrec != 0) {
                fclose(stream);
                sprintf(errtxt, "%s: %s (comment line in data part)",
                        SC_ERROR_UFS_TXT, linetabfile);
                return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s",
                                             errtxt);
            }
        } else if (isdigit(str[0]) || isspace(str[0]) || str[0] == '-') {
            if (nrec == 0) {
                /* Number of values per line */
                ncol0 = sscanf(str, "%le %le %le %d %d %d %d", &lam, &flux,
                               &trans, &feat, &syst, &agroup, &bgroup);
                if (ncol0 == 0) {
                    fclose(stream);
                    sprintf(errtxt, "%s: %s (empty line)",
                        SC_ERROR_UFS_TXT, linetabfile);
                    return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s",
                                                 errtxt);
                } else if (ncol0 < ncolmin) {
                    fclose(stream);
                    sprintf(errtxt, "%s: %s (too low number of columns)",
                        SC_ERROR_UFS_TXT, linetabfile);
                    return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s",
                                                 errtxt);
                }
            }
            nrec++;
        } else {
            fclose(stream);
            sprintf(errtxt, "%s: %s (unexpected first character at line)",
                    SC_ERROR_UFS_TXT, linetabfile);
            return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s",
                                         errtxt);
        }
    }
    rewind(stream);

    /* No data points */
    if (nrec == 0) {
        fclose(stream);
        sprintf(errtxt, "%s: %s (no data)",
                SC_ERROR_UFS_TXT, linetabfile);
        return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s", errtxt);
    }

    /* Resize output table */
    cpl_table_set_size(groups, nrec);

    /* Skip header */
    for (i = 0; i < nhead; i++) {
        if (fgets(str, SC_LENLINE+2, stream)) {};
    }

    /* Read line data from file and write it to CPL table */

    for (i = 0; i < nrec; i++) {

        ncol = fscanf(stream, "%le %le %le %d %d %d %d",
                      &lam, &flux, &trans, &feat, &syst, &agroup, &bgroup);

        if (ncol != ncol0) {
            cpl_table_set_size(groups, 0);
            fclose(stream);
            sprintf(errtxt, "%s: %s (unexpected number of values at line)",
                    SC_ERROR_UFS_TXT, linetabfile);
            return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s",
                                         errtxt);
        }

        cpl_table_set(groups, "lambda", i, lam);
        cpl_table_set(groups, "flux", i, flux);
        cpl_table_set(groups, "trans", i, trans);
        cpl_table_set(groups, "feat", i, feat);
        cpl_table_set(groups, "system", i, syst);
        cpl_table_set(groups, "groupA", i, agroup);
        if (syst > 0 && bgroup == 0) {
            /* Negative numbers for molecular lines without B group number */
            cpl_table_set(groups, "groupB", i, -syst);
        } else {
            cpl_table_set(groups, "groupB", i, bgroup);
        }

        if (flux < 0.) {
            flux = 0.;
            isoutrange = CPL_TRUE;
        }

        if (trans < 0.) {
            trans = 0.;
            isoutrange = CPL_TRUE;
        }

        if (trans > 1.) {
            trans = 1.;
            isoutrange = CPL_TRUE;
        }

        if (feat < 1) {
            feat = 1;
            isoutrange = CPL_TRUE;
        }

        if (agroup < 0) {
            agroup = 0;
            isoutrange = CPL_TRUE;
        }

        if (bgroup < 0) {
            bgroup = 0;
            isoutrange = CPL_TRUE;
        }

    }

    fclose(stream);

    if (isoutrange == CPL_TRUE) {
        cpl_msg_warning(cpl_func, "%s: Input value(s) out of range -> "
                        "Take lowest/highest allowed value(s)",
                        linetabname);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_lines_getmodelbins(cpl_parameterlist *parlist)
{
    /*!
     * Gets season and time bins for airglow model from observing date in
     * years and UT observing time in s and writes them into the general
     * parameter list. Moreover, the solar radio flux is obtained from year
     * and month which are also taken from observing date. The flux value is
     * read from an ASCII file in the data folder. If the required year and
     * month cannot be found, the local file is updated by a one downloaded
     * from the Canadian reference web server. The most recent entry is taken,
     * if the required month is not in the list. If the procedure fails, an
     * average value of 130 sfu is written into the general parameter list.
     * The entire procedure for obtaining the solar radio flux is not carried
     * out if a positive value for the \e solflux parameter is already in the
     * parameter list.
     *
     * \b INPUT:
     * \param parlist  input CPL parameter list
     *
     * \b OUTPUT:
     * \param parlist  parameter list with sky model bins
     *
     * \b ERRORS:
     * - Ivalid object value(s)
     * - see ::sc_lines_readsolflux
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameter *p, *q;
    char errtxt[SC_MAXLEN];
    int year = 0, month = 0, day = 0, hh = 0, mm = 0;
    int season = 0, timebin = 0;
    double ss = 0., fracyear = 0., ut = 0., nlen = 0.;
    double tlim1 = 0., tlim2 = 0.;
    double nstart[6] = {0.98, 0.48, -0.44, -0.57, -0.23, 0.29};
    double nend[6] = {8.50, 9.29, 9.77, 10.07, 9.61, 8.60};

    /* Get date in years and UT time in s from parameter list */
    p = cpl_parameterlist_find(parlist, "date_val");
    fracyear = cpl_parameter_get_double(p);
    p = cpl_parameterlist_find(parlist, "time_val");
    ut = cpl_parameter_get_double(p) / 3600.;

    /* Check for existence of date and time information */
    if (fracyear < 0 || ut < 0) {
        sprintf(errtxt, "%s: cpl_parameterlist *parlist (invalid date_val "
                "and/or time_val)", SC_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }

    /* Derive season bin */
    sc_basic_fracyear2date(&year, &month, &day, &hh, &mm, &ss, &fracyear);
    if (month == 12) {
        season = 1;
    } else {
        season = floor(month / 2) + 1;
    }

    /* Derive time bin */
    if (ut > 16.) {
        ut -= 24.;
    }
    nlen = nend[season-1] - nstart[season-1];
    tlim1 = nstart[season-1] + nlen / 3;
    tlim2 = nend[season-1] - nlen / 3;
    if (ut < tlim1) {
        timebin = 1;
    } else if (ut >= tlim2) {
        timebin = 3;
    } else {
        timebin = 2;
    }

    /* Write season and time bins into parameter list */
    p = cpl_parameter_new_value("season", CPL_TYPE_INT, "", "", season);
    cpl_parameterlist_append(parlist, p);
    p = cpl_parameter_new_value("timebin", CPL_TYPE_INT, "", "", timebin);
    cpl_parameterlist_append(parlist, p);

    /* Write year and month into parameter list */
    p = cpl_parameter_new_value("year", CPL_TYPE_INT, "", "", year);
    cpl_parameterlist_append(parlist, p);
    p = cpl_parameter_new_value("month", CPL_TYPE_INT, "", "", month);
    cpl_parameterlist_append(parlist, p);

    /* Take solar radio flux from parameter list if given */
    q = cpl_parameterlist_find(parlist, "solflux");
    if (cpl_parameter_get_double(q) >= 0) {
        cpl_msg_info(cpl_func, "SOLFLUX = %g sfu",
                     cpl_parameter_get_double(q));
        return CPL_ERROR_NONE;
    }

    /* Get solar radio flux for year and month from file in sysdata/ */
    p = cpl_parameterlist_find(parlist, "soldatsource");
    cpl_parameter_set_string(p, "LOCAL");
    status = sc_lines_readsolflux(parlist);

    /* If reading fails, take flux from file on web server */
    q = cpl_parameterlist_find(parlist, "solflux");
    if (cpl_parameter_get_double(q) < 0) {
        cpl_parameter_set_string(p, "WEB");
        status = sc_lines_readsolflux(parlist);
    }

    /* If the procedure fails, take average value of 130 sfu */
    if (cpl_parameter_get_double(q) < 0) {
        cpl_parameter_set_string(p, "NONE");
        sc_lines_readsolflux(parlist);
    }

    /* Write info message */
    q = cpl_parameterlist_find(parlist, "solflux");
    cpl_msg_info(cpl_func, "SOLFLUX = %g sfu", cpl_parameter_get_double(q));

    return status;
}


cpl_error_code sc_lines_readsolflux(cpl_parameterlist *parlist)
{
    /*!
     * Reads solar radio flux for given year and month from an ASCII file in
     * the sysdata/ folder. If the file source is set to 'WEB' (instead of
     * 'LOCAL'), the local file is substituted by the one from the Canadian
     * reference web site before the reading is started. If the programme
     * fails to find the correct flux value, the \e solflux parameter is set
     * to -1. If no file source is given, the average value of the solar
     * cycles 19 to 23, i.e. 130 sfu, is assumed.
     *
     * \b INPUT:
     * \param parlist  input CPL parameter list
     *
     * \b OUTPUT:
     * \param parlist  CPL parameter list with updated solar radio flux
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     */

    FILE *stream;
    cpl_parameter *p;
    char soldatsource[SC_MAXLEN], soldaturl[SC_MAXLEN];
    char soldatname[SC_MAXLEN], soldatwebfile[SC_MAXLEN];
    char basedir[SC_MAXLEN], datadir[SC_MAXLEN];
    char soldatfile[SC_MAXLEN], sys[SC_MAXLEN];
    char errtxt[SC_MAXLEN], str[SC_LENLINE+2];
    int year = 0, month = 0, i = 0, y = 0, m = 0, ncol = 0, ncolmin = 5;
    double obsflux = 0., adjflux = 0., absflux = 0., solflux = -1.;

    /* Get year and month from parameter list */
    p = cpl_parameterlist_find(parlist, "year");
    year = cpl_parameter_get_int(p);
    p = cpl_parameterlist_find(parlist, "month");
    month = cpl_parameter_get_int(p);

    /* Set solar radio flux to -1 to indicate if month is not found */
    p = cpl_parameterlist_find(parlist, "solflux");
    cpl_parameter_set_double(p, -1.);

    /* Get source of the solar radio flux data from parameter list */
    p = cpl_parameterlist_find(parlist, "soldatsource");
    strncpy(soldatsource, cpl_parameter_get_string(p), SC_MAXLEN);

    /* Take average solar radio flux if no data source */
    if ((strcmp(soldatsource, "LOCAL") != 0 &&
         strcmp(soldatsource, "WEB") != 0) ||
        (strcmp(soldatsource, "WEB") == 0 &&
         ((year == 1947 && month == 1) || year < 1947))) {
        p = cpl_parameterlist_find(parlist, "solflux");
        cpl_parameter_set_double(p, 130.);
        cpl_msg_warning(cpl_func, "No source for solar radio flux available"
                        " -> Take average of 130 sfu");
        return CPL_ERROR_NONE;
    }

    /* Merge web address and name of solar radio data file */
    p = cpl_parameterlist_find(parlist, "soldaturl");
    strncpy(soldaturl, cpl_parameter_get_string(p), SC_MAXLEN);
    p = cpl_parameterlist_find(parlist, "soldatname");
    strncpy(soldatname, cpl_parameter_get_string(p), SC_MAXLEN);
    sprintf(soldatwebfile, "ftp://%s/%s", soldaturl, soldatname);

    /* Get local data folder */
    p = cpl_parameterlist_find(parlist, "inst_dir");
    strncpy(basedir, cpl_parameter_get_string(p), SC_MAXLEN);
    p = cpl_parameterlist_find(parlist, "data_dir");
    sc_basic_abspath(datadir, cpl_parameter_get_string(p), basedir);

    /* Merge local path and name of solar radio data file */
    sprintf(soldatfile, "%s%s", datadir, soldatname);

    /* Download file from Canadian web server if "WEB" is desired as source */
    if (strcmp(soldatsource, "WEB") == 0) {
        cpl_msg_info(cpl_func,
                     "Download solar radio flux data from web server");
        sprintf(sys, "wget -q --directory-prefix=%s %s", datadir,
                soldatwebfile);
        if (system(sys)) {};
        sprintf(sys, "mv %s.1 %s", soldatfile, soldatfile);
        if (system(sys)) {};
    }

    /* Write info message */
    cpl_msg_info(cpl_func, "Get solar radio flux for %d/%d from %s", month,
                 year, soldatfile);

    /* Open file with sky line data */
    if ((stream = fopen(soldatfile, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", SC_ERROR_FOF_TXT, soldatfile);
        return cpl_error_set_message(cpl_func, SC_ERROR_FOF, "%s", errtxt);
    }

    /* Skip header */
    for (i = 0; i < 2; i++) {
        if (fgets(str, SC_LENLINE+2, stream)) {};
    }

    /* Read file and find solar radio flux for given year and month */

    while (fgets(str, SC_LENLINE+2, stream) != NULL) {

        if (isdigit(str[0]) || isspace(str[0])) {

            /* Read values */
            ncol = sscanf(str, "%d %d %le %le %le",
                          &y, &m, &obsflux, &adjflux, &absflux);

            /* Check number of values per line */
            if (ncol == 0) {
                fclose(stream);
                sprintf(errtxt, "%s: %s (empty line)",
                        SC_ERROR_UFS_TXT, soldatfile);
                return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s",
                                             errtxt);
            } else if (ncol < ncolmin) {
                fclose(stream);
                sprintf(errtxt, "%s: %s (too low number of columns)",
                        SC_ERROR_UFS_TXT, soldatfile);
                return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s",
                                             errtxt);
            }

            /* Check validity of values */
            if (y < 1947 || m < 1 || m > 12 || obsflux < 0.) {
                fclose(stream);
                sprintf(errtxt, "%s: %s (invalid value(s))",
                        SC_ERROR_UFS_TXT, soldatfile);
                return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s",
                                             errtxt);
            }

            /* Find solar radio flux for given year and month */
            if (y == year && m == month) {
                solflux = obsflux;
                break;
            }

        } else {

            /* No digit or space at the beginning of the line */
            fclose(stream);
            sprintf(errtxt, "%s: %s (unexpected first character at line)",
                    SC_ERROR_UFS_TXT, soldatfile);
            return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s",
                                         errtxt);

        }

    }

    fclose(stream);

    /* If the required month is not in the file on the web server, take the
       last entry */
    if(strcmp(soldatsource, "WEB") == 0) {
        if ((month != 1 && m < month-1 && y == year) ||
            (month != 1 && y < year) ||
            (month == 1 && m < 12 && y == year-1) ||
            (month == 1 && y < year-1)) {
            /* Warning message if last file entry is outdated */
            cpl_msg_warning(cpl_func, "Solar radio fluxes until %d/%d only",
                            m, y);
        }
        if (obsflux != solflux) {
            solflux = obsflux;
        }
    }

    /* Write solar radio flux into parameter list */
    p = cpl_parameterlist_find(parlist, "solflux");
    cpl_parameter_set_double(p, floor(solflux + 0.5));

    return CPL_ERROR_NONE;
}


cpl_error_code sc_lines_readvarpar(cpl_array *varpar,
                                   cpl_parameterlist *parlist)
{
    /*!
     * Reads data related to airglow scaling from an ASCII file in the
     * sysdata/ folder and fills a CPL array with feature-related correction
     * factors. The resulting factors are related to the standard strengths of
     * airglow emission features of the upper atmosphere given by Hanuschik
     * (2003) and Rousselot et al. (2000). The scaling factors depend on the
     * emission layer width (related to airmass), the monthly-averaged solar
     * flux in sfu, the season of the year, and the time of the day.
     *
     * \b INPUT:
     * \param parlist  input CPL parameter list
     *
     * \b OUTPUT:
     * \param varpar   correction factors for variability classes as CPL array
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     * - Invalid object value(s)
     */

    FILE *stream;
    cpl_parameter *p;
    scpar x[SC_MAXPAR], m[SC_MAXPAR];
    char basedir[SC_MAXLEN], datadir[SC_MAXLEN];
    char vardatname[SC_MAXLEN], vardatfile[SC_MAXLEN];
    char errtxt[SC_MAXLEN];
    int n = SC_MAXPAR, nfeat = 0, nseason = 0, ntime = 0, season = 0;
    int timebin = 0, i = 0, j = 0;
    double alt = 0., solflux = 0., height = 0., scale = 0., cons = 0.;
    double slope = 0., mean = 0., z = 0., cvr = 0., csol = 0.;

    /* Check existence of airglow scaling parameter file */
    p = cpl_parameterlist_find(parlist, "inst_dir");
    strncpy(basedir, cpl_parameter_get_string(p), SC_MAXLEN);
    p = cpl_parameterlist_find(parlist, "data_dir");
    sc_basic_abspath(datadir, cpl_parameter_get_string(p), basedir);
    p = cpl_parameterlist_find(parlist, "vardatname");
    strncpy(vardatname, cpl_parameter_get_string(p), SC_MAXLEN);
    sprintf(vardatfile, "%s%s", datadir, vardatname);
    if ((stream = fopen(vardatfile, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", SC_ERROR_FOF_TXT, vardatfile);
        return cpl_error_set_message(cpl_func, SC_ERROR_FOF, "%s", errtxt);
    }

    /* Write info message */
    cpl_msg_info(cpl_func, "Read line variability model file %s", vardatfile);

    /* Read airglow scaling parameter file */

    /* Read numbers for data set size */
    sc_basic_readline(stream, x, &n);
    nfeat = x[0].i;
    sc_basic_readline(stream, x, &n);
    nseason = x[0].i;
    ntime = x[1].i;
    //nbin = (nseason + 1) * (ntime + 1);

    /* Return NULL vectors in the case of one or more zero for the data set
       size */
    if (nfeat == 0 || nseason == 0 || ntime == 0) {
        fclose(stream);
        sprintf(errtxt, "%s: %s (nfeat == 0 || nseason == 0 || ntime == 0)",
                SC_ERROR_UFS_TXT, vardatfile);
        return cpl_error_set_message(cpl_func, SC_ERROR_UFS, "%s", errtxt);
    }

    /* Get and check airglow model parameters */
    p = cpl_parameterlist_find(parlist, "telalt_val");
    alt = cpl_parameter_get_double(p);
    p = cpl_parameterlist_find(parlist, "solflux");
    solflux = cpl_parameter_get_double(p);
    p = cpl_parameterlist_find(parlist, "season");
    season = cpl_parameter_get_int(p);
    p = cpl_parameterlist_find(parlist, "timebin");
    timebin = cpl_parameter_get_int(p);
    if (alt < 0. || alt > 90. || solflux < 0. ||
        season < 0 || season > nseason || timebin < 0 || timebin > ntime) {
        fclose(stream);
        sprintf(errtxt, "%s: cpl_parameterlist *parlist (telalt, solflux, "
                "season, and/or timebin out of range) ", SC_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }

    /* Set size of varpar array */
    cpl_array_set_size(varpar, nfeat);

    /* Read data for each feature, extract bin data for selected time, and
       fill CPL array with feature-specific flux correction factors */

    for (i = 0; i < nfeat; i++) {

        /* Time-independent parameters */
        sc_basic_readline(stream, x, &n); // skip molmass
        sc_basic_readline(stream, x, &n); // skip temp
        sc_basic_readline(stream, x, &n);
        height = x[0].d;
        sc_basic_readline(stream, x, &n);
        scale = x[0].d;
        sc_basic_readline(stream, x, &n);
        cons = x[0].d;
        slope = x[1].d;

        /* Mean value for selected time */
        for (j = 0; j < ntime + 1; j++) {
            sc_basic_readline(stream, m, &n);
            if (j == timebin) {
                mean = m[season].d;
            }
        }

        /* Standard deviation for selected time (skipped) */
        for (j = 0; j < ntime + 1; j++) {
            sc_basic_readline(stream, m, &n);
        }

        /* Emission line brightness depending on airmass and layer height
           (van Rhijn 1921) */
        z = (90. - alt) * CPL_MATH_RAD_DEG;
        cvr = 1 / sqrt(1 - pow((SC_ERAD / (SC_ERAD + height)) * sin(z), 2));

        /* Influence of solar radio flux [sfu] */
        csol = slope * solflux + cons;

        /* Set feature-specific correction factors */
        cpl_array_set(varpar, i, scale * cvr * csol * mean);

    }

    /* Close file */
    fclose(stream);

    return CPL_ERROR_NONE;
}


cpl_error_code sc_lines_scalelines(cpl_table *groups, const cpl_array *varpar,
                                   const cpl_parameterlist *parlist)
{
    /*!
     * Modifies line fluxes by multiplying factors depending on the
     * variability class. The factors originate from ::sc_lines_readvarpar.
     * The routine also considers the flux reduction due to atmospheric
     * absorption for the given airmass-corrected transmission spectrum (part
     * of the input line table). The airmass is calculated using the given
     * telescope altitude angle and the formula of Rozenberg (1966).
     *
     * \b INPUT:
     * \param groups   CPL table with airglow line information
     * \param varpar   correction factors for variability classes as CPL array
     * \param parlist  input CPL parameter list
     *
     * \b OUTPUT:
     * \param groups   line table with corrected fluxes
     *
     * \b ERRORS:
     * - none
     */

    const cpl_parameter *p;
    cpl_boolean isoutrange = CPL_FALSE;
    int nlin = 0, nfeat = 0, i = 0, feat = 0;
    const double *corr = NULL;
    double z = 0., xz = 0., trans = 0., flux = 0.;

    /* Get number of lines and variability classes */
    nlin = cpl_table_get_nrow(groups);
    nfeat = cpl_array_get_size(varpar);

    /* Get pointer to CPL array of flux correction factors */
    corr = cpl_array_get_data_double_const(varpar);

    /* Get and convert altitude in deg into zenith distance in rad */
    p = cpl_parameterlist_find_const(parlist, "telalt_val");
    z = (90. - cpl_parameter_get_double(p)) * CPL_MATH_RAD_DEG;

    /* Use airmass approximation of Rozenberg (1966) for molecular
       absorption correction */
    xz = 1 / (cos(z) + 0.025 * exp(-11 * cos(z)));

    /* Write info message */
    cpl_msg_info(cpl_func, "Adapt line list");

    /* Correct line fluxes depending on line type and reduce fluxes due to
       airmass-dependent atmospheric absorption */
    for (i = 0; i < nlin; i++) {
        feat = cpl_table_get(groups, "feat", i, NULL);
        if (feat > nfeat) {
            feat = 1;
            isoutrange = CPL_TRUE;
        }
        trans = pow(cpl_table_get(groups, "trans", i, NULL), xz);
        flux = cpl_table_get(groups, "flux", i, NULL);
        cpl_table_set(groups, "flux", i, flux * corr[feat-1] * trans);
    }

    /* Print warning if the read line type is too high
       -> Value is substituted by 1 */
    if (isoutrange == CPL_TRUE) {
        cpl_msg_warning(cpl_func, "Variability class too high -> "
                        "Take class 1");
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sc_lines_vactoair(cpl_table *groups,
                                 const cpl_parameterlist *parlist)
{
    /*!
     * Converts vacuum wavelengths to air wavelengths by using the formula of
     * Edlen (1966). No action is performed if the wavelengths of the input
     * spectrum are already for vacuum.
     *
     * \b INPUT:
     * \param groups   line list with vacuum wavelengths in \f$\mu{\rm m}\f$
     * \param parlist  general CPL parameter list
     *
     * \b OUTPUT:
     * \param groups   line list with vacuum wavelengths in \f$\mu{\rm m}\f$
     *                 if required
     *
     * \b ERRORS:
     * - Invalid object value(s)
     */

    const cpl_parameter *p;
    char vac_air[4], errtxt[SC_MAXLEN];
    int nrow = 0, i = 0;
    double *lam, sig2 = 0., n = 0.;

    /* Input spectrum: vacuum or air wavelengths? */
    p = cpl_parameterlist_find_const(parlist, "vac_air");
    strcpy(vac_air, cpl_parameter_get_string(p));

    /* Check label for wavelength type */
    if (strncmp(vac_air, "vac", 3) != 0 && strncmp(vac_air, "air", 3) != 0) {
        sprintf(errtxt, "%s: cpl_parameterlist *parlist (vac_air neither "
                "'vac' nor 'air')", SC_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SC_ERROR_IOV, "%s", errtxt);
    }

    /* No conversion in the case of vacuum line wavelengths in the spectrum */
    if (strncmp(vac_air, "vac", 3) == 0) {
        return CPL_ERROR_NONE;
    }

    /* Get number of rows */
    nrow = cpl_table_get_nrow(groups);

    /* Get pointer to wavelengths */
    lam = cpl_table_get_data_double(groups, "lambda");

    /* Get refractive index and convert wavelengths */
    for (i = 0; i < nrow; i++) {
         sig2 = pow(lam[i], -2);
         n = 8342.13 + 2406030. / (130. - sig2) + 15997. / (38.9 - sig2);
         n = 1. + 1e-8 * n;
         lam[i] /= n;
    }

    return CPL_ERROR_NONE;
}

/**@}*/
