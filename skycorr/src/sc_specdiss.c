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
 * \file sc_specdiss.c
 *
 * \brief programme to find emission lines in spectra
 *
 * \author Wolfgang Kausch, Stefan Noll, & ESO In-Kind Team Innsbruck
 * \since  01 Feb 2011
 * \date   17 Sep 2013
 *
 */


/*****************************************************************************
 *                                  INCLUDES                                 *
 ****************************************************************************/

#include <sc_specdiss.h>


/*****************************************************************************
 *                                  DEFINES                                  *
 ****************************************************************************/

#define EXIST F_OK        /* check for file existence,  unistd.h,  access() */
#define EXECUTABLE X_OK   /* check for execute permissions                  */


/*****************************************************************************
 *                                    CODE                                   *
 ****************************************************************************/

cpl_error_code sc_specdiss_find_emissionlines(cpl_table *input_spectrum,
                                              cpl_table *linetab,
                                              cpl_parameterlist *parlist,
                                              const cpl_table *groups)
{
    /*!
     * \callgraph
     *
     * This routine is aimed to detect line peaks in an input spectrum to
     * separate pixels belonging to lines and to the continuum. The lines are
     * finally separated in isolated/blended lines, the first being used to
     * determine a rough estimate of the FWHM of the spectrum.
     *
     * The input spectrum has to be a CPL table containing (at least) the
     * following columns:
     *
     *      col #1: lambda
     *      col #2: flux
     *
     * The wavelength dimension must be [micron], the flux is arbitrary. An
     * additional column labelled "class" is added to the spectrum table,
     * which is used to classify the single pixels in the following way:
     *
     *      0 = Continuum
     *      1 = Line pixel
     *      2 = Line peak
     *      3 = Isolated line peak
     *
     * The second input parameter must be an empty CPL table, but already
     * initialised with sc_specdiss_init_linetab(). This table is filled with
     * information about the detected lines
     * (see ::sc_specdiss_find_isolatedlines).
     *
     * \b INPUT:
     * \param input_spectrum  table containing input spectra
     * \param linetab         empty table to be filled with line information
     * \param parlist         CPL parameter list containing information of the
     *                        parameter file
     * \param groups          CPL table with airglow line information
     *
     * \b OUTPUT:
     *
     * \b ERRORS:
     * - CPL_ERROR_NONE:           no error occurred
     * - CPL_ERROR_ILLEGAL_INPUT:  missing column "lambda" or "flux" in input
     *                             spectrum, zero length input spectrum
     */

    cpl_error_code err_code=CPL_ERROR_NONE;
    cpl_error_code status = CPL_ERROR_NONE;
    char err_msg[SC_MAXLEN];         /* String containing error message  */

    int in_nrow=0;                      /* Length of input spectrum         */
    int loop=0;                         /* Loop variables                   */
    double value=0.;                    /* Simple spectrum value            */

    double sig = 0.;                    /* sigma of 1st derivative          */
    double mean = 0.;                   /* mean of abs. val. of 1st deriv.  */
    double maxrat = 4.0;                /* max. sig/mean for peaks          */
    double deriv1_lim=0.;               /* 1st derivative peak criterion    */

    float lim_fac=1;                    /* threshold: lim_fac * deriv1_lim  */
    double flux_low_lim=0., flux_up_lim=0.; /* limits to search for         */

    /* Histogram stuff */
    //cpl_table *hist1;               /* CPL table,  2cols 'bins',  'counts'*/
    //int n_bins=0;                     /* number of bins to be used        */
    //int bin_scale=30; /* # of histogram bins n_bins = #rows / bin_scale   */
                        /* This approach is required to obtain a roughly    */
                        /* equal sampling of the histogram                  */
                        /* NOTE: n_bins is an int!                          */
    //char column_name[SC_MAXLEN];      /* Column name to be used for hist  */
    char col1[SC_MAXLEN], col2[SC_MAXLEN];       /* col name string         */
    char plot_title[SC_MAXLEN];

    cpl_table *derivatives;       /* Containing derivatives of inspec       */
    cpl_table *tmptable;
    //cpl_parameterlist *plottags;
    cpl_parameter *pptr;

/*--------------------------------------------------------------------------*/
/*------------------------------- INITIALISING -----------------------------*/
/*--------------------------------------------------------------------------*/

    /* checking error state */

    /* Checking input spectrum */
    if ( (cpl_table_has_column(input_spectrum, "lambda")!= 1 ) ||
         (cpl_table_has_column(input_spectrum, "flux")  != 1 ) )
    {
        sprintf(err_msg, "Either column \"lambda\" or \"flux\" missing in "
                        "input spectrum CPL table.\n"
                        "Cannot continue. Emergency stop.");
        cpl_msg_info(cpl_func, "-------------------------------------------");
        cpl_msg_info(cpl_func, "Incorrect input spectrum CPL table "
                              "structure:");
        cpl_table_dump_structure(input_spectrum, NULL);
        return cpl_error_set_message(cpl_func,  CPL_ERROR_ILLEGAL_INPUT,
                                    "%s",  err_msg);
    }
    in_nrow=cpl_table_get_nrow(input_spectrum);
    if (in_nrow == 0)
    {
        sprintf(err_msg, "Input spectrum has length zero.\n"
                        "Cannot continue. Emergency stop.");
        cpl_msg_info(cpl_func, "-------------------------------------------");
        cpl_msg_info(cpl_func, "Incorrect input spectrum CPL table "
                              "structure:");
        cpl_table_dump_structure(input_spectrum, NULL);
        return cpl_error_set_message(cpl_func,  CPL_ERROR_ILLEGAL_INPUT,
                                    "%s",  err_msg);
    }

/*--------------------------------------------------------------------------*/
/* ------------------------- Stats + Derivatives ---------------------------*/
/*--------------------------------------------------------------------------*/

    /* Creating histogram of spectral values */
    //n_bins=ceil(in_nrow/bin_scale);
    //hist1=cpl_table_new(n_bins);    /* first derivative  */

    //cpl_table_new_column(hist1,  "bins",  CPL_TYPE_DOUBLE);
    //cpl_table_new_column(hist1,  "counts",  CPL_TYPE_INT);

    /* Creating derivatives table */
    derivatives=cpl_table_new(in_nrow);
    cpl_table_duplicate_column(derivatives, "lambda",
                               input_spectrum, "lambda");
    cpl_table_duplicate_column(derivatives, "flux",
                               input_spectrum, "flux");
    cpl_table_duplicate_column(derivatives, "weight",
                               input_spectrum, "weight");
    cpl_table_new_column(derivatives,  "first",  CPL_TYPE_DOUBLE);
    cpl_table_fill_column_window(derivatives, "first", 0, in_nrow, 0.);
    //cpl_table_new_column(derivatives,  "second",  CPL_TYPE_DOUBLE);

    /* Calculating first derivative:                                        */
    /* Defined simply as the flux difference between two consecutive pixels:*/
    /*      \Delta(flux)/pixel                                              */
    /* This simple approach turned out to be the most efficient method      */
    cpl_table_set_double(derivatives, "first", 0, 0.);
    for(loop=0;loop<in_nrow-1;loop++)
    {
        if (cpl_table_get_double(input_spectrum, "weight", loop, NULL)
            != 0. &&
            cpl_table_get_double(input_spectrum, "weight", loop+1, NULL)
            != 0. &&
            cpl_table_get_double(input_spectrum, "lambda", loop, NULL)
            <= SC_THERMIRLIM)
        {
            value=cpl_table_get_double(input_spectrum, "flux", loop+1, NULL)-
                  cpl_table_get_double(input_spectrum, "flux", loop, NULL);
            cpl_table_set_double(derivatives, "first", loop+1, value);
        }
    }

    /* Creating hists/stats on input spectrum/derivatives */
    //sprintf(column_name, "%s", "first");
    //sc_specdiss_create_hist(hist1, derivatives, column_name, n_bins);
    //deriv1_lim=cpl_table_get_column_stdev(derivatives, "first");
    cpl_table_unselect_all(derivatives);
    cpl_table_or_selected_double(derivatives, "first", CPL_NOT_EQUAL_TO, 0.);
    tmptable = cpl_table_extract_selected(derivatives);
    cpl_table_select_all(derivatives);
    cpl_table_abs_column(tmptable, "first");
    sig=cpl_table_get_column_stdev(tmptable, "first");
    mean=cpl_table_get_column_mean(tmptable, "first");
    cpl_table_delete(tmptable);
    if (sig / mean > maxrat) {
        deriv1_lim = maxrat * mean;
    } else {
        deriv1_lim = sig;
    }

    /* Creating histograms */
//     plottags=cpl_parameterlist_new();
//     sc_setplottags_hist(plottags, "first derivative", "{/Symbol D}flux",
//                         "counts per bin", parlist);
//     sc_plot_hist(hist1, plottags);
//     cpl_parameterlist_delete(plottags);

/*--------------------------------------------------------------------------*/
/* ------------------------------ Line search ------------------------------*/
/*--------------------------------------------------------------------------*/

    /* Searching for line peaks                                             */
    flux_low_lim=-lim_fac*deriv1_lim;
    flux_up_lim=lim_fac*deriv1_lim;
    sprintf(col1, "%s", "first");
    sprintf(col2, "%s", "peaks");
    sc_specdiss_find_localmaxima(derivatives, col2, col1, flux_low_lim,
                                 flux_up_lim);

    /* Column "range" in CPL table derivatives denote the pixels regions    */
    /* belonging to lines (line pixels: value = 1)                          */
    /* --> copying to input spectrum table as column "class"                 */
    /* Column "peaks" in CPL table derivatives denote the line peaks        */
    /* (line peaks: value = 2)                          */

    if (cpl_table_has_column(input_spectrum, "class") == 0) {
        cpl_table_duplicate_column(input_spectrum, "class", derivatives,
                                   "range");
    } else {
        cpl_table_copy_data_int(input_spectrum, "class",
                          cpl_table_get_data_int_const(derivatives, "range"));
    }

    for(loop=0;loop<in_nrow;loop++)
    {
        if(cpl_table_get_int(derivatives, "peaks", loop, NULL)==2)
        {
            cpl_table_set_int(input_spectrum, "class", loop, 2);
        }
    }

    /* Checking for airglow lines                   */
    /* Read + merging linelist and input spectrum   */
    status = sc_specdiss_merge_speclinelist(input_spectrum, groups, parlist);
    if (status != CPL_ERROR_NONE) {
        //cpl_table_delete(hist1);
        cpl_table_delete(derivatives);
        return status;
    }

    /* Identify airglow lines in spectrum */
    status = sc_specdiss_identify_airglowlines(input_spectrum);

//     /* Plotting spectrum with identified line peaks */
//     plottags=cpl_parameterlist_new();
//     sc_setplottags_single_spec_lines(plottags, "SINFONI",*/
//                                 "{/Symbol l} [micron]", "Flux", parlist);
//     sc_plot_single_spec_with_lines(input_spectrum, plottags);
//     cpl_parameterlist_delete(plottags);

    /* Identifying isolated line peaks for FWHM determination with the
       "class" column */
    sc_specdiss_find_isolatedlines(linetab, input_spectrum, parlist);

    /* Plotting spectrum with identified line peaks */
    pptr = cpl_parameterlist_find(parlist, "output_name");
    sprintf(plot_title,"%s",cpl_parameter_get_string(pptr));

//     plottags=cpl_parameterlist_new();
//     sc_setplottags_single_spec_lines(plottags, plot_title,
//                                 "{/Symbol l} [micron]", "Flux", parlist);
//     sc_plot_single_spec_with_lines(input_spectrum, plottags);
//     cpl_parameterlist_delete(plottags);

    /* Delete temporary column in input spectrum */
    cpl_table_erase_column(input_spectrum, "linelist");

    /* Cleaning */
    //cpl_table_delete(hist1);
    cpl_table_delete(derivatives);

    return err_code;
}

/*--------------------------------------------------------------------------*/

cpl_error_code sc_specdiss_find_isolatedlines(cpl_table *linetab,
                                              cpl_table *input_spectrum,
                                              cpl_parameterlist *parlist)
{
    /*!
     * \callgraph
     *
     * This programme searches for isolated lines (required for determining
     * the FWHM) in an input spectrum.
     *
     * The following two criteria to identify an isolated line are applied:
     *
     *      (1) Distance to neighbouring line: larger / equal
     *              [min_line_dist_fac] * [measured FWHM] (in px);
     *          For first line: = distance to spectrum start (in px),
     *              last  line: = distance to spec end
     *      (2) Symmetry of line shape: # pixels left of peak = # pixels on
     *          the right hand side (+/- 1 px)
     *
     * The input spectrum must be a CPL table containing (at least) the
     * following columns:
     *
     *      col #1: lambda
     *      col #2: flux
     *      col #3: class
     *
     * The column "class" finally contains information on the pixel type
     * in following way:
     *
     *      0 = Continuum
     *      1 = Line pixel
     *      2 = Line peak
     *      3 = Isolated line peak (used for FWHM determination)
     *
     *
     * Additionally, the CPL table linetab is filled with the following
     * information:
     *
     *      col #1:  width          - line width [pixels]
     *      col #2:  peak_loc       - Peak location [pixels]
     *      col #3:  peak_lam       - Peak location [wavelength]
     *      col #4:  peak_flux      - Flux of line peak pixel
     *      col #5:  px_bef         - number of line pixels before line peak
     *      col #6:  px_aft         - number of line pixels after line peak
     *      col #7:  line_px_start  - first pixel of line
     *      col #8:  line_px_end    - last pixel of line
     *      col #9:  dist_bef       - distance to previous line [pixel]
     *      col #10: dist_aft       - distance to subsequent line [pixels]
     *      col #11: isol_flag      - flag indicating an isol. line
     *                                (0=no / 1=yes)
     *      col #12: fwhm           - Full Width Half Maximum column (to be
     *                                filled later)
     *      col #13: line_list      - Line list value:
     *                                      =0 : not present
     *                                      >0 : # of linelist lines in pixel
     *
     * \note In the case of a wavelength-dependent line width (varfwhm = 1),
     *       the measured fwhm and peak_flux are modified in the way as the
     *       line would be moved to the centre of the wavelength range of the
     *       input spectrum.
     *
     * \b INPUT:
     * \param linetab         table to be filled with information on the lines
     * \param input_spectrum  table containing input spectra
     * \param parlist         CPL parameter list containing information of the
     *                        parameter file
     *
     * \b OUTPUT:
     * \param linetab         information about the detected lines is stored
     *                        to this table
     * \param input_spectrum  spectrum table col "class" is filled
     *
     * \b ERRORS:
     * - CPL_ERROR_NONE:      no error occurred
     */

    cpl_error_code err_code=CPL_ERROR_NONE;

    int in_nrow=0;              /* Number of columns / rows in input spec   */
    int loop=0;                 /* Loop variables                           */
    int iteration=0;            /* # of current iteration                   */
    int peak_flag=0, peak_loc=0; /* Flag if line peak is found within line  */
                                 /* peak_loc=location of line peak in pix!  */
    int line_width=0, px_bef_peak=0, px_aft_peak=-1;
    int line_dist=0, line_count=0;

    int line_list_val=0;
    double min_line_dist=0,     /* Minimum distance between two lines       */
           min_line_dist_0=0;

    double fwhm=0;              /* initial FWHM                             */
    int varfwhm=0;              /* wavelength-dependent FWHM?               */
    double meanlam=0;           /* mean wavelength [micron]                 */

    int linetabsize=0;
    int line_start_flag=0, line_px_start=0;
    double peak_lambda=0., peak_flux=0.;
    int found=0, index=0;       /* found flag / index for writing isol lines*/
    int dist_bef=0, dist_aft=0, symm=0; /* used for definition of criteria  */

    cpl_parameter *p;           /* CPL parameters used to read parlist      */

/*--------------------------------------------------------------------------*/
/* ------------------------------ INITIALISING -----------------------------*/
/*--------------------------------------------------------------------------*/

    /* Message of iteration */
    p = cpl_parameterlist_find(parlist,"iteration");
    iteration = cpl_parameter_get_int(p);
    iteration++;
    cpl_parameter_set_int(p,iteration);

    /* Initial FWHM */
    p = cpl_parameterlist_find(parlist, "fwhm");
    fwhm = cpl_parameter_get_double(p);

    /* Variable FWHM (linear change with wavelength)? */
    p = cpl_parameterlist_find(parlist, "varfwhm");
    varfwhm = cpl_parameter_get_int(p);

    /* Get mean wavelength */
    p = cpl_parameterlist_find(parlist, "meanlam");
    meanlam = cpl_parameter_get_double(p);

    /* Minimal line dist */
    p = cpl_parameterlist_find(parlist,"min_line_dist_fac");
    min_line_dist_0  = fwhm * cpl_parameter_get_double(p);

    /* Creating histogram of spectral values */
    in_nrow=cpl_table_get_nrow(input_spectrum);

    line_dist=0;
    for(loop=1;loop<in_nrow-1;loop++)
    {
        /* Initialising search over line region */
        line_width=0;
        px_bef_peak=0;    /* # of pixels before/after the line peak         */
        px_aft_peak=-1;   /* default -1, because px of peak is counted      */
        peak_loc=0;
        line_start_flag=0;
        line_list_val=0;
        line_dist++;

        /* Investigating region with a line (==value 1 in col "class")       */
        while(cpl_table_get_int(input_spectrum, "class", loop, NULL)>=1)
        {
            /* Determining first pixel of line */
            if (line_start_flag==0)
            {
                line_start_flag=1;
                line_px_start=loop;
            }
            line_width++;
            found=1;
            /* Checking location of line peak */
            peak_flag=cpl_table_get_int(input_spectrum, "class", loop,
                                        NULL);
            if(peak_flag==2)
            {
                peak_loc=loop;
                peak_lambda=cpl_table_get_double(input_spectrum, "lambda",
                                                 loop, NULL);
                peak_flux=cpl_table_get_double(input_spectrum, "flux",
                                                 loop, NULL);
                line_list_val=cpl_table_get_int(input_spectrum, "linelist",
                                                 loop, NULL);

            }
            /* Counting px before/after peak for symmetry criterion         */
            if(peak_loc==0)      /* peak found if peak_loc > 0 */
            {
                px_bef_peak++;
            }
            else
            {
                px_aft_peak++;
            }
            loop++;
        }

        /* Correct peak_flux if FWHM is wavelength-dependent */
        if (varfwhm == 1) {
            peak_flux *= peak_lambda / meanlam;
        }

        /* Filling line property table */
        if (found==1 && peak_flux > 0.)
        {
            linetabsize++;
            cpl_table_set_size(linetab, linetabsize);
            cpl_table_set_int(linetab, "width", linetabsize-1, line_width);
            cpl_table_set_int(linetab, "peak_loc", linetabsize-1, peak_loc);
            cpl_table_set_double(linetab, "peak_lam", linetabsize-1,
                                 peak_lambda);
            cpl_table_set_double(linetab, "peak_flux", linetabsize-1,
                                 peak_flux);
            cpl_table_set_int(linetab, "px_bef", linetabsize-1, px_bef_peak);
            cpl_table_set_int(linetab, "px_aft", linetabsize-1, px_aft_peak);
            cpl_table_set_int(linetab, "line_px_start", linetabsize-1,
                                 line_px_start);
            cpl_table_set_int(linetab, "line_px_end", linetabsize-1, loop-1);
            if (iteration == 1) {
                cpl_table_set_double(linetab, "fwhm", linetabsize-1, 0.);
            } else {
                cpl_table_set_double(linetab, "fwhm", linetabsize-1, fwhm);
            }
            cpl_table_set_int(linetab, "line_list", linetabsize-1,
                              line_list_val);
            found=0;
            line_dist=0;
        }
    }

    /* Search for isolated lines:
     *
     * Applied criteria:
     *
     *      (1) Distance to neighbouring line >=
     *              <min_line_dist_fac> * <measured FWHM>   [pixel]
     *          For first line: = distance to spectrum start (in px),
     *              last  line: = distance to spec end
     *      (2) # pixels left of peak = # pixels on the right hand side
     *          (+/- 1 px)
     *
     */

    for(loop=0;loop<linetabsize;loop++)
    {
        peak_lambda=cpl_table_get_double(linetab, "peak_lam", loop, NULL);

        /* First isolated line: Distance before=dist to spectrum start      */
        if(loop==0)
        {
            dist_bef=cpl_table_get_int(linetab, "line_px_start", loop, NULL);
            cpl_table_set_int(linetab, "dist_bef", loop, dist_bef);
        }
        else
        {
            dist_bef=cpl_table_get_int(linetab, "line_px_start", loop, NULL)-
                     cpl_table_get_int(linetab, "line_px_end", loop-1, NULL);
            cpl_table_set_int(linetab, "dist_bef", loop, dist_bef);
        }

        if(loop==linetabsize-1)
        {
            dist_aft=in_nrow-cpl_table_get_int(linetab, "line_px_end", loop,
                                               NULL);
            cpl_table_set_int(linetab, "dist_aft", loop, dist_aft);
        }
        else
        {
            dist_aft=cpl_table_get_int(linetab,"line_px_start", loop+1, NULL)-
                     cpl_table_get_int(linetab, "line_px_end", loop, NULL);
            cpl_table_set_int(linetab, "dist_aft", loop, dist_aft);
        }

        symm=abs(cpl_table_get_int(linetab, "px_bef", loop, NULL)-
                                   cpl_table_get_int(linetab, "px_aft", loop,
                                   NULL));

        /* Get min_line_dist for given wavelength depending on varfwhm flag */
        if (varfwhm == 1) {
            min_line_dist = min_line_dist_0 * peak_lambda / meanlam;
        } else {
            min_line_dist = min_line_dist_0;
        }

        /* Applying criteria */
        if( ( dist_bef   >= min_line_dist ) &&
            ( dist_aft   >= min_line_dist ) &&
            ( symm <= 1 ) )
        {
            index = cpl_table_get_int(linetab, "peak_loc", loop, NULL);
            cpl_table_set_int(input_spectrum, "class", index, 3);
            cpl_table_set_int(linetab, "isol_flag", loop, 1);
            line_count++;
        }
        else
        {
            cpl_table_set_int(linetab, "isol_flag", loop, 0);
        }
    }

    return err_code;
}

/*--------------------------------------------------------------------------*/

/* Internal routines */

cpl_error_code sc_specdiss_find_localmaxima(cpl_table *input_table,
                                            char *newcolname,
                                            char *colname,
                                            const double lower_threshold,
                                            const double upper_threshold)
{
    /*!
     * \callgraph
     *
     * Routine to search for local maxima below and above two thresholds
     * within a column table. It is used to identify line peaks within a range
     * of pixels belonging to a single line.
     *
     * \b INPUT:
     * \param input_table      table containing input spectra
     * \param newcolname       table column to be added to input_table
     * \param colname          name of column to be used for search
     * \param lower_threshold  lower threshold for search
     * \param upper_threshold  upper threshold for search
     *
     * \b OUTPUT:
     * \param input_table      an additional column is added to the input
     *                         table
     * \param newcolname       table column to be added to input_table
     *                         containing search results
     *
     * \b ERRORS:
     * - CPL_ERROR_NONE:           no error occurred
     */

    int nrow=0;             /* # of rows of input table                     */
    int loop=0;             /* Used for loopings                            */
    double value=0.;        /* Table value                                  */
    double loc_max_val=0.;  /* temporary local maximum                      */
    int found_flag=0;
    int loc_max=0;
    char col1[SC_MAXLEN], col2[SC_MAXLEN];

    cpl_table *temporary;   /* internal temporary table containing info not */
                            /* to be stored permanently                     */
    int range_flag=0;

/*--------------------------------------------------------------------------*/

    nrow=cpl_table_get_nrow(input_table);
    temporary=cpl_table_new(nrow);

    /* all entries fullfilling threshold conditions */
    cpl_table_new_column(temporary, "found", CPL_TYPE_INT);
    cpl_table_new_column(temporary, "range", CPL_TYPE_INT);
    cpl_table_new_column(temporary, "range_red", CPL_TYPE_INT);
    cpl_table_new_column(temporary, "range_blue", CPL_TYPE_INT);

    /* Creating new col to be added to input table & initially filled with 0*/
    cpl_table_new_column(input_table, newcolname, CPL_TYPE_INT);
    cpl_table_fill_column_window_int(input_table, newcolname, 0, nrow, 0);

    /* Searching for values */
    cpl_table_duplicate_column(temporary, "target_col", input_table, colname);
/*    sc_specdiss_find_valuesabove(temporary, "above_thresh",
                        "target_col", upper_threshold);*/
    sprintf(col1, "%s", "target_col");
    sprintf(col2, "%s", "found");

    /* Searching for pixels above/below the thresh (= +1/-1, 0 otherwise)   */
    sc_specdiss_find_valuesoutside(temporary, col2, col1, lower_threshold,
                                   upper_threshold);

    /* Searching for line range (=range of pixels belonging to one line)    */
    /* Line range is limited by a change of the values in column "found"    */
    /* from -1 --> 0 or 1 and 1 --> 0 or -1                                 */

    /* Applying sweep technique to avoid false line range identifications   */
    /* caused by asymmetric lines. Sweep technique = search direction first */
    /* towards red end, towards end afterwards blue                         */

    /* Run towards red end of spectrum */
    for (loop=0;loop<nrow-1;loop++)
    {
        found_flag=cpl_table_get_int(temporary, "found", loop, NULL);
        if ( found_flag == 1 )
        {
            range_flag=1;
        }
        if ((cpl_table_get_int(temporary, "found", loop, NULL) ==-1 &&
             cpl_table_get_int(temporary, "found", loop+1, NULL) == 0) ||
            (cpl_table_get_int(temporary, "found", loop, NULL) ==-1 &&
             cpl_table_get_int(temporary, "found", loop+1, NULL) == 1))
        {
            range_flag=0;
        }
        cpl_table_set_int(temporary, "range_red", loop, range_flag);
    }

    /* Run towards blue end of spectrum */
    for ( loop = nrow-2; loop > 0; loop-- )
    {
        found_flag = cpl_table_get_int(temporary, "found", loop, NULL);
        if ( found_flag == -1 )
        {
            range_flag = 1;
        }
        if (( cpl_table_get_int(temporary, "found", loop+1, NULL) == 1 &&
              cpl_table_get_int(temporary, "found", loop, NULL) == 0 ) ||
            ( cpl_table_get_int(temporary, "found", loop+1, NULL) == 1 &&
              cpl_table_get_int(temporary, "found", loop, NULL) == -1 ))
        {
            range_flag = 0;
        }
        cpl_table_set_int(temporary, "range_blue", loop, range_flag);
    }
    cpl_table_set_int(temporary, "range_blue", loop, 0); /* 1st row         */

    /* Merging the two line searches                                        */
    cpl_table_new_column(input_table, "range", CPL_TYPE_INT);

    for (loop=0;loop<nrow;loop++)
    {
        if (( cpl_table_get_int(temporary, "range_red",  loop, NULL) == 1 ) &&
            ( cpl_table_get_int(temporary, "range_blue", loop, NULL) == 1 ) &&
            ( cpl_table_get_double(input_table, "lambda", loop, NULL)
              <= SC_THERMIRLIM ))
        {
            cpl_table_set_int(input_table, "range", loop, 1);
        }
        else
        {
            cpl_table_set_int(input_table, "range", loop, 0);
        }
    }

    /* Searching for local maxima in line ranges */
    loop=0;
    loc_max_val=0;

    while (loop < nrow)
    {
        range_flag = cpl_table_get_int(input_table, "range", loop, NULL);
        if (range_flag == 0)
        {
            loc_max=0;
            loc_max_val=0;
        }
        else
        {
            while ( cpl_table_get_int(input_table, "range", loop, NULL) == 1)
            {
                value=cpl_table_get_double(input_table, "flux", loop, NULL);
                if (value > loc_max_val)
                {
                    loc_max=loop;
                    loc_max_val=value;
                }
                loop++;
            }
            loop--;
            cpl_table_set_int(input_table, newcolname, loc_max, 2);
        }
        loop++;
    }

    cpl_table_delete(temporary);
    return CPL_ERROR_NONE;

}

/*--------------------------------------------------------------------------*/

cpl_error_code sc_specdiss_find_valuesoutside(cpl_table *input_table,
                                              char *newcolname,
                                              char *colname,
                                              const double lowlim,
                                              const double uplim)
{
    /*!
     * Routine to search for entries in a table column fulfilling the
     * following criteria:
     *
     *      entry <= lowlim && entry >= uplim
     *
     * An additional column (2nd input parameter) is added to the input table,
     * rows fulfilling the condition are flagged:
     *
     *      value = -1,  if entry <= lowlim
     *      value =  1,  if entry >= uplim
     *      value =  0,  otherwise
     *
     * \b INPUT:
     * \param input_table  table containing input spectra
     * \param newcolname   column to be added to input table containing search
     *                     results
     * \param lowlim       lower threshold
     * \param uplim        upper threshold
     * \param colname      column to be used for search
     *
     * \b OUTPUT:
     * \param input_table  an additional column is added to the input table
     * \param newcolname   table column to be added to input_table containing
     *                     search results
     *
     * \b ERRORS:
     * - CPL_ERROR_NONE:   no error occurred
     *
     */

    int nrow=0;         /* # of rows of input table                         */
    int loop=0;         /* Used for loopings                                */
    double value=0.;    /* Table value                                      */

    nrow=cpl_table_get_nrow(input_table);

    for(loop=0;loop<nrow;loop++)
    {
        cpl_table_set_int(input_table, newcolname, loop, 0); /* default=0   */
        value=cpl_table_get(input_table, colname, loop, NULL);
        if( (value <= lowlim) )
        {
            cpl_table_set_int(input_table, newcolname, loop, -1);
        }
        if ((value >= uplim) )
        {
            cpl_table_set_int(input_table, newcolname, loop, 1);
        }

    }

    return CPL_ERROR_NONE;

}

/*--------------------------------------------------------------------------*/

int sc_specdiss_count_lines(const cpl_table *intable)
{
    /*!
     * Routine to count detected line peaks,  i.e. counting the # of rows
     * flagged with a value > 0 in column "peaks" of intable
     *
     * \b INPUT:
     * \param intable  input CPL table
     *
     * \b OUTPUT:
     * \return         number of detected line peaks
     */

    int n_lines=0;
    int loop=0;
    int in_nrow=0;

    in_nrow=cpl_table_get_nrow(intable);

    for(loop=0;loop<in_nrow-1;loop++)
    {
        n_lines=n_lines+cpl_table_get_int(intable, "peaks", loop, NULL);
    }
    return n_lines;
}

/*--------------------------------------------------------------------------*/

cpl_error_code sc_specdiss_create_hist(cpl_table *hist,
                                       const cpl_table *input_table,
                                       char *col_name,  const int n_bins)
{
    /*!
     * \callgraph
     *
     *  Routine to create a histogram binning with n_bins from a choosable
     *  column of an input CPL table. The histogram data is stored into a CPL
     *  table
     *
     * \b INPUT:
     * \param hist         CPL table to be filled with histogram information
     * \param input_table  table containing column "col_name"
     * \param col_name     column of input table to be used for binning
     * \param n_bins       # of histogram bins
     *
     * \b OUTPUT:
     * \param hist         CPL table with histogram information (bins/counts)
     *
     * \b ERRORS:
     * - CPL_ERROR_NONE:   no error occurred
     */

    int tab_nrow=0;                     /* length of table                  */
    double col_min=0., col_max=0.;      /* min/max value of column          */
    int loop=0, loop2=0;                /* for loopings                     */
    int count=0;
    double bin_loop=0.;
    double data_point=0.;
    double lowlim=0., uplim=0.;

    cpl_table *bin_limits;

/*--------------------------------------------------------------------------*/

    tab_nrow=cpl_table_get_nrow(input_table);

    col_min=cpl_table_get_column_min(input_table, col_name);
    col_max=cpl_table_get_column_max(input_table, col_name);

    /* Creating equal n_bins histogram bins between col_min and col_max     */
    bin_limits=cpl_table_new(n_bins);
    cpl_table_new_column(bin_limits, "lowlim", CPL_TYPE_DOUBLE);
    cpl_table_new_column(bin_limits, "uplim", CPL_TYPE_DOUBLE);

    for (loop=0;loop<=n_bins-1;loop++)
    {
        bin_loop=col_min+loop*(col_max-col_min)/n_bins;
        cpl_table_set_double(bin_limits, "lowlim", loop, bin_loop);
        bin_loop=col_min+(loop+1)*(col_max-col_min)/n_bins;
        cpl_table_set_double(bin_limits, "uplim", loop, bin_loop);
    }
    /* Creating hist */
    for (loop=0;loop<=n_bins-1;loop++)
    {
        lowlim=cpl_table_get_double(bin_limits, "lowlim", loop, NULL);
        uplim=cpl_table_get_double(bin_limits, "uplim", loop, NULL);

        /* Searching for values in bin */
        count=0;
        for (loop2=0;loop2<=tab_nrow-1;loop2++)
        {
            data_point=cpl_table_get(input_table, col_name, loop2, NULL);
            if( ( data_point >= lowlim ) && ( data_point < uplim ) )
            {
                count++;
            }
        }
    cpl_table_set(hist, "bins", loop, lowlim+(uplim-lowlim)/2);
    cpl_table_set(hist, "counts", loop, count);
    }

    /* Cleaning */
    cpl_table_delete(bin_limits);
    return CPL_ERROR_NONE;
}

/*--------------------------------------------------------------------------*/

void sc_specdiss_init_linetab(cpl_table *linetab)
{
    /*!
     * Routine to initialise CPL table containing all information about
     * the detected lines. The following columns are added:
     *
     *      col #1:  width          - line width [pixels]
     *      col #2:  peak_loc       - Peak location [pixels]
     *      col #3:  peak_lam       - Peak location [wavelength]
     *      col #4:  peak_flux      - Flux of line peak pixel
     *      col #5:  px_bef         - number of line pixels before line peak
     *      col #6:  px_aft         - number of line pixels after line peak
     *      col #7:  line_px_start  - first pixel of line
     *      col #8:  line_px_end    - last pixel of line
     *      col #9:  dist_bef       - distance to previous line [pixel]
     *      col #10: dist_aft       - distance to subsequent line [pixels]
     *      col #11: isol_flag      - flag indicating an isol. line
     *                                (0=no / 1=yes)
     *      col #12: fwhm           - Full Width Half Maximum column
     *      col #13: line_list      - Line list value:
     *                                      =0 : not present
     *                                      >0 : # of linelist lines in pixel
     *
     *
     * \b INPUT:
     * \param linetab  empty CPL table
     *
     * \b OUTPUT:
     * \param linetab  table with initialised columns for information on
     *                 detected lines
     */

    cpl_table_new_column(linetab,  "width",  CPL_TYPE_INT);
    cpl_table_new_column(linetab,  "peak_loc",  CPL_TYPE_INT);
    cpl_table_new_column(linetab,  "peak_lam",  CPL_TYPE_DOUBLE);
    cpl_table_new_column(linetab,  "peak_flux",  CPL_TYPE_DOUBLE);
    cpl_table_new_column(linetab,  "px_bef",  CPL_TYPE_INT);
    cpl_table_new_column(linetab,  "px_aft",  CPL_TYPE_INT);
    cpl_table_new_column(linetab,  "line_px_start",  CPL_TYPE_INT);
    cpl_table_new_column(linetab,  "line_px_end",  CPL_TYPE_INT);
    cpl_table_new_column(linetab,  "dist_bef",  CPL_TYPE_INT);
    cpl_table_new_column(linetab,  "dist_aft",  CPL_TYPE_INT);
    cpl_table_new_column(linetab,  "isol_flag",  CPL_TYPE_INT);
    cpl_table_new_column(linetab,  "fwhm",  CPL_TYPE_DOUBLE);
    cpl_table_new_column(linetab,  "line_list",  CPL_TYPE_INT);
}

/*--------------------------------------------------------------------------*/

cpl_error_code sc_specdiss_merge_speclinelist(cpl_table *spec,
                                              const cpl_table *groups,
                                             const cpl_parameterlist *parlist)
{
    /*!
     * Routine to read line list and merge line list with input spectrum. The
     * line list is read from a file, all lines of the list are marked
     * in the column "linelist" of the input spectrum table in the
     * following way:
     *
     *      1 = line list pixel
     *      2 = line list peak
     *
     * \b INPUT:
     * \param spec     CPL table with input spectrum
     * \param groups   CPL table with airglow line information
     * \param parlist  parameter list containing file name of line list
     *
     * \b OUTPUT:
     * \param spec     CPL table with initialised columns for information on
     *                 detected lines
     */


    double fwhm = 0, fwhm_crit = 0, fwhm_crit_fac = 0.1, /* FWHM criterion  */
           fwhm_crit_0 = 0.,
           line_flux_sum = 0.,
           pxscl = 0.,                          /* pixel scale [micron/px]  */
           meanlam = 0.;                        /* mean lambda [micron]     */

    double lam_spec0 = 0., lam_spec1 = 0.,      /* neighbouring criterion   */
           lam_spec2 = 0., lam_list = 0.,       /* wavelength in line list  */
           lam_lim1 = 0., lam_lim2 = 0.;


    int listlength=0;          /* length of original line list              */
    int line_arr_size=0;       /* length of selected line list              */

    int i=0, j=0, k=0,         /* loop variables                            */
        n_loops = 10;          /* maximum number of loops/iterations        */

    int peakcount1 = 2;        /* # lines of linelist in pixel              */

    int varfwhm = 0;           /* wavelength-dependent FWHM?                */

    cpl_array *line_arr;       /* line array                                */

    const cpl_parameter *p;    /* CPL parameter                             */

/*--------------------------------------------------------------------------*/

    /* Initialising                                                         */

    /* Initialising line array */
    listlength = cpl_table_get_nrow(groups);
    n_loops = cpl_table_get_nrow(spec);
    if (cpl_table_has_column(spec, "linelist") != 1)
    {
        cpl_table_new_column(spec,  "linelist",  CPL_TYPE_INT);
    }

    /* Initialising new column */
    for (i = 0;i < n_loops; i++)
    {
        cpl_table_set_int(spec,"linelist", i, 0);
    }

    /* Calculating FWHM                                                     */
    p = cpl_parameterlist_find_const(parlist, "fwhm");
    fwhm = cpl_parameter_get_double(p);

    pxscl = (cpl_table_get_double(spec,"lambda",n_loops-1,NULL) -
             cpl_table_get_double(spec,"lambda",0,NULL) ) / n_loops;
    fwhm_crit_0 = fwhm_crit_fac * fwhm * pxscl;

    /* Variable FWHM (linear change with wavelength)? */
    p = cpl_parameterlist_find_const(parlist, "varfwhm");
    varfwhm = cpl_parameter_get_int(p);

    /* Get mean wavelength */
    p = cpl_parameterlist_find_const(parlist, "meanlam");
    meanlam = cpl_parameter_get_double(p);

    /* Search & count lines */
    peakcount1=0;
    /* Counting linelist lines per spectrum pixel */
    for (i = 1; i < n_loops-1; i++) {
        const double * pdlam;
        cpl_size nbad;

        /* Reset peak counters */
        peakcount1=0;

        /* Counting lines from line list belonging to pixel #i              */
        /* Get wavelength borders of pixels #i-1, #i, #i+1                  */
        lam_spec0 = cpl_table_get_double(spec, "lambda", i-1, NULL);
        lam_spec1 = cpl_table_get_double(spec, "lambda", i,   NULL);
        lam_spec2 = cpl_table_get_double(spec, "lambda", i+1, NULL);

        /* Half distance to neighboured pixels                              */
        lam_lim1 = lam_spec1 - (lam_spec1 - lam_spec0) / 2;
        lam_lim2 = lam_spec1 + (lam_spec2 - lam_spec1) / 2;

        /* Get fwhm_crit for given wavelength depending on varfwhm flag */
        if (varfwhm == 1) {
            fwhm_crit = fwhm_crit_0 * lam_spec1 / meanlam;
        } else {
            fwhm_crit = fwhm_crit_0;
        }

        /* Creating list of lines within current pixel */
        line_arr_size = 0;
        line_arr = cpl_array_new(line_arr_size, CPL_TYPE_DOUBLE);
        line_flux_sum = 0;
        /* TODO probably better to use selections here */
        pdlam = cpl_table_get_data_double_const(groups, "lambda");
        nbad = cpl_table_count_invalid(groups, "lambda");

        /* Looping over linelist */
        for (j = 0; j < listlength; j++)
        {
            /* get current wavelength in line list */
            if (nbad == 0) {
                lam_list = pdlam[j];
            }
            else {
                lam_list = cpl_table_get_double(groups, "lambda", j, NULL);
            }

            if (lam_list < lam_spec0)
            {
                /* skip lines at low end of wavelength range in spec */
                continue;
            }

            /* Counting the lines from the line list  belonging to pixel #i */
            /* Selection criterion: distance closer to pixel #i than        */
            /* to pixels #i-1 or #i+1                                       */
            if ( ( lam_list > lam_lim1 ) && ( lam_list <= lam_lim2 ) )
            {
                    line_arr_size++;
                    cpl_array_set_size(line_arr,line_arr_size);
                    cpl_array_set_double(line_arr,line_arr_size-1,lam_list);
                    line_flux_sum += cpl_table_get_double(groups,
                                                          "flux",j,NULL);

                    peakcount1++;
            }
        }

        /* Checking significance of lines                                   */
        /* Check distance between lines: If 2 or more lines are too closely */
        /* neighboured (delta_lambda < fwhm_crit) --> assumed as being a    */
        /* single one                                                       */
        if (line_arr_size > 1)
        {
            k=0;
            while (k < line_arr_size-1)
            {
                if ( cpl_array_get_double(line_arr,k+1,NULL) -
                     cpl_array_get_double(line_arr,k,NULL) < fwhm_crit )
                {
                    peakcount1 --;
                }
                k++;
            }

        }

        cpl_array_delete(line_arr);
        cpl_table_set_int(spec, "linelist", i, peakcount1);
    }

    return CPL_ERROR_NONE;

}

/*--------------------------------------------------------------------------*/

cpl_error_code sc_specdiss_identify_airglowlines(cpl_table *spec)
{
    /*!
     * Identifies airglow lines in spectrum and updates 'linelist' column
     *
     * \b INPUT:
     * \param spec  CPL table with spectrum
     *
     * \b OUTPUT:
     * \param spec  spectrum table with updated 'linelist' column
     */

    int nrow=0;
    int loop=0;

    nrow = cpl_table_get_nrow(spec);

    for (loop = 0; loop < nrow; loop++)
    {
        if( ( cpl_table_get_int(spec, "class", loop, NULL) >= 2 ) &&
            ( cpl_table_get_int(spec, "linelist", loop, NULL) == 5 ) )
            {
                cpl_table_set_int(spec, "linelist", loop, 3);
            }
    }

    return CPL_ERROR_NONE;

}

/*--------------------------------------------------------------------------*/

cpl_error_code sc_specdiss_find_valuesabove(cpl_table *input_table,
                                            char *newcolname,
                                            char *colname,
                                            const double lowlim)
{
    /*!
     *  Routine to search for values within a column table:
     *      lowlim <= VALUES
     *
     * An additional column (2nd input parameter) is added to the input table,
     * rows fulfilling the condition are flagged (value=1, =0 otherwise).
     *
     * \b INPUT:
     * \param input_table  table containing input spectra
     * \param newcolname   column to be added to input table containing search
     *                     results
     * \param lowlim       lower threshold
     * \param colname      column to be used for search
     *
     * \b OUTPUT:
     * \param input_table  an additional column is added to the input table
     * \param newcolname   table column to be added to input_table containing
     *                     search results
     *
     * \b ERRORS:
     * - CPL_ERROR_NONE:   no error occurred
     *
     */

    int nrow=0;         /* # of rows of input table */
    int loop=0;         /* Used for loopings        */
    double value=0.;    /* Table value              */

    nrow=cpl_table_get_nrow(input_table);
    cpl_table_new_column(input_table, newcolname, CPL_TYPE_INT);

    for(loop=0;loop<nrow;loop++)
    {
        cpl_table_set_int(input_table, newcolname, loop, 0); /* default=0   */
        value=cpl_table_get(input_table, colname, loop, NULL);
        if( (value>=lowlim) )
        {
            cpl_table_set_int(input_table, newcolname, loop, 1);
        }
    }

    return CPL_ERROR_NONE;

}

/**@}*/



