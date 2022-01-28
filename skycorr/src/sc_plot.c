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
 * \file sc_plot.c
 *
 * plotting library for molecfit
 *
 * \author Wolfgang Kausch & ESO In-Kind Team Innsbruck
 * \since   02 Feb 2011
 * \date    28 Aug 2013
 *
 */


/*****************************************************************************
 *                                  INCLUDES                                 *
 ****************************************************************************/

#include <sc_plot.h>

/*****************************************************************************
 *                                  DEFINES                                  *
 ****************************************************************************/

#define EXIST W_OK          /* check for file existence, unistd.h, access() */

/*****************************************************************************
 *                                 FUNCTIONS                                 *
 ****************************************************************************/

/*****************************************************************************
 *                                    CODE                                   *
 ****************************************************************************/

/*****************************************************************************
 *                                                                           *
 *  Routine to plot a single spectrum                                        *
 *                                                                           *
 ****************************************************************************/

void sc_setplottags_single_spec(cpl_parameterlist *plottags,
                                const char *x_column,
                                const char *y_column,
                                const char *title,
                                const char *x_label,
                                const char *y_label,
                                const cpl_parameterlist *parlist)
{
/*!
 * \callgraph
 *
 * Routine to set plot tags for single spectrum plot.
 * The following input tags are required:
 *
 * \b INPUT:
 * \param plottags  Parameter list to be filled with tags
 * \param x_column  CPL table column name used for the abscissa
 * \param y_column  CPL table column name used for the ordinate
 * \param title     Plot title
 * \param x_label   Label string for x-axis
 * \param y_label   Label string for y-axis
 * \param parlist   Parameter list containing information from the parameter
 *                  file
 * \b OUTPUT:
 * \param plottags  Filled parameter list containing tags
 *
 * Base directory, output directory, output name, and gnuplot terminal type
 * are added automatically from parameter list.
 * NOTE: NULL is not accepted, use " " instead!
 *
 */

    cpl_parameter *p;
    const cpl_parameter *cp;
    char terminal[SC_MAXLEN], basedir[SC_MAXLEN];
    char outdir[SC_MAXLEN], outname[SC_MAXLEN];

    p = cpl_parameter_new_value("plot_title",CPL_TYPE_STRING,"","",title);
    cpl_parameterlist_append(plottags, p);
    p = cpl_parameter_new_value("xlabel",CPL_TYPE_STRING,"","",x_label);
    cpl_parameterlist_append(plottags, p);
    p = cpl_parameter_new_value("ylabel",CPL_TYPE_STRING,"","",y_label);
    cpl_parameterlist_append(plottags, p);
    p = cpl_parameter_new_value("xcol",CPL_TYPE_STRING,"","",x_column);
    cpl_parameterlist_append(plottags, p);
    p = cpl_parameter_new_value("ycol",CPL_TYPE_STRING,"","",y_column);
    cpl_parameterlist_append(plottags, p);

    cp = cpl_parameterlist_find_const(parlist,"plot_type");
    sprintf(terminal,"%s",cpl_parameter_get_string(cp));
    p = cpl_parameter_new_value("terminal",CPL_TYPE_STRING,"","",terminal);
    cpl_parameterlist_append(plottags, p);

    cp = cpl_parameterlist_find_const(parlist, "inst_dir");
    sprintf(basedir,"%s",cpl_parameter_get_string(cp));
    p = cpl_parameter_new_value("inst_dir",CPL_TYPE_STRING,"","",basedir);
    cpl_parameterlist_append(plottags, p);

    cp = cpl_parameterlist_find_const(parlist, "output_dir");
    sprintf(outdir,"%s",cpl_parameter_get_string(cp));
    p = cpl_parameter_new_value("output_dir",CPL_TYPE_STRING,"","",outdir);
    cpl_parameterlist_append(plottags, p);

    cp = cpl_parameterlist_find_const(parlist, "output_name");
    sprintf(outname,"%s",cpl_parameter_get_string(cp));
    p = cpl_parameter_new_value("output_name",CPL_TYPE_STRING,"","",outname);
    cpl_parameterlist_append(plottags, p);

}

/****************************************************************************/
/****************************************************************************/

cpl_error_code sc_plot_single_spec(const cpl_table *spec,
                                   cpl_parameterlist *plottags)
{
/*!
 * \callgraph
 *
 * This program creates a plot from two choosable columns of a CPL table. It
 * expects the following plot tags:
 *
 *      - column name used for the x-axis
 *      - column name used for the y-axis
 *      - plot title
 *      - label for x-axis
 *      - label for y-axis
 *
 * These plot tags have to be set via sc_setplottags_single_spec() routine.
 * Additionally, the gnuplot terminal type, and the directory structure
 * given in the parameter file are added automatically.
 *
 * \b INPUT:
 * \param spec       2 column cpl table containing spectrum
 *                   (wavelength [micron], radiance flux/transmission)
 * \param plottags   CPL parameterlist containing plot tags
 *
 * \b ERRORS:
 * - CPL_ERROR_NONE:           no error occurred
 *
 */

    FILE *specfile;                 /* output ASCII spec file ptr           */
    FILE *gnufile;                  /* gnuplot driverfile                   */
    char filename[SC_MAXLEN]="sci_single_spec.dat";            /* filename  */
    char gnuname1[SC_MAXLEN]="single_gnufile_wxt.gnu";
    char gnuname2[SC_MAXLEN]="single_gnufile_ps.gnu";
    char gnuname3[SC_MAXLEN]="single_gnufile_x11.gnu";
    char ps_filename[SC_MAXLEN];

    char tmpdir[SC_MAXLEN];
    char tmpfilename[SC_MAXLEN];
//     char spectype[SC_MAXLEN];
    char plot_type[SC_MAXLEN];      /* plot type selection (plot_creation)  */

    char system_call[SC_MAXLEN];    /* system call                          */

    int len=0;                      /* table length                         */
    int run=0;                      /* runnning variable                    */
    //int ncol=0;                   /* number of columns in inputspectable  */
    //int dummy=0;                  /* dummy return value for system calls  */
    double lambda=0;                /* wavelength of spectrum               */
    double y_value=0;               /* y-value of plot (RAD/TRA)            */
    int dir_exist_flag=0;           /* Checking for existence of tmp dir    */
    /* Plot limits */
    double plot_xmin=0., plot_xmax=0.;

    /* plot tags */
    char x_label[SC_MAXLEN], y_label[SC_MAXLEN];                  /* labels */
    char x_col[SC_MAXLEN], y_col[SC_MAXLEN];
    char title[SC_MAXLEN];                                        /* title  */

    /* directory + filename for output ps file */
    cpl_parameter *basedirpar, *outdirpar, *filenamepar, *par;

/*--------------------------------------------------------------------------*/
/* ------------------------------ INITIALISING -----------------------------*/
/*--------------------------------------------------------------------------*/

    /* Reading plot tags */
    par=cpl_parameterlist_find(plottags,"plot_title");
    sprintf(title,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"xlabel");
    sprintf(x_label,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"ylabel");
    sprintf(y_label,"%s",cpl_parameter_get_string(par));

    par=cpl_parameterlist_find(plottags,"xcol");
    sprintf(x_col,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"ycol");
    sprintf(y_col,"%s",cpl_parameter_get_string(par));

    par=cpl_parameterlist_find(plottags,"terminal");
    sprintf(plot_type,"%s",cpl_parameter_get_string(par));

    basedirpar = cpl_parameterlist_find(plottags, "inst_dir");
    outdirpar = cpl_parameterlist_find(plottags, "output_dir");
    filenamepar = cpl_parameterlist_find(plottags, "output_name");

    /* Checking table (=input spectra) properties */
    len=cpl_table_get_nrow(spec);   /* length of spectrum                   */
    //ncol=cpl_table_get_ncol(spec);  /* # of columns                       */
    plot_xmin=cpl_table_get_double(spec,x_col,0,NULL);
    plot_xmax=cpl_table_get_double(spec,x_col,len-1,NULL);

    /* Creating tmp-directory */
    sprintf(tmpdir,"__tmpDIRtmp__");
    dir_exist_flag=access(tmpdir,EXIST);
    if (dir_exist_flag == 0){
        cpl_msg_warning(cpl_func,"Directory %s already exists!",tmpdir);
    }
    else
    {
        if (mkdir(tmpdir,0777)) {};
    }

    /* writing .dat file containing spectrum information */
    sprintf(tmpfilename,"%s/%s",tmpdir,filename);
    specfile = fopen(tmpfilename,"w");
    for (run=0;run<len;run++)
    {
        lambda=cpl_table_get_double(spec,x_col,run,NULL);
        y_value=cpl_table_get_double(spec,y_col,run,NULL);
        fprintf(specfile,"%5.6g\t%5.6g\n",lambda,y_value);
    }
    fclose(specfile);

    /* Creating wxt terminal gnuplot driver file */
    if ( (strcmp(plot_type,"W") == 0)  || (strcmp(plot_type,"WX")  == 0) ||
         (strcmp(plot_type,"XW")  == 0) )
    {
        sc_basic_initstring(tmpfilename,SC_MAXLEN);
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname1);
        gnufile = fopen(tmpfilename,"w");
        fprintf(gnufile,"set term wxt\n");
        fprintf(gnufile,"set termoption enhanced  \n");
        fprintf(gnufile,"set nokey\n");
        fprintf(gnufile,"set tmargin 2\n");
        fprintf(gnufile,"set bmargin 5\n");
        fprintf(gnufile,"set lmargin 13\n");
        fprintf(gnufile,"set rmargin 3\n");
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin,plot_xmax);
        fprintf(gnufile,"set title \"%s\"\n",title);
        fprintf(gnufile,"set xlabel \"%s\"\n",x_label);
        fprintf(gnufile,"set ylabel \"%s\"\n",y_label);
        fprintf(gnufile,"plot '%s/%s' using 1:2 with lines\n",tmpdir,
                filename);
        fprintf(gnufile,"\n");
        fclose(gnufile);

        /* Calling gnuplot */
        sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname1);
        if (system(system_call)) {};

        /* Cleaning */
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname1);
        if (remove(tmpfilename)) {};
    }

    /* Creating postscript terminal gnuplot driver file */
    sprintf(ps_filename,"%s/%s/%s_%s_singleplot.ps",
            cpl_parameter_get_string(basedirpar),
            cpl_parameter_get_string(outdirpar),
            cpl_parameter_get_string(filenamepar),title);

    sc_basic_initstring(tmpfilename,SC_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,gnuname2);
    gnufile = fopen(tmpfilename,"w");
    fprintf(gnufile,"set term postscript enhanced color\n");
    fprintf(gnufile,"set output \"%s\"\n",ps_filename);
    fprintf(gnufile,"set nokey\n");
    fprintf(gnufile,"set tmargin 2\n");
    fprintf(gnufile,"set bmargin 5\n");
    fprintf(gnufile,"set lmargin 13\n");
    fprintf(gnufile,"set rmargin 3\n");
    fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin,plot_xmax);
    fprintf(gnufile,"set title \"%s\"\n",title);
    fprintf(gnufile,"set xlabel \"%s\"\n",x_label);
    fprintf(gnufile,"set ylabel \"%s\"\n",y_label);
    fprintf(gnufile,"plot '%s/%s' using 1:2 with lines\n",tmpdir,
            filename);
    fprintf(gnufile,"\n");
    fclose(gnufile);

    /* Calling gnuplot */
    sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname2);
    if (system(system_call)) {};

    /* Cleaning */
    sprintf(tmpfilename,"%s/%s",tmpdir,gnuname2);
    if (remove(tmpfilename)) {};

    /* Creating x11 terminal gnuplot driver file */
    if ( (strcmp(plot_type,"X") == 0)  || (strcmp(plot_type,"WX")  == 0) ||
         (strcmp(plot_type,"XW") == 0) )
    {
        sc_basic_initstring(tmpfilename,SC_MAXLEN);
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname3);
        gnufile = fopen(tmpfilename,"w");
        fprintf(gnufile,"set term x11\n");
        fprintf(gnufile,"set termoption enhanced  \n");
        fprintf(gnufile,"set nokey\n");
        fprintf(gnufile,"set tmargin 2\n");
        fprintf(gnufile,"set bmargin 5\n");
        fprintf(gnufile,"set lmargin 13\n");
        fprintf(gnufile,"set rmargin 3\n");
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin,plot_xmax);
        fprintf(gnufile,"set title \"%s\"\n",title);
        fprintf(gnufile,"set xlabel \"%s\"\n",x_label);
        fprintf(gnufile,"set ylabel \"%s\"\n",y_label);
        fprintf(gnufile,"plot '%s/%s' using 1:2 with lines\n",tmpdir,
                                                              filename);
        fprintf(gnufile,"\n");
        fclose(gnufile);

        /* Calling gnuplot */
        sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname3);
        if (system(system_call)) {};

        /* Cleaning */
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname3);
        if (remove(tmpfilename)) {};
    }

    /* Cleaning */
    sc_basic_initstring(tmpfilename,SC_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,filename);
    if (remove(tmpfilename)) {};
    if (rmdir(tmpdir)) {};

    return CPL_ERROR_NONE;
}


/*****************************************************************************
 *                                                                           *
 *  Routine to plot two spectra                                              *
 *                                                                           *
 ****************************************************************************/

void sc_setplottags_double_spec(cpl_parameterlist *plottags,
                                const char *x_column1,
                                const char *y_column1,
                                const char *title1,
                                const char *x_label1,
                                const char *y_label1,
                                const char *x_column2,
                                const char *y_column2,
                                const char *title2,
                                const char *x_label2,
                                const char *y_label2,
                                const cpl_parameterlist *parlist)
{
/*!
 * \callgraph
 *
 * Routine to set plot tags for double spectrum plot. The two spectra are
 * plotted on two different panels, lying upon another. Hence, tags for the
 * individual data sets have to be given. Spectrum 1 is plotted in the upper
 * Spectrum 2 in the lower panel.
 *
 * The following input tags are required:
 *
 * \b INPUT:
 * \param plottags   Parameter list to be filled with tags
 * \param x_column1  CPL table column name used for the abscissa (Spectrum 1)
 * \param y_column1  CPL table column name used for the ordinate (Spectrum 1)
 * \param title1     Plot title (Spectrum 1)
 * \param x_label1   Label string for x-axis (Spectrum 1)
 * \param y_label1   Label string for y-axis (Spectrum 1)
 * \param x_column2  CPL table column name used for the abscissa (Spectrum 2)
 * \param y_column2  CPL table column name used for the ordinate (Spectrum 2)
 * \param title2     Plot title (Spectrum 2)
 * \param x_label2   Label string for x-axis (Spectrum 2)
 * \param y_label2   Label string for y-axis (Spectrum 2)
 * \param parlist    Parameter list containing information from the parameter
 *                   file
 * \b OUTPUT:
 * \param plottags  Filled parameter list containing tags
 *
 * Base directory, output directory, output name, and gnuplot terminal type
 * are added automatically from parameter list.
 * NOTE: NULL is not accepted, use " " instead!
 *
 */

    cpl_parameter *p;
    const cpl_parameter *cp;
    char terminal[SC_MAXLEN], basedir[SC_MAXLEN];
    char outdir[SC_MAXLEN], outname[SC_MAXLEN];

    p = cpl_parameter_new_value("title1",CPL_TYPE_STRING,"","",title1);
    cpl_parameterlist_append(plottags, p);
    p = cpl_parameter_new_value("xlabel1",CPL_TYPE_STRING,"","",x_label1);
    cpl_parameterlist_append(plottags, p);
    p = cpl_parameter_new_value("ylabel1",CPL_TYPE_STRING,"","",y_label1);
    cpl_parameterlist_append(plottags, p);
    p = cpl_parameter_new_value("xcol1",CPL_TYPE_STRING,"","",x_column1);
    cpl_parameterlist_append(plottags, p);
    p = cpl_parameter_new_value("ycol1",CPL_TYPE_STRING,"","",y_column1);
    cpl_parameterlist_append(plottags, p);

    p = cpl_parameter_new_value("title2",CPL_TYPE_STRING,"","",title2);
    cpl_parameterlist_append(plottags, p);
    p = cpl_parameter_new_value("xlabel2",CPL_TYPE_STRING,"","",x_label2);
    cpl_parameterlist_append(plottags, p);
    p = cpl_parameter_new_value("ylabel2",CPL_TYPE_STRING,"","",y_label2);
    cpl_parameterlist_append(plottags, p);
    p = cpl_parameter_new_value("xcol2",CPL_TYPE_STRING,"","",x_column2);
    cpl_parameterlist_append(plottags, p);
    p = cpl_parameter_new_value("ycol2",CPL_TYPE_STRING,"","",y_column2);
    cpl_parameterlist_append(plottags, p);

    cp = cpl_parameterlist_find_const(parlist,"plot_type");
    sprintf(terminal,"%s",cpl_parameter_get_string(cp));
    p = cpl_parameter_new_value("terminal",CPL_TYPE_STRING,"","",terminal);
    cpl_parameterlist_append(plottags, p);

    cp = cpl_parameterlist_find_const(parlist, "inst_dir");
    sprintf(basedir,"%s",cpl_parameter_get_string(cp));
    p = cpl_parameter_new_value("inst_dir",CPL_TYPE_STRING,"","",basedir);
    cpl_parameterlist_append(plottags, p);

    cp = cpl_parameterlist_find_const(parlist, "output_dir");
    sprintf(outdir,"%s",cpl_parameter_get_string(cp));
    p = cpl_parameter_new_value("output_dir",CPL_TYPE_STRING,"","",outdir);
    cpl_parameterlist_append(plottags, p);

    cp = cpl_parameterlist_find_const(parlist, "output_name");
    sprintf(outname,"%s",cpl_parameter_get_string(cp));
    p = cpl_parameter_new_value("output_name",CPL_TYPE_STRING,"","",outname);
    cpl_parameterlist_append(plottags, p);

}

/****************************************************************************/
/****************************************************************************/

cpl_error_code sc_plot_double_spec(const cpl_table *spec1,
                                   const cpl_table *spec2,
                                   cpl_parameterlist *plottags)
{
/*!
 * \callgraph
 *
 * This program plots two spectra, e.g. object and sky spectra. Input
 * parameters must be CPL tables. The two spectra are plotted on two different
 * panels, lying upon another. Spectrum 1 is plotted in the upper Spectrum 2
 * in the lower panel. Hence, required tags for the individual data sets have
 * to be set using sc_setplottags_double_spec().
 *
 * \b INPUT:
 * \param spec1    Spectrum 1
 * \param spec2    Spectrum 2
 * \param plottags parameter list containing tags
 *
 * \b ERRORS:
 * - CPL_ERROR_NONE:           no error occurred
 */

    //cpl_errorstate err_state;
    cpl_error_code err_code=CPL_ERROR_NONE;

    FILE *specfile;                 /* output ASCII spec file ptr           */
    FILE *gnufile;                  /* gnuplot driverfile                   */
    char filename1[SC_MAXLEN]="plot_spec1.dat";     /* input spec 1         */
    char filename2[SC_MAXLEN]="plot_spec2.dat";     /* input spec 2         */
    char gnuname1[SC_MAXLEN] ="double_gnufile_wxt.gnu";
    char gnuname2[SC_MAXLEN] ="double_gnufile_ps.gnu";
    char gnuname3[SC_MAXLEN] ="double_gnufile_x11.gnu";

    char tmpdir[SC_MAXLEN];
    char tmpfilename[SC_MAXLEN];
    char ps_filename[SC_MAXLEN];    /* Name of postscript file              */
    char plot_type[SC_MAXLEN];      /* plot type selection (plot_creation)  */

    char system_call[SC_MAXLEN];

    int len1=0, len2=0;             /* table lengtes                        */
    //int ncol1=0, ncol2=0;         /* number of columns in inputspectables */
    //int dummy=0;                  /* dummy return value for system calls  */
    int run=0;                      /* runnning variable                    */
    int gridcheckflag=0;            /* flag for checking lambda grid        */
    int dir_exist_flag=0;           /* Checking for existence of tmp dir    */

    /* Plot limits */
    double plot_xmin1=0., plot_xmax1=0., plot_xmin2=0., plot_xmax2=0.;

    double lambda=0;                /* wavelength of spectrum               */
    double y_value=0;               /* y-value of plot (RAD/TRA)            */

    /* Plot tags */
    char title1[SC_MAXLEN];                                  /* title       */
    char x_label1[SC_MAXLEN],y_label1[SC_MAXLEN];            /* labels      */
    char x_col1[SC_MAXLEN], y_col1[SC_MAXLEN];       /* columns to plot     */
    char title2[SC_MAXLEN];                                  /* title       */
    char x_label2[SC_MAXLEN],y_label2[SC_MAXLEN];            /* labels      */
    char x_col2[SC_MAXLEN], y_col2[SC_MAXLEN];       /* columns to plot     */

    /* directory + filename for output ps file */
    cpl_parameter *basedirpar, *outdirpar, *filenamepar, *par;

/*--------------------------------------------------------------------------*/
/* ------------------------------ INITIALISING -----------------------------*/
/*--------------------------------------------------------------------------*/

    /* Reading plot tags */
    par=cpl_parameterlist_find(plottags,"title1");
    sprintf(title1,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"xlabel1");
    sprintf(x_label1,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"ylabel1");
    sprintf(y_label1,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"xcol1");
    sprintf(x_col1,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"ycol1");
    sprintf(y_col1,"%s",cpl_parameter_get_string(par));

    par=cpl_parameterlist_find(plottags,"title2");
    sprintf(title2,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"xlabel2");
    sprintf(x_label2,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"ylabel2");
    sprintf(y_label2,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"xcol2");
    sprintf(x_col2,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"ycol2");
    sprintf(y_col2,"%s",cpl_parameter_get_string(par));

    par=cpl_parameterlist_find(plottags,"terminal");
    sprintf(plot_type,"%s",cpl_parameter_get_string(par));
    basedirpar = cpl_parameterlist_find(plottags, "inst_dir");
    outdirpar = cpl_parameterlist_find(plottags, "output_dir");
    filenamepar = cpl_parameterlist_find(plottags, "output_name");

    /* checking error state */
    //err_state = cpl_errorstate_get();

    /* CHECKING INPUT ------------------------------------------------------*/
    /* Input spectra will be checked against
         - same size (#cols, #rows)
         - same CPL table structure
         - same wavelength grid
         - type of input specta
    */
    /* Checking table (=input spectra) properties */
    len1=cpl_table_get_nrow(spec1);   /* length of spectrum 1               */
    //ncol1=cpl_table_get_ncol(spec1);  /* # of columns spec 1              */

    len2=cpl_table_get_nrow(spec2);   /* length of spectrum 2               */
    //ncol2=cpl_table_get_ncol(spec2);  /* # of columns spec 2              */

    plot_xmin1=cpl_table_get_double(spec1,x_col1,0,NULL);
    plot_xmax1=cpl_table_get_double(spec1,x_col1,len1-1,NULL);
    plot_xmin2=cpl_table_get_double(spec2,x_col2,0,NULL);
    plot_xmax2=cpl_table_get_double(spec2,x_col2,len2-1,NULL);

    /* Checking length of input CPL tables */
    if (len1 != len2)
    {
        cpl_msg_warning(cpl_func,"Input spectra do not have the same size.");
    }

    /* Checking structure of input CPL tables */
    if ( cpl_table_compare_structure(spec1,spec2) )
    {
        cpl_msg_warning(cpl_func,"Input spectra do not have the same "
                                 "structure.");
    }

    /* Checking if wavelength grid is identical */
    run=0;
    for (run=0;run<len1;run++)
    {
        if (cpl_table_get_double(spec1,x_col1,run,NULL) !=
            cpl_table_get_double(spec2,x_col2,run,NULL) )
        {
            gridcheckflag=1;
        }
    }
    if (gridcheckflag)
    {
        cpl_msg_warning(cpl_func,"Input spectra have differing wavelength "
                                 "grids.");
    }

    /* Checking / Creating tmp-directory */
    sprintf(tmpdir,"__tmpDIRtmp__");
    dir_exist_flag=access(tmpdir,EXIST);
    if (dir_exist_flag == 0){
        cpl_msg_warning(cpl_func,"Directory %s already exists!",tmpdir);
    }
    else
    {
        if (mkdir(tmpdir,0777)) {};
    }

    /* writing .dat files containing spectrum information */
    sc_basic_initstring(tmpfilename,SC_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,filename1);
    specfile = fopen(tmpfilename,"w"); /* Science spectrum */
    for (run=0;run<len1;run++)
    {
        lambda=cpl_table_get_double(spec1,x_col1,run,NULL);
        y_value=cpl_table_get_double(spec1,y_col1,run,NULL);
        fprintf(specfile,"%5.6g\t%5.6g\n",lambda,y_value);
    }
    fclose(specfile);

    sc_basic_initstring(tmpfilename,SC_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,filename2);
    specfile = fopen(tmpfilename,"w"); /* Fit spectrum */
    for (run=0;run<len2;run++)
    {
        lambda=cpl_table_get_double(spec2,x_col2,run,NULL);
        y_value=cpl_table_get_double(spec2,y_col2,run,NULL);
        fprintf(specfile,"%5.6g\t%5.6g\n",lambda,y_value);
    }
    fclose(specfile);

/****************************************************************************/

    /* Creating wxt terminal gnuplot file */
    if ( (strcmp(plot_type,"W") == 0)  || (strcmp(plot_type,"WX")  == 0) ||
         (strcmp(plot_type,"XW")  == 0) )
    {
        sc_basic_initstring(tmpfilename,SC_MAXLEN);
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname1);
        gnufile = fopen(tmpfilename,"w");
        fprintf(gnufile,"set term wxt\n");
        fprintf(gnufile,"set termoption enhanced  \n");
        fprintf(gnufile,"# Plotting\n");
        fprintf(gnufile,"set nokey\n");
        fprintf(gnufile,"set tmargin 2\n");
        fprintf(gnufile,"set bmargin 5\n");
        fprintf(gnufile,"set lmargin 12\n");
        fprintf(gnufile,"set rmargin 3\n");
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin1,plot_xmax1);
        fprintf(gnufile,"set termoption font \"Times,7\"\n");
        fprintf(gnufile,"set multiplot layout 2,1 \n");
        fprintf(gnufile,"set title \"%s\"\n",title1);
        fprintf(gnufile,"set xlabel \"%s\"\n",x_label1);
        fprintf(gnufile,"set ylabel \"%s\" offset 1,0\n",y_label1);
        fprintf(gnufile,"set style data boxes\n");
        fprintf(gnufile,"plot '%s/%s' using 1:2 title "
                        " \"%s\" with lines \n",tmpdir,filename1,title1);
        fprintf(gnufile,"set nokey\n");
        fprintf(gnufile,"set ylabel \"%s\" offset 1,0\n",y_label2);
        fprintf(gnufile,"set xlabel \"%s\"\n",x_label2);
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin2,plot_xmax2);

        fprintf(gnufile,"set title \"%s\"\n",title2);
        fprintf(gnufile,"plot '%s/%s' "
                        "using 1:2 title \"%s\" with lines \n",
                        tmpdir,filename2,title2);
        fprintf(gnufile,"unset multiplot    \n");
        fclose(gnufile);

        /* Calling gnuplot */
        sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname1);
        if (system(system_call)) {};

        /* Cleaning */
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname1);
        if (remove(tmpfilename)) {};
    }

    /* Creating postscript terminal gnuplot driver file */
    sprintf(ps_filename,"%s/%s/%s_doubleplot.ps",
            cpl_parameter_get_string(basedirpar),
            cpl_parameter_get_string(outdirpar),
            cpl_parameter_get_string(filenamepar));

    sc_basic_initstring(tmpfilename,SC_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,gnuname2);
    gnufile = fopen(tmpfilename,"w");

    fprintf(gnufile,"set term postscript enhanced color\n");
    fprintf(gnufile,"set output \"%s\"\n", ps_filename);
    fprintf(gnufile,"# Plotting\n");
    fprintf(gnufile,"set nokey\n");
    fprintf(gnufile,"set tmargin 2\n");
    fprintf(gnufile,"set bmargin 5\n");
    fprintf(gnufile,"set lmargin 12\n");
    fprintf(gnufile,"set rmargin 3\n");
    fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin1,plot_xmax1);
    fprintf(gnufile,"set termoption font \"Times,12\"\n");
    fprintf(gnufile,"set multiplot layout 2,1 \n");
    fprintf(gnufile,"set title \"%s\"\n",title1);
    fprintf(gnufile,"set xlabel \"%s\"\n",x_label1);
    fprintf(gnufile,"set ylabel \"%s\" offset 1,0\n",y_label1);
    fprintf(gnufile,"set style data boxes\n");
    fprintf(gnufile,"plot '%s/%s' using 1:2 title "
                    " \"%s\" with lines \n",tmpdir,filename1,title1);
    fprintf(gnufile,"set nokey\n");
    fprintf(gnufile,"set ylabel \"%s\" offset 1,0\n",y_label2);
    fprintf(gnufile,"set xlabel \"%s\"\n",x_label2);
    fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin2,plot_xmax2);

    fprintf(gnufile,"set title \"%s\"\n",title2);
    fprintf(gnufile,"plot '%s/%s' "
                    "using 1:2 title \"%s\" with lines \n",
                    tmpdir,filename2,title2);
    fprintf(gnufile,"unset multiplot    \n");
    fclose(gnufile);

    /* Calling gnuplot */
    sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname2);
    if (system(system_call)) {};

    /* Cleaning */
    sprintf(tmpfilename,"%s/%s",tmpdir,gnuname2);
    if (remove(tmpfilename)) {};

    /* Creating x11 terminal gnuplot file */
    if ( (strcmp(plot_type,"X") == 0)  || (strcmp(plot_type,"WX")  == 0) ||
         (strcmp(plot_type,"XW")  == 0) )
    {
        sc_basic_initstring(tmpfilename,SC_MAXLEN);
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname3);
        gnufile = fopen(tmpfilename,"w");
        fprintf(gnufile,"set term x11\n");
        fprintf(gnufile,"set termoption enhanced  \n");
        fprintf(gnufile,"# Plotting\n");
        fprintf(gnufile,"set nokey\n");
        fprintf(gnufile,"set tmargin 2\n");
        fprintf(gnufile,"set bmargin 5\n");
        fprintf(gnufile,"set lmargin 12\n");
        fprintf(gnufile,"set rmargin 3\n");
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin1,plot_xmax1);
        fprintf(gnufile,"set termoption font \"Times,7\"\n");
        fprintf(gnufile,"set multiplot layout 2,1 \n");
        fprintf(gnufile,"set title \"%s\"\n",title1);
        fprintf(gnufile,"set xlabel \"%s\"\n",x_label1);
        fprintf(gnufile,"set ylabel \"%s\" offset 1,0\n",y_label1);
        fprintf(gnufile,"set style data boxes\n");
        fprintf(gnufile,"plot '%s/%s' using 1:2 title "
                        " \"%s\" with lines \n",tmpdir,filename1,title1);
        fprintf(gnufile,"set nokey\n");
        fprintf(gnufile,"set ylabel \"%s\" offset 1,0\n",y_label2);
        fprintf(gnufile,"set xlabel \"%s\"\n",x_label2);
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin2,plot_xmax2);
        fprintf(gnufile,"set title \"%s\"\n",title2);
        fprintf(gnufile,"plot '%s/%s' "
                        "using 1:2 title \"%s\" with lines \n",
                        tmpdir,filename2,title2);
        fprintf(gnufile,"unset multiplot    \n");
        fclose(gnufile);

        /* Calling gnuplot */
        sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname3);
        if (system(system_call)) {};

        /* Cleaning */
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname3);
        if (remove(tmpfilename)) {};
    }

    /* Cleaning */
    sc_basic_initstring(tmpfilename,SC_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,filename1);
    if (remove(tmpfilename)) {};
    sprintf(tmpfilename,"%s/%s",tmpdir,filename2);
    if (remove(tmpfilename)) {};
    if (rmdir(tmpdir)) {};

    return err_code;
}


/*****************************************************************************
 *                                                                           *
 *  Routine to overplot two spectra and calculate residual                   *
 *                                                                           *
 ****************************************************************************/

void sc_setplottags_overplot_spec(cpl_parameterlist *plottags,
                                const char *x_column,
                                const char *y_column1,
                                const char *y_column2,
                                const char *specname1,
                                const char *specname2,
                                const char *title,
                                const char *x_label,
                                const char *y_label,
                                const cpl_parameterlist *parlist)
{
/*!
 * \callgraph
 *
 * Routine to set plot tags for spectra overplot routine. The two spectra are
 * plotted in one panel (upper), the residual of the two input spectra
 * (Spectrum #1 - Spectrum #2) is plotted in the lower panel.
 *
 * The following input tags are required:
 *
 * \b INPUT:
 * \param plottags   Parameter list to be filled with tags
 * \param x_column   CPL table column name used for the abscissa
 * \param y_column1  CPL table column name used for the ordinate (Spectrum 1)
 * \param y_column2  CPL table column name used for the ordinate (Spectrum 2)
 * \param specname1  Designator of spectrum 1 (used in legend)
 * \param specname2  Designator of spectrum 2 (used in legend)
 * \param title      Plot title
 * \param x_label    Label string for x-axis
 * \param y_label    Label string for y-axis
 * \param parlist    Parameter list containing information from the parameter
 *                   file
 * \b OUTPUT:
 * \param plottags  Filled parameter list containing tags
 *
 * Base directory, output directory, output name, and gnuplot terminal type
 * are added automatically from parameter list.
 * NOTE: NULL is not accepted, use " " instead!
 *
 */

    cpl_parameter *p;
    const cpl_parameter *cp;
    char terminal[SC_MAXLEN], basedir[SC_MAXLEN];
    char outdir[SC_MAXLEN], outname[SC_MAXLEN];

    p = cpl_parameter_new_value("xcol",CPL_TYPE_STRING,"","",x_column);
    cpl_parameterlist_append(plottags, p);
    p = cpl_parameter_new_value("ycol1",CPL_TYPE_STRING,"","",y_column1);
    cpl_parameterlist_append(plottags, p);
    p = cpl_parameter_new_value("ycol2",CPL_TYPE_STRING,"","",y_column2);
    cpl_parameterlist_append(plottags, p);

    p = cpl_parameter_new_value("specname1",CPL_TYPE_STRING,"","",specname1);
    cpl_parameterlist_append(plottags, p);
    p = cpl_parameter_new_value("specname2",CPL_TYPE_STRING,"","",specname2);
    cpl_parameterlist_append(plottags, p);
    p = cpl_parameter_new_value("title",CPL_TYPE_STRING,"","",title);
    cpl_parameterlist_append(plottags, p);
    p = cpl_parameter_new_value("xlabel",CPL_TYPE_STRING,"","",x_label);
    cpl_parameterlist_append(plottags, p);
    p = cpl_parameter_new_value("ylabel",CPL_TYPE_STRING,"","",y_label);
    cpl_parameterlist_append(plottags, p);

    cp=cpl_parameterlist_find_const(parlist,"plot_type");
    sprintf(terminal,"%s",cpl_parameter_get_string(cp));
    p = cpl_parameter_new_value("terminal",CPL_TYPE_STRING,"","",terminal);
    cpl_parameterlist_append(plottags, p);

    cp=cpl_parameterlist_find_const(parlist, "inst_dir");
    sprintf(basedir,"%s",cpl_parameter_get_string(cp));
    p = cpl_parameter_new_value("inst_dir",CPL_TYPE_STRING,"","",basedir);
    cpl_parameterlist_append(plottags, p);

    cp=cpl_parameterlist_find_const(parlist, "output_dir");
    sprintf(outdir,"%s",cpl_parameter_get_string(cp));
    p = cpl_parameter_new_value("output_dir",CPL_TYPE_STRING,"","",outdir);
    cpl_parameterlist_append(plottags, p);

    cp=cpl_parameterlist_find_const(parlist, "output_name");
    sprintf(outname,"%s",cpl_parameter_get_string(cp));
    p = cpl_parameter_new_value("output_name",CPL_TYPE_STRING,"","",outname);
    cpl_parameterlist_append(plottags, p);

}



/****************************************************************************/
/****************************************************************************/

cpl_error_code sc_overplot_spec(const cpl_table *spec,
                                cpl_parameterlist *plottags)
{
/*!
 * \callgraph
 *
 * This program is dedicated to produce plots for a direct comparison of
 * spectra. The input must be a CPL table containing ONE column used as
 * abscissa ('x_column' in plot tags), but TWO individual columns must be
 * present to be used as individual ordinates ('y_column1' and 'y_column2',
 * see plottags). This is necessary, as it must be ensured that both spectra
 * are defined on the same wavelength grid.
 *
 * The comparison is made with two plot panels lying above each other. The
 * upper panel of the plot is a direct overplot of the two spectra, whereas in
 * the lower panel the Residual plot (Spectrum #1 - Spectrum #2) is given.
 *
 * \b INPUT:
 * \param spec        Table containing both spectra (output of molecfit)
 * \param plottags    parameter list containing tags
 *
 * \b ERRORS:
 * - CPL_ERROR_NONE:           no error occurred
 */

    cpl_error_code err_code=CPL_ERROR_NONE;

    FILE *specfile;                 /* output ASCII spec file ptr           */
    FILE *gnufile;                  /* gnuplot driverfile                   */
    char filename1[SC_MAXLEN]="sci_spec.dat";     /* input spec 1           */
    char filename2[SC_MAXLEN]="fit_spec.dat";     /* input spec 2           */
    char filename3[SC_MAXLEN]="resi_spec.dat";    /* residual spec          */
    char gnuname1[SC_MAXLEN]="overplot_gnufile_wxt.gnu";
    char gnuname2[SC_MAXLEN]="overplot_gnufile_ps.gnu";
    char gnuname3[SC_MAXLEN]="overplot_gnufile_x11.gnu";
    char ps_filename[SC_MAXLEN];    /* Name of postscript file              */
    char tmpdir[SC_MAXLEN];
    char tmpfilename[SC_MAXLEN];
    char system_call[SC_MAXLEN];    /* string for system call               */
    char plot_type[SC_MAXLEN];      /* plot type selection (plot_creation)  */
    int  dir_exist_flag=0;          /* Checking for existence of tmp dir    */

    /* Plot limits */
    double plot_xmin=0., plot_xmax=0., plot_ymin=0., plot_ymax=0., tmp1=0.,
      tmp2=0., dy=0.;
    char x_label[SC_MAXLEN],y_label[SC_MAXLEN];                   /* labels */
    char title[SC_MAXLEN];                                        /* title  */
    char x_col[SC_MAXLEN], y_col1[SC_MAXLEN], y_col2[SC_MAXLEN];
    char specname1[SC_MAXLEN], specname2[SC_MAXLEN];      /*name of spectra */

    int len1=0;                     /* table length                         */
    int run=0;                      /* runnning variable                    */
    double y_value=0,               /* y-value of plot                      */
           ymin=0;                  /* minimum residual                     */

    /* directory + filename for output ps file */
    cpl_parameter *basedirpar, *outdirpar, *filenamepar, *par;

/*--------------------------------------------------------------------------*/
/* ------------------------------ INITIALISING -----------------------------*/
/*--------------------------------------------------------------------------*/

    /* Reading plot tags */
    par=cpl_parameterlist_find(plottags,"xcol");
    sprintf(x_col,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"ycol1");
    sprintf(y_col1,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"ycol2");
    sprintf(y_col2,"%s",cpl_parameter_get_string(par));

    par=cpl_parameterlist_find(plottags,"specname1");
    sprintf(specname1,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"specname2");
    sprintf(specname2,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"title");
    sprintf(title,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"xlabel");
    sprintf(x_label,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"ylabel");
    sprintf(y_label,"%s",cpl_parameter_get_string(par));

    par=cpl_parameterlist_find(plottags,"terminal");
    sprintf(plot_type,"%s",cpl_parameter_get_string(par));
    basedirpar = cpl_parameterlist_find(plottags, "inst_dir");
    outdirpar = cpl_parameterlist_find(plottags, "output_dir");
    filenamepar = cpl_parameterlist_find(plottags, "output_name");

    /* CHECKING INPUT ------------------------------------------------------*/
    /* Checking table (=input spectra) properties                           */
    len1=cpl_table_get_nrow(spec);   /* length of spectrum                  */

    /* Checking / Creating tmp-directory */
    sprintf(tmpdir,"__tmpDIRtmp__");
    dir_exist_flag=access(tmpdir,EXIST);
    if (dir_exist_flag == 0){
        cpl_msg_warning(cpl_func,"Directory %s already exists!",tmpdir);
    }
    else
    {
        if (mkdir(tmpdir,0777)) {};
    }

    /* writing .dat files containing spectrum information                   */
    /* Science spectrum */
    sc_basic_initstring(tmpfilename,SC_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,filename1);
    specfile = fopen(tmpfilename,"w");
    for (run=0;run<len1;run++)
    {
        fprintf(specfile,"%5.6g\t%5.6g\n",
                cpl_table_get_double(spec,x_col,run,NULL),
                cpl_table_get_double(spec,y_col1,run,NULL));
    }
    fclose(specfile);

    /* Fit spectrum */
    sc_basic_initstring(tmpfilename,SC_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,filename2);
    specfile = fopen(tmpfilename,"w");
    for (run=0;run<len1;run++)
    {
        fprintf(specfile,"%5.6g\t%5.6g\n",
                cpl_table_get_double(spec,x_col,run,NULL),
                cpl_table_get_double(spec,y_col2,run,NULL));
    }
    fclose(specfile);

    /* Residual spectrum: Sci - Fit */
    sc_basic_initstring(tmpfilename,SC_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,filename3);
    specfile = fopen(tmpfilename,"w");
    for (run=0;run<len1;run++)
    {
        /* Consider pixels with non-zero weight only */
        if (cpl_table_get_double(spec,"weight",run,NULL) > 0) {
            y_value=cpl_table_get_double(spec,y_col1,run,NULL)-
                    cpl_table_get_double(spec,y_col2,run,NULL);
        } else {
            y_value=0;
        }
        fprintf(specfile,"%5.6g\t%5.6g\n",
                cpl_table_get_double(spec,x_col,run,NULL),y_value);
        /* Get minimum for plot limits */
        if (y_value < ymin) {
            ymin = y_value;
        }
    }
    fclose(specfile);

    /* Setting plot limits */
    plot_xmin=cpl_table_get_double(spec,x_col,0,NULL);
    plot_xmax=cpl_table_get_double(spec,x_col,len1-1,NULL);
    tmp1=cpl_table_get_column_min(spec,y_col1);
    tmp2=cpl_table_get_column_min(spec,y_col2);
    if (tmp1 <= tmp2)
    {
        plot_ymin = tmp1;
    }
    else
    {
        plot_ymin = tmp2;
    }
    if (ymin < plot_ymin) {
        plot_ymin = ymin;
    }
    tmp1=cpl_table_get_column_max(spec,y_col1);
    tmp2=cpl_table_get_column_max(spec,y_col2);
    if (tmp1 >= tmp2)
    {
        plot_ymax = tmp1;
    }
    else
    {
        plot_ymax = tmp2;
    }
    /* Adding + 5% for better looking plot  */
    dy = plot_ymax - plot_ymin;
    plot_ymin = plot_ymin - 0.05 *dy;
    plot_ymax = plot_ymax + 0.05 *dy;


/****************************************************************************/

    /* Creating wxt terminal gnuplot file */
    if ( (strcmp(plot_type,"W") == 0)  || (strcmp(plot_type,"WX")  == 0) ||
         (strcmp(plot_type,"XW")  == 0) )
    {
        sc_basic_initstring(tmpfilename,SC_MAXLEN);
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname1);
        gnufile = fopen(tmpfilename,"w");

        fprintf(gnufile,"set term wxt\n");
        fprintf(gnufile,"set termoption enhanced  \n");
        fprintf(gnufile,"# Plotting\n");
        fprintf(gnufile,"set key at screen 0.76, 0.55 autotitle column box "
                        "samplen 1 left\n");
        fprintf(gnufile,"set tmargin 0\n");
        fprintf(gnufile,"set bmargin 5\n");
        fprintf(gnufile,"set lmargin 12\n");
        fprintf(gnufile,"set rmargin 14\n");
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin,plot_xmax);
        fprintf(gnufile,"set yrange [ %g : %g ]\n",plot_ymin,plot_ymax);
        fprintf(gnufile,"unset title\n");
        fprintf(gnufile,"set termoption font \"Times,7\"\n");
        fprintf(gnufile,"set multiplot layout 2,1 title \"%s\" \n",title);
        fprintf(gnufile,"set xlabel \"%s\"\n",x_label);
        fprintf(gnufile,"set ylabel \"%s\" offset 1,0\n",y_label);
        fprintf(gnufile,"set style data boxes\n");
        fprintf(gnufile,"plot '%s/%s' using 1:2 title \"%s\" "
                "with lines lt -1,  '%s/%s' "
                "using 1:2 title \"%s\" with lines lt 8\n",
                tmpdir,filename1,specname1,tmpdir,filename2,specname2);
        fprintf(gnufile,"set key at screen 0.73, 0.1 autotitle column box "
                "samplen 1 left\n");
        fprintf(gnufile,"set ylabel \"Residual (%s-%s)\" offset 1,0\n",
                specname1,specname2);
        fprintf(gnufile,"set y2tics border out scale 1,0.5 nomirror norotate"
                "  offset character 0, 0, 0 autofreq \n");
        fprintf(gnufile,"set ytics nomirror\n");
        fprintf(gnufile,"set y2tics nomirror textcolor lt 2\n");
        fprintf(gnufile,"set ylabel \"Residual (input-best-fit sky)\" offset"
                " 1,0\n");
        fprintf(gnufile,"set y2label \"Residual (input-best-fit sky)\" "
                "textcolor lt 2\n");

        fprintf(gnufile,"set xlabel \"%s\"\n",x_label);
        fprintf(gnufile,"plot '__tmpDIRtmp__/resi_spec.dat' using 1:2 title"
                " \"original scaling\" with lines lw 2 lt -1 "
                "axes x1y1,'__tmpDIRtmp__/resi_spec.dat' using 1:2 title "
                "\"optimal scaling\" with lines lt 2 axes x2y2\n");
/*        fprintf(gnufile,"plot '%s/%s' using 1:2 "
                        "with lines lt -1\n",tmpdir,filename3);*/
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin,plot_xmax);
        fprintf(gnufile,"set yrange [ %g : %g ]\n",plot_ymin,plot_ymax);
        fprintf(gnufile,"set autoscale y2\n");
        fprintf(gnufile,"set title \n");
        fprintf(gnufile,"unset multiplot    \n");
        fclose(gnufile);

        /* Calling gnuplot */
        sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname1);
        if (system(system_call)) {};

        /* Cleaning */
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname1);
        if (remove(tmpfilename)) {};
    }

    /* Creating postscript terminal gnuplot driver file */

    if ( strncmp(cpl_parameter_get_string(outdirpar), "/", 1) == 0 )
    {
        sprintf(ps_filename,"%s/%s_overplot.ps",
                        cpl_parameter_get_string(outdirpar),
                        cpl_parameter_get_string(filenamepar));
    }
    else
    {
        sprintf(ps_filename,"%s/%s/%s_overplot.ps",
                        cpl_parameter_get_string(basedirpar),
                        cpl_parameter_get_string(outdirpar),
                        cpl_parameter_get_string(filenamepar));
    }

    sc_basic_initstring(tmpfilename,SC_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,gnuname2);
    gnufile = fopen(tmpfilename,"w");

    fprintf(gnufile,"set term postscript enhanced color\n");
    fprintf(gnufile,"set output \"%s\"\n",ps_filename);
    fprintf(gnufile,"# Plotting\n");
    fprintf(gnufile,"set key at screen 0.70, 0.55 autotitle column box "
                    "samplen 1 left\n");
    fprintf(gnufile,"set tmargin 0\n");
    fprintf(gnufile,"set bmargin 5\n");
    fprintf(gnufile,"set lmargin 12\n");
    fprintf(gnufile,"set rmargin 14\n");
    fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin,plot_xmax);
    fprintf(gnufile,"set yrange [ %g : %g ]\n",plot_ymin,plot_ymax);
    fprintf(gnufile,"unset title\n");
    fprintf(gnufile,"set termoption font \"Times,12\"\n");
    fprintf(gnufile,"set multiplot layout 2,1 title \"%s\" \n",title);
    fprintf(gnufile,"set xlabel \"%s\"\n",x_label);
    fprintf(gnufile,"set ylabel \"%s\" offset 1,0\n",y_label);
    fprintf(gnufile,"set style data boxes\n");
    fprintf(gnufile,"plot '%s/%s' using 1:2 title \"%s\" "
            "with lines lt -1,  '%s/%s' "
            "using 1:2 title \"%s\" with lines lt 8\n",
            tmpdir,filename1,specname1,tmpdir,filename2,specname2);

    fprintf(gnufile,"set key at screen 0.66, 0.07 autotitle column box "
            "samplen 1 left\n");
    fprintf(gnufile,"set ylabel \"Residual (%s-%s)\" offset 1,0\n",
            specname1,specname2);
    fprintf(gnufile,"set y2tics border out scale 1,0.5 nomirror norotate"
            "  offset character 0, 0, 0 autofreq \n");
    fprintf(gnufile,"set style line 1 lt 2 lc rgb \"red\" lw 3\n");
    fprintf(gnufile,"set style line 2 lt 1 lc rgb \"green\" lw 1\n");
    fprintf(gnufile,"set ytics nomirror\n");
    fprintf(gnufile,"set y2tics nomirror textcolor lt 2\n");
    fprintf(gnufile,"set ylabel \"Residual (input-best-fit sky)\" offset"
            " 1,0\n");
    fprintf(gnufile,"set y2label \"Residual (input-best-fit sky)\" "
            "textcolor lt 2\n");

    fprintf(gnufile,"set xlabel \"%s\"\n",x_label);
    fprintf(gnufile,"plot '__tmpDIRtmp__/resi_spec.dat' using 1:2 title"
            " \"original scaling\" with lines lw 2 lt -1 "
            "axes x1y1,'__tmpDIRtmp__/resi_spec.dat' using 1:2 title "
            "\"optimal scaling\" with lines ls 2 axes x2y2\n");
/*        fprintf(gnufile,"plot '%s/%s' using 1:2 "
                    "with lines lt -1\n",tmpdir,filename3);*/
    fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin,plot_xmax);
    fprintf(gnufile,"set yrange [ %g : %g ]\n",plot_ymin,plot_ymax);
    fprintf(gnufile,"set autoscale y2\n");
    fprintf(gnufile,"set title \n");
    fprintf(gnufile,"unset multiplot    \n");
    fclose(gnufile);
/*
    fprintf(gnufile,"# Plotting\n");
    fprintf(gnufile,"set key at screen 0.77, 0.55 autotitle column box "
                    "samplen 1 left\n");
    fprintf(gnufile,"set tmargin 0\n");
    fprintf(gnufile,"set bmargin 5\n");
    fprintf(gnufile,"set lmargin 12\n");
    fprintf(gnufile,"set rmargin 3\n");
    fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin,plot_xmax);
    fprintf(gnufile,"set yrange [ %g : %g ]\n",plot_ymin,plot_ymax);
    fprintf(gnufile,"unset title\n");
    fprintf(gnufile,"set termoption font \"Times,12\"\n");
    fprintf(gnufile,"set multiplot layout 2,1 title \"%s\" \n",title);
    fprintf(gnufile,"set xlabel \"%s\"\n",x_label);
    fprintf(gnufile,"set ylabel \"%s\" offset 1,0\n",y_label);
    fprintf(gnufile,"set style data boxes\n");
    fprintf(gnufile,"plot '%s/%s' using 1:2 title \"%s\" "
            "with lines lt -1,  '%s/%s' "
            "using 1:2 title \"%s\" with lines lt 8\n",
            tmpdir,filename1,specname1,tmpdir,filename2,specname2);
    fprintf(gnufile,"set nokey\n");
    fprintf(gnufile,"set ylabel \"Residual (%s-%s)\" offset 1,0\n",
            specname1,specname2);
    fprintf(gnufile,"set xlabel \"%s\"\n",x_label);
    fprintf(gnufile,"plot '%s/%s' using 1:2 "
                    "with lines lt -1\n",tmpdir,filename3);
    fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin,plot_xmax);
    fprintf(gnufile,"set yrange [ %g : %g ]\n",plot_ymin,plot_ymax);
    fprintf(gnufile,"set title \n");
    fprintf(gnufile,"unset multiplot    \n");
    fclose(gnufile);*/

    /* Calling gnuplot */
    sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname2);
    if (system(system_call)) {};

    /* Cleaning */
    sprintf(tmpfilename,"%s/%s",tmpdir,gnuname2);
    if (remove(tmpfilename)) {};

    /* Creating x11 terminal gnuplot file */
    if ( (strcmp(plot_type,"X") == 0)  || (strcmp(plot_type,"WX")  == 0) ||
         (strcmp(plot_type,"XW")  == 0) )
    {
        sc_basic_initstring(tmpfilename,SC_MAXLEN);
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname3);
        gnufile = fopen(tmpfilename,"w");

        fprintf(gnufile,"set term x11\n");
        fprintf(gnufile,"set termoption enhanced  \n");
        fprintf(gnufile,"# Plotting\n");
        fprintf(gnufile,"set key at screen 0.72, 0.55 autotitle column box "
                        "samplen 1 left\n");
        fprintf(gnufile,"set tmargin 0\n");
        fprintf(gnufile,"set bmargin 5\n");
        fprintf(gnufile,"set lmargin 12\n");
        fprintf(gnufile,"set rmargin 14\n");
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin,plot_xmax);
        fprintf(gnufile,"set yrange [ %g : %g ]\n",plot_ymin,plot_ymax);
        fprintf(gnufile,"unset title\n");
        fprintf(gnufile,"set termoption font \"Times,15\"\n");
        fprintf(gnufile,"set multiplot layout 2,1 title \"%s\" \n",title);
        fprintf(gnufile,"set xlabel \"%s\"\n",x_label);
        fprintf(gnufile,"set ylabel \"%s\" offset 1,0\n",y_label);
        fprintf(gnufile,"set style data boxes\n");
        fprintf(gnufile,"plot '%s/%s' using 1:2 title \"%s\" "
                "with lines lt -1,  '%s/%s' "
                "using 1:2 title \"%s\" with lines lt 8\n",
                tmpdir,filename1,specname1,tmpdir,filename2,specname2);
        fprintf(gnufile,"set key at screen 0.68, 0.1 autotitle column box "
                "samplen 1 left\n");
        fprintf(gnufile,"set ylabel \"Residual (%s-%s)\" offset 1,0\n",
                specname1,specname2);
        fprintf(gnufile,"set y2tics border out scale 1,0.5 nomirror norotate"
                "  offset character 0, 0, 0 autofreq \n");
        fprintf(gnufile,"set ytics nomirror\n");
        fprintf(gnufile,"set y2tics nomirror textcolor lt 3\n");
        fprintf(gnufile,"set ylabel \"Residual (input-best-fit sky)\" offset"
                " 1,0\n");
        fprintf(gnufile,"set y2label \"Residual (input-best-fit sky)\" "
                "textcolor lt 3\n");

        fprintf(gnufile,"set xlabel \"%s\"\n",x_label);
        fprintf(gnufile,"plot '__tmpDIRtmp__/resi_spec.dat' using 1:2 title"
                " \"original scaling\" with lines lw 2 lt -1 "
                "axes x1y1,'__tmpDIRtmp__/resi_spec.dat' using 1:2 title "
                "\"optimal scaling\" with lines lt 3 axes x2y2\n");
/*        fprintf(gnufile,"plot '%s/%s' using 1:2 "
                        "with lines lt -1\n",tmpdir,filename3);*/
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin,plot_xmax);
        fprintf(gnufile,"set yrange [ %g : %g ]\n",plot_ymin,plot_ymax);
        fprintf(gnufile,"set autoscale y2\n");
        fprintf(gnufile,"set title \n");
        fprintf(gnufile,"unset multiplot    \n");
        fclose(gnufile);
//         fprintf(gnufile,"set termoption enhanced  \n");
//         fprintf(gnufile,"# Plotting\n");
//         fprintf(gnufile,"set key at screen 0.77, 0.55 autotitle column box "
//                         "samplen 1 left\n");
//         fprintf(gnufile,"set tmargin 0\n");
//         fprintf(gnufile,"set bmargin 5\n");
//         fprintf(gnufile,"set lmargin 12\n");
//         fprintf(gnufile,"set rmargin 3\n");
//         fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin,plot_xmax);
//         fprintf(gnufile,"set yrange [ %g : %g ]\n",plot_ymin,plot_ymax);
//         fprintf(gnufile,"unset title\n");
//         fprintf(gnufile,"set termoption font \"Times,7\"\n");
//         fprintf(gnufile,"set multiplot layout 2,1 title \"%s\" \n",title);
//         fprintf(gnufile,"set xlabel \"%s\"\n",x_label);
//         fprintf(gnufile,"set ylabel \"%s\" offset 1,0\n",y_label);
//         fprintf(gnufile,"set style data boxes\n");
//         fprintf(gnufile,"plot '%s/%s' using 1:2 title \"%s\" "
//                 "with lines lt -1,  '%s/%s' "
//                 "using 1:2 title \"%s\" with lines lt 8\n",
//                 tmpdir,filename1,specname1,tmpdir,filename2,specname2);
//         fprintf(gnufile,"set nokey\n");
//         fprintf(gnufile,"set ylabel \"Residual (%s-%s)\" offset 1,0\n",
//                 specname1,specname2);
//         fprintf(gnufile,"set xlabel \"%s\"\n",x_label);
//         fprintf(gnufile,"plot '%s/%s' using 1:2 "
//                         "with lines lt -1\n",tmpdir,filename3);
//         fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin,plot_xmax);
//         fprintf(gnufile,"set yrange [ %g : %g ]\n",plot_ymin,plot_ymax);
//         fprintf(gnufile,"set title \n");
//         fprintf(gnufile,"unset multiplot    \n");
//         fclose(gnufile);

        /* Calling gnuplot */
        sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname3);
        if (system(system_call)) {};

        /* Cleaning */
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname3);
        if (remove(tmpfilename)) {};
    }

    /* Cleaning */
    sc_basic_initstring(tmpfilename,SC_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,filename1);
    if (remove(tmpfilename)) {};
    sprintf(tmpfilename,"%s/%s",tmpdir,filename2);
    if (remove(tmpfilename)) {};
    sprintf(tmpfilename,"%s/%s",tmpdir,filename3);
    if (remove(tmpfilename)) {};
    if (rmdir(tmpdir)) {};

    return err_code;
}

/*****************************************************************************
 *                                                                           *
 *  Routine to plot a histogram                                              *
 *                                                                           *
 ****************************************************************************/

void sc_setplottags_hist(cpl_parameterlist *plottags, const char *title,
                         const char *x_label, const char *y_label,
                         const cpl_parameterlist *parlist)
{
/*!
 * \callgraph
 *
 * Routine to set plot tags for histogram plots.
 *
 * The following input tags are required:
 *
 * \b INPUT:
 * \param plottags   Parameter list to be filled with tags
 * \param title      Plot title
 * \param x_label    Label string for x-axis
 * \param y_label    Label string for y-axis
 * \param parlist    Parameter list containing information from the parameter
 *                   file
 * \b OUTPUT:
 * \param plottags  Filled parameter list containing tags
 *
 * Base directory, output directory, output name, and gnuplot terminal type
 * are added automatically from parameter list.
 * NOTE: NULL is not accepted, use " " instead!
 *
 */

    cpl_parameter *p;
    const cpl_parameter *cp;
    char terminal[SC_MAXLEN], basedir[SC_MAXLEN];
    char outdir[SC_MAXLEN], outname[SC_MAXLEN];

    p = cpl_parameter_new_value("plot_title",CPL_TYPE_STRING,"","",title);
    cpl_parameterlist_append(plottags, p);
    p = cpl_parameter_new_value("xlabel",CPL_TYPE_STRING,"","",x_label);
    cpl_parameterlist_append(plottags, p);
    p = cpl_parameter_new_value("ylabel",CPL_TYPE_STRING,"","",y_label);
    cpl_parameterlist_append(plottags, p);

    cp = cpl_parameterlist_find_const(parlist,"plot_type");
    sprintf(terminal,"%s",cpl_parameter_get_string(cp));
    p = cpl_parameter_new_value("terminal",CPL_TYPE_STRING,"","",terminal);
    cpl_parameterlist_append(plottags, p);

    cp = cpl_parameterlist_find_const(parlist, "inst_dir");
    sprintf(basedir,"%s",cpl_parameter_get_string(cp));
    p = cpl_parameter_new_value("inst_dir",CPL_TYPE_STRING,"","",basedir);
    cpl_parameterlist_append(plottags, p);

    cp = cpl_parameterlist_find_const(parlist, "output_dir");
    sprintf(outdir,"%s",cpl_parameter_get_string(cp));
    p = cpl_parameter_new_value("output_dir",CPL_TYPE_STRING,"","",outdir);
    cpl_parameterlist_append(plottags, p);

    cp = cpl_parameterlist_find_const(parlist, "output_name");
    sprintf(outname,"%s",cpl_parameter_get_string(cp));
    p = cpl_parameter_new_value("output_name",CPL_TYPE_STRING,"","",outname);
    cpl_parameterlist_append(plottags, p);

}

/****************************************************************************/
/****************************************************************************/

cpl_error_code sc_plot_hist(const cpl_table *histdat,
                            cpl_parameterlist *plottags)
{
/*!
 * \callgraph
 *
 * This program creates a histogram from a 2-column CPL table. This table must
 * have two columns labelled as "bins" and "counts" (see create_hist() ).
 *
 * \b INPUT:
 * \param histdat       2 column cpl table
 * \param plottags      parameter list containing info from parameter file
 *
 * \b ERRORS:
 * - CPL_ERROR_NONE: no error occurred
 * - CPL_ERROR_ILLEGAL_INPUT
 */


    //cpl_errorstate err_state;
    cpl_error_code err_code=CPL_ERROR_NONE;

    FILE *specfile;                 /* output ASCII spec file ptr           */
    FILE *gnufile;                  /* gnuplot driverfile                   */
    char filename1[SC_MAXLEN]="hist_data.dat";     /* data for histogram    */
    char gnuname1[SC_MAXLEN]="hist_plot_wxt.gnu";
    char gnuname2[SC_MAXLEN]="hist_plot_ps.gnu";
    char gnuname3[SC_MAXLEN]="hist_plot_x11.gnu";
    char ps_filename[SC_MAXLEN];    /* Name of postscript file              */
    char tmpdir[SC_MAXLEN];
    char tmpfilename[SC_MAXLEN];
    char system_call[SC_MAXLEN];    /* string for system call               */
    char plot_type[SC_MAXLEN];      /* plot type selection (plot_creation)  */

    int  dir_exist_flag=0;          /* Checking for existence of tmp dir    */
    //int  dummy=0;                 /* dummy return value for system calls  */
    int run=0;                      /* runnning variable                    */

    double binlims=0.;
    int cts=0;

    /* Plot tags */
    char plot_title[SC_MAXLEN];
    char x_colname[SC_MAXLEN], y_colname[SC_MAXLEN];
    char err_msg[SC_MAXLEN];        /* error message to be returned         */

    int  len1=0, ncol1=0;           /* Dimension of cpl-table 'histdat'     */

    /* directory + filename for output ps file */
    cpl_parameter *basedirpar, *outdirpar, *filenamepar, *par;

    cpl_array *column_names1; /* Col. names of input tables */


/*--------------------------------------------------------------------------*/
/* ------------------------------ INITIALISING -----------------------------*/
/*--------------------------------------------------------------------------*/
    /* Reading plot tags */
    par=cpl_parameterlist_find(plottags,"plot_title");
    sprintf(plot_title,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"xlabel");
    sprintf(x_colname,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"ylabel");
    sprintf(y_colname,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"terminal");
    sprintf(plot_type,"%s",cpl_parameter_get_string(par));

    /* checking error state */
    //err_state = cpl_errorstate_get();

    /* CHECKING INPUT ------------------------------------------------------*/
    /* Checking table (=input data) properties                              */
    len1=cpl_table_get_nrow(histdat);   /* # of bins                        */
    ncol1=cpl_table_get_ncol(histdat);  /* # of columns spectrum            */
    column_names1 = cpl_table_get_column_names(histdat); /* Col names       */

    /* Checking whether 2 columns for input file exist */
    if (ncol1 != 2)
    {
        sprintf(err_msg,"Number of columns not equal 2 in data file...");
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                            "%s", err_msg);
    }

    /* Checking / Creating tmp-directory */
    sprintf(tmpdir,"__tmpDIRtmp__");
    dir_exist_flag=access(tmpdir,EXIST);
    if (dir_exist_flag == 0){
        cpl_msg_warning(cpl_func,"Directory %s already exists!",tmpdir);
    }
    else
    {
        if (mkdir(tmpdir,0777)) {};
    }

    /* writing .dat files containing spectrum information                   */
    /* Science spectrum */
    sc_basic_initstring(tmpfilename,SC_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,filename1);
    specfile = fopen(tmpfilename,"w");
    for (run=0;run<len1;run++)
    {
        binlims=cpl_table_get_double(histdat,"bins",run,NULL);
        cts=cpl_table_get_int(histdat,"counts",run,NULL);
        fprintf(specfile,"%5.3g\t%i\n",binlims,cts);
    }
    fclose(specfile);

    /* Checking plot options */

    /* Creating wxt terminal gnuplot driver file */
    if ( (strcmp(plot_type,"W") == 0)  || (strcmp(plot_type,"WX")  == 0) ||
         (strcmp(plot_type,"XW")  == 0) )
    {
        sc_basic_initstring(tmpfilename,SC_MAXLEN);
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname1);
        gnufile = fopen(tmpfilename,"w");
        fprintf(gnufile,"set term wxt\n");
        fprintf(gnufile,"set termoption enhanced\n");
        fprintf(gnufile,"set title \"%s\"\n",plot_title);
        fprintf(gnufile,"set termoption font \"Times,9\"\n");
        fprintf(gnufile,"set xlabel '%s'\n",x_colname);
        fprintf(gnufile,"set ylabel '%s'\n",y_colname);
        fprintf(gnufile,"set boxwidth 0.75 absolute\n");
        fprintf(gnufile,"set style fill solid 1.00 border -1\n");
        fprintf(gnufile,"set style histogram rowstacked\n");
        fprintf(gnufile,"set style data histograms\n");
        fprintf(gnufile,"plot '%s/%s' u 1:2 smooth frequency with histeps"
                        " t \"%s\"\n",tmpdir,filename1,y_colname);
        fprintf(gnufile,"unset xlabel\n");
        fprintf(gnufile,"unset ylabel\n");
        fprintf(gnufile,"\n");
        fclose(gnufile);

        /* Calling gnuplot */
        sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname1);
        if (system(system_call)) {};

        /* Cleaning */
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname1);
        if (remove(tmpfilename)) {};
    }

    /* Creating postscript terminal gnuplot driver file */
    basedirpar = cpl_parameterlist_find(plottags, "inst_dir");
    outdirpar = cpl_parameterlist_find(plottags, "output_dir");
    filenamepar = cpl_parameterlist_find(plottags, "output_name");

    sprintf(ps_filename,"%s/%s/%s_histogram.ps",
            cpl_parameter_get_string(basedirpar),
            cpl_parameter_get_string(outdirpar),
            cpl_parameter_get_string(filenamepar));

    sc_basic_initstring(tmpfilename,SC_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,gnuname2);
    gnufile = fopen(tmpfilename,"w");
    fprintf(gnufile,"set term postscript enhanced color\n");
    fprintf(gnufile,"set output \"%s\"\n",ps_filename);
    fprintf(gnufile,"set title \"%s\"\n",plot_title);
    fprintf(gnufile,"set termoption font \"Times,9\"\n");
    fprintf(gnufile,"set xlabel '%s'\n",x_colname);
    fprintf(gnufile,"set ylabel '%s'\n",y_colname);
    fprintf(gnufile,"set boxwidth 0.9 absolute\n");
    fprintf(gnufile,"set style fill solid 1.00 border -1\n");
    fprintf(gnufile,"set style data histograms\n");
    fprintf(gnufile,"plot '%s/%s' u 1:2 smooth frequency with histeps"
                    " t \"%s\"\n",tmpdir,filename1,y_colname);
    fprintf(gnufile,"unset xlabel\n");
    fprintf(gnufile,"unset ylabel\n");
    fprintf(gnufile,"\n");
    fclose(gnufile);

    /* Calling gnuplot */
    sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname2);
    if (system(system_call)) {};

    /* Cleaning */
    sprintf(tmpfilename,"%s/%s",tmpdir,gnuname2);
    if (remove(tmpfilename)) {};

    /* Creating x11 terminal gnuplot driver file */
    if ( (strcmp(plot_type,"X") == 0)  || (strcmp(plot_type,"WX")  == 0) ||
         (strcmp(plot_type,"XW")  == 0) )
    {
        sc_basic_initstring(tmpfilename,SC_MAXLEN);
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname3);
        gnufile = fopen(tmpfilename,"w");
        fprintf(gnufile,"set term x11\n");
        fprintf(gnufile,"set termoption enhanced\n");

        fprintf(gnufile,"set title \"%s\"\n",plot_title);
        fprintf(gnufile,"set termoption font \"Times,9\"\n");
        fprintf(gnufile,"set xlabel '%s'\n",x_colname);
        fprintf(gnufile,"set ylabel '%s'\n",y_colname);
        fprintf(gnufile,"set boxwidth 0.75 absolute\n");
        fprintf(gnufile,"set style fill solid 1.00 border -1\n");
        fprintf(gnufile,"set style histogram rowstacked\n");
        fprintf(gnufile,"set style data histograms\n");
        fprintf(gnufile,"plot '%s/%s' u 1:2 smooth frequency with histeps"
                        " t \"%s\"\n",tmpdir,filename1,y_colname);
        fprintf(gnufile,"unset xlabel\n");
        fprintf(gnufile,"unset ylabel\n");
        fprintf(gnufile,"\n");
        fclose(gnufile);

        /* Calling gnuplot */
        sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname3);
        if (system(system_call)) {};

        /* Cleaning */
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname3);
        if (remove(tmpfilename)) {};
    }

    /* Cleaning */
    sc_basic_initstring(tmpfilename,SC_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,filename1);
    if (remove(tmpfilename)) {};
    if (rmdir(tmpdir)) {};

    cpl_array_delete(column_names1);

    return err_code;

}

/*****************************************************************************
 *                                                                           *
 *  Routine to plot a single spectrum with lines                             *
 *                                                                           *
 ****************************************************************************/



void sc_setplottags_single_spec_lines(cpl_parameterlist *plottags,
                                      const char *title,
                                      const char *x_label,
                                      const char *y_label,
                                      const cpl_parameterlist *parlist)
{
/*!
 * \callgraph
 *
 * Routine to set plot tags for single spectrum plot with detected lines.
 * Required input parameters are:
 *
 * \b INPUT:
 * \param plottags  Parameter list to be filled with tags
 * \param title     Plot title
 * \param x_label   Label string for x-axis
 * \param y_label   Label string for y-axis
 * \param parlist   Parameter list containing information from the parameter
 *                  file
 * \b OUTPUT:
 * \param plottags  Filled parameter list containing tags
 *
 * Base directory, output directory, output name, and gnuplot terminal type
 * are added automatically from parameter list.
 * NOTE: NULL is not accepted, use " " instead!
 *
 */

    cpl_parameter *p;
    const cpl_parameter *cp;
    char terminal[SC_MAXLEN], basedir[SC_MAXLEN];
    char outdir[SC_MAXLEN], outname[SC_MAXLEN];

    p = cpl_parameter_new_value("plot_title",CPL_TYPE_STRING,"","",title);
    cpl_parameterlist_append(plottags, p);
    p = cpl_parameter_new_value("xlabel",CPL_TYPE_STRING,"","",x_label);
    cpl_parameterlist_append(plottags, p);
    p = cpl_parameter_new_value("ylabel",CPL_TYPE_STRING,"","",y_label);
    cpl_parameterlist_append(plottags, p);

    cp = cpl_parameterlist_find_const(parlist,"plot_type");
    sprintf(terminal,"%s",cpl_parameter_get_string(cp));
    p = cpl_parameter_new_value("terminal",CPL_TYPE_STRING,"","",terminal);
    cpl_parameterlist_append(plottags, p);

    cp = cpl_parameterlist_find_const(parlist, "inst_dir");
    sprintf(basedir,"%s",cpl_parameter_get_string(cp));
    p = cpl_parameter_new_value("inst_dir",CPL_TYPE_STRING,"","",basedir);
    cpl_parameterlist_append(plottags, p);

    cp = cpl_parameterlist_find_const(parlist, "output_dir");
    sprintf(outdir,"%s",cpl_parameter_get_string(cp));
    p = cpl_parameter_new_value("output_dir",CPL_TYPE_STRING,"","",outdir);
    cpl_parameterlist_append(plottags, p);

    cp = cpl_parameterlist_find_const(parlist, "output_name");
    sprintf(outname,"%s",cpl_parameter_get_string(cp));
    p = cpl_parameter_new_value("output_name",CPL_TYPE_STRING,"","",outname);
    cpl_parameterlist_append(plottags, p);

}



/****************************************************************************/
/****************************************************************************/

cpl_error_code sc_plot_single_spec_with_lines(const cpl_table *spec,
                                              cpl_parameterlist *plottags)
{
/*!
 * \callgraph
 *
 * This program plots a single spectrum with detected lines from a CPL table.
 * One column (labelled "lambda") must contain wavelength information,
 * another one ("flux") gives the flux, and a the third column labelled
 * "class" must contain information whether the corresponding pixel belongs to
 * the continuum (value=0), a line (=1), a line peak (=2) or is the line peak
 * of an isolated line (=3). Line peaks are marked by a marker, i.e. a black
 * line, isolated peaks by a blue line. In addition, the wavelength of the
 * corresponding line peak pixel is added.
 *
 * \b INPUT:
 * \param spec          Cpl table containing spectrum
 * \param plottags      parameter list containing tags
 *
 * \b ERRORS:
 * - CPL_ERROR_NONE:           no error occurred
 *
 */

    FILE *specfile;                 /* output ASCII spec file ptr           */
    FILE *gnufile;                  /* gnuplot driverfile                   */
    char filename[SC_MAXLEN]="sci_single_plot_with_lines_spec.dat";
    char gnuname1[SC_MAXLEN]="single_plot_with_lines_gnufile_wxt.gnu";
    char gnuname2[SC_MAXLEN]="single_plot_with_lines_gnufile_ps.gnu";
    char gnuname3[SC_MAXLEN]="single_plot_with_lines_gnufile_x11.gnu";
    char ps_filename[SC_MAXLEN];

    char tmpdir[SC_MAXLEN];
    char tmpfilename[SC_MAXLEN];
    char plot_type[SC_MAXLEN];      /* plot type selection (plot_creation)  */

    char system_call[SC_MAXLEN];    /* system call                          */

    int len=0;                      /* table length                         */
    int run=0;                      /* runnning variable                    */
    int line_run=0;                 /* running number for lines to plot     */
    int n_lines=0;                  /* number of lines in spectrum          */
    //int ncol=0;                   /* number of columns in inputspectable  */
    //int dummy=0;                  /* dummy return value for system calls  */
    double lambda=0.;               /* wavelength of spectrum               */
    double line=0.;                 /* Line                                 */
    double y_value=0.;              /* y-value of plot (RAD/TRA)            */
    int dir_exist_flag=0;           /* Checking for existence of tmp dir    */
    /* Plot limits, line marker limits */
    double plot_xmin=0., plot_xmax=0.;
    double plot_ymin=0., plot_ymax=0.;
    double plot_dist_min=0., plot_dist_min_local=0., plot_dist_max=0.;
    double plot_label=0.;
    double corrfac=0.;              /* Correction factor for line labels    */

    /* plot tags */
    char plot_title[SC_MAXLEN];
    char x_colname[SC_MAXLEN], y_colname[SC_MAXLEN];

    /* directory + filename for output ps file */
    cpl_parameter *basedirpar, *outdirpar, *filenamepar, *par;

/*--------------------------------------------------------------------------*/
/* ------------------------------ INITIALISING -----------------------------*/
/*--------------------------------------------------------------------------*/

    /* Reading plot tags */
    par=cpl_parameterlist_find(plottags,"plot_title");
    sprintf(plot_title,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"xlabel");
    sprintf(x_colname,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"ylabel");
    sprintf(y_colname,"%s",cpl_parameter_get_string(par));
    par=cpl_parameterlist_find(plottags,"terminal");
    sprintf(plot_type,"%s",cpl_parameter_get_string(par));

    /* Checking table (=input spectra) properties */
    len=cpl_table_get_nrow(spec);   /* length of spectrum                   */
    //ncol=cpl_table_get_ncol(spec);  /* # of columns                       */
    plot_xmin=cpl_table_get_double(spec,"lambda",0,NULL);
    plot_xmax=cpl_table_get_double(spec,"lambda",len-1,NULL);
    corrfac=(plot_xmax-plot_xmin)/100;

    plot_ymin=cpl_table_get_column_min(spec,"flux");
    plot_ymax=1.5*cpl_table_get_column_max(spec,"flux");
    /* Line markers / Line labels */
    /* Length of line marker */
    plot_dist_max=plot_ymax-(plot_ymax-plot_ymin)/4;
    plot_dist_min=(plot_ymax-plot_ymin)/50;
    /* Location of labels */
    plot_label=plot_dist_max+plot_dist_min;//(plot_ymax-plot_ymin)/10;

/*    plot_dist_min=1.1*cpl_table_get_column_max(spec,"flux");*/

    /* Creating tmp-directory */
    sprintf(tmpdir,"__tmpDIRtmp__");
    dir_exist_flag=access(tmpdir,EXIST);
    if (dir_exist_flag == 0){
        cpl_msg_warning(cpl_func,"Directory %s already exists!",tmpdir);
    }
    else
    {
        if (mkdir(tmpdir,0777)) {};
    }

    /* writing .dat file containing spectrum information */
    sprintf(tmpfilename,"%s/%s",tmpdir,filename);
    specfile = fopen(tmpfilename,"w");
    for (run=0;run<len;run++)
    {
        lambda=cpl_table_get_double(spec,"lambda",run,NULL);
        y_value=cpl_table_get_double(spec,"flux",run,NULL);
        fprintf(specfile,"%5.6g\t%5.6g\n",lambda,y_value);
    }
    fclose(specfile);

    /* Counting lines */
    for(run=0;run<len-1;run++)
    {
            n_lines=n_lines+cpl_table_get_int(spec,"class",run,NULL);
    }
    if (n_lines == 0)
    {
        cpl_msg_warning(cpl_func,"Input spectrum contains no lines "
                                    "to plot!");
    }

    /* Creating wxt terminal gnuplot driver file */
    if ( (strcmp(plot_type,"W") == 0)  || (strcmp(plot_type,"WX")  == 0) ||
         (strcmp(plot_type,"XW")  == 0) )
    {
        sc_basic_initstring(tmpfilename,SC_MAXLEN);
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname1);
        gnufile = fopen(tmpfilename,"w");
        fprintf(gnufile,"set term wxt\n");
        fprintf(gnufile,"set termoption enhanced  \n");
        fprintf(gnufile,"set nokey\n");
        fprintf(gnufile,"set tmargin 2\n");
        fprintf(gnufile,"set bmargin 5\n");
        fprintf(gnufile,"set lmargin 13\n");
        fprintf(gnufile,"set rmargin 3\n");
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin,plot_xmax);
        fprintf(gnufile,"set yrange [ %g : %g ]\n",plot_ymin,plot_ymax);
        fprintf(gnufile,"set title \"%s\"\n",plot_title);
        fprintf(gnufile,"set xlabel '%s'\n",x_colname);
        fprintf(gnufile,"set ylabel '%s'\n",y_colname);
        fprintf(gnufile,"set parametric\n");
/*        fprintf(gnufile,"set trange [ %g : %g ]\n",plot_dist_min,
                                                   plot_dist_max);*/
        for(run=0;run<len-1;run++)
        {
            if (cpl_table_get_int(spec,"class",run,NULL) == 2)
            {
                plot_dist_min_local=cpl_table_get_double(spec,"flux",run,NULL)
                                    +plot_dist_min;
                line=cpl_table_get_double(spec,"lambda",run,NULL);
                fprintf(gnufile,"line%i=%g\n",line_run,line);
                fprintf(gnufile,"set arrow from %g,%g to %g,%g nohead "
                                "lw 0.5 lc -1\n",line,plot_dist_min_local,
                                line,plot_dist_max);
                fprintf(gnufile,"line%i=%g\n",line_run,line+corrfac);
                fprintf(gnufile,"set label \"{/=4 %2.5g}\" at "
                                "line%i,%g rotate\n",line,line_run,plot_label);
                line_run++;
            }
            if (cpl_table_get_int(spec,"class",run,NULL) == 3)
            {
                plot_dist_min_local=cpl_table_get_double(spec,"flux",run,NULL)
                                    +plot_dist_min;
                line=cpl_table_get_double(spec,"lambda",run,NULL);
                fprintf(gnufile,"line%i=%g\n",line_run,line);
                fprintf(gnufile,"set arrow from %g,%g to %g,%g nohead lw 1 "
                                "lc 3\n",line,plot_dist_min_local,
                                line,plot_dist_max);
                fprintf(gnufile,"line%i=%g\n",line_run,line+corrfac);
                fprintf(gnufile,"set label \"{/=4 %2.5g}\" at "
                                "line%i,%g rotate\n",line,line_run,plot_label);
                line_run++;
            }

        }
        fprintf(gnufile,"plot '%s/%s' using 1:2 with lines\n",
                tmpdir,filename);

        fprintf(gnufile,"\n");
        fclose(gnufile);

        /* Calling gnuplot */
        sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname1);
        if (system(system_call)) {};

        /* Cleaning */
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname1);
        if (remove(tmpfilename)) {};
    }

    /* Creating postscript terminal gnuplot driver file */
    basedirpar = cpl_parameterlist_find(plottags, "inst_dir");
    outdirpar = cpl_parameterlist_find(plottags, "output_dir");
    filenamepar = cpl_parameterlist_find(plottags, "output_name");

    sprintf(ps_filename,"%s/%s/%s_with_lines.ps",
            cpl_parameter_get_string(basedirpar),
            cpl_parameter_get_string(outdirpar),
            cpl_parameter_get_string(filenamepar));


    sc_basic_initstring(tmpfilename,SC_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,gnuname2);
    gnufile = fopen(tmpfilename,"w");
    fprintf(gnufile,"set term postscript enhanced color\n");
    fprintf(gnufile,"set output \"%s\"\n",ps_filename);
    fprintf(gnufile,"set nokey\n");
    fprintf(gnufile,"set tmargin 2\n");
    fprintf(gnufile,"set bmargin 5\n");
    fprintf(gnufile,"set lmargin 13\n");
    fprintf(gnufile,"set rmargin 3\n");
    fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin,plot_xmax);
    fprintf(gnufile,"set yrange [ %g : %g ]\n",plot_ymin,plot_ymax);
    fprintf(gnufile,"set title \"%s\"\n",plot_title);
    fprintf(gnufile,"set xlabel '%s'\n",x_colname);
    fprintf(gnufile,"set ylabel '%s'\n",y_colname);
    fprintf(gnufile,"set parametric\n");
    for(run=0;run<len-1;run++)
    {
        if (cpl_table_get_int(spec,"class",run,NULL) == 2)
        {
            plot_dist_min_local=cpl_table_get_double(spec,"flux",run,NULL)
                                +plot_dist_min;
            line=cpl_table_get_double(spec,"lambda",run,NULL);
            fprintf(gnufile,"line%i=%g\n",line_run,line);
            fprintf(gnufile,"set arrow from %g,%g to %g,%g nohead lw 0.5"
                            " lc -1\n",line,plot_dist_min_local,line,
                            plot_dist_max);
            fprintf(gnufile,"line%i=%g\n",line_run,line-corrfac/2);
            fprintf(gnufile,"set label \"{/=8    %2.5g}\" at "
                            "line%i,%g rotate\n",line,line_run,plot_label);
            line_run++;
        }
        if (cpl_table_get_int(spec,"class",run,NULL) == 3)
        {
            plot_dist_min_local=cpl_table_get_double(spec,"flux",run,NULL)
            +plot_dist_min;
            line=cpl_table_get_double(spec,"lambda",run,NULL);
            fprintf(gnufile,"line%i=%g\n",line_run,line);
            fprintf(gnufile,"set arrow from %g,%g to %g,%g nohead lw 1 "
                            "lc 3\n",line,plot_dist_min_local,
                            line,plot_dist_max);
            fprintf(gnufile,"line%i=%g\n",line_run,line-corrfac/2);
            fprintf(gnufile,"set label \"{/=8 %2.5g}\" at "
                            "line%i,%g rotate\n",line,line_run,plot_label);
            line_run++;
        }

    }
    fprintf(gnufile,"plot '%s/%s' using 1:2 with lines\n",
            tmpdir,filename);

    fprintf(gnufile,"\n");
    fclose(gnufile);

    /* Calling gnuplot */
    sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname2);
    if (system(system_call)) {};

    /* Cleaning */
    sprintf(tmpfilename,"%s/%s",tmpdir,gnuname2);
    if (remove(tmpfilename)) {};

    /* Creating x11 terminal gnuplot driver file */
    if ( (strcmp(plot_type,"X") == 0)  || (strcmp(plot_type,"WX")  == 0) ||
         (strcmp(plot_type,"XW")  == 0) )
    {
        sc_basic_initstring(tmpfilename,SC_MAXLEN);
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname3);
        gnufile = fopen(tmpfilename,"w");
        fprintf(gnufile,"set term x11\n");
        fprintf(gnufile,"set termoption enhanced  \n");
        fprintf(gnufile,"set nokey\n");
        fprintf(gnufile,"set tmargin 2\n");
        fprintf(gnufile,"set bmargin 5\n");
        fprintf(gnufile,"set lmargin 13\n");
        fprintf(gnufile,"set rmargin 3\n");
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin,plot_xmax);
        fprintf(gnufile,"set yrange [ %g : %g ]\n",plot_ymin,plot_ymax);
        fprintf(gnufile,"set title \"%s\"\n",plot_title);
        fprintf(gnufile,"set xlabel '%s'\n",x_colname);
        fprintf(gnufile,"set ylabel '%s'\n",y_colname);
        fprintf(gnufile,"set parametric\n");
        for(run=0;run<len-1;run++)
        {
            if (cpl_table_get_int(spec,"class",run,NULL) == 2)
            {
                plot_dist_min_local=cpl_table_get_double(spec,"flux",run,NULL)
                                    +plot_dist_min;
                line=cpl_table_get_double(spec,"lambda",run,NULL);
                fprintf(gnufile,"line%i=%g\n",line_run,line);
                fprintf(gnufile,"set arrow from %g,%g to %g,%g nohead lw 0.5 "
                                "lc -1\n",line,plot_dist_min_local,line,
                                plot_dist_max);
                fprintf(gnufile,"line%i=%g\n",line_run,line);
                fprintf(gnufile,"set label \"{/=4 %2.5g}\" at "
                                "line%i,%g rotate\n",line,line_run,plot_label);
                line_run++;
            }
            if (cpl_table_get_int(spec,"class",run,NULL) == 3)
            {
                plot_dist_min_local=cpl_table_get_double(spec,"flux",run,NULL)
                                    +plot_dist_min;
                line=cpl_table_get_double(spec,"lambda",run,NULL);
                fprintf(gnufile,"line%i=%g\n",line_run,line);
                fprintf(gnufile,"set arrow from %g,%g to %g,%g nohead lw 1 "
                                "lc 3\n",line,plot_dist_min_local,
                                line,plot_dist_max);
                fprintf(gnufile,"line%i=%g\n",line_run,line);
                fprintf(gnufile,"set label \"{/=4 %2.5g}\" at "
                                "line%i,%g rotate\n",line,line_run,plot_label);
                line_run++;
            }

        }
        fprintf(gnufile,"plot '%s/%s' using 1:2 with lines\n",
                tmpdir,filename);
        fprintf(gnufile,"\n");
        fclose(gnufile);

        /* Calling gnuplot */
        sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname3);
        if (system(system_call)) {};

        /* Cleaning */
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname3);
        if (remove(tmpfilename)) {};
    }

    /* Cleaning */
    sc_basic_initstring(tmpfilename,SC_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,filename);
    if (remove(tmpfilename)) {};
    if (rmdir(tmpdir)) {};

    return CPL_ERROR_NONE;
}

/**@}*/


