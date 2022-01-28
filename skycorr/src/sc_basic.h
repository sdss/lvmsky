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
 * \defgroup sky_correction Sky Correction Code
 *
 * This module provides functions for the correction of sky emission.
 */

/**@{*/

/*!
 * \file sc_basic.h
 *
 * Basic header for sky correction code
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  15 Feb 2011
 * \date   12 Jun 2014
 */

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************************
 *                                INCLUDES                                   *
 ****************************************************************************/

/* Config header */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* C standard libraries */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <errno.h>
#include <unistd.h>
#include <assert.h>
#include <sys/dir.h>
#include <sys/stat.h>
#include <sys/wait.h>

/* CPL library */

#include <cpl.h>

/* MPFIT */

#include <mpfit.h>

/*****************************************************************************
 *                                DEFINES                                    *
 ****************************************************************************/

#ifndef SC_BASIC_H
#define SC_BASIC_H

/* Basic definitions */

/*! Macro for finding the minimum of two values */
#define SC_MIN(a, b) ((a) > (b) ? (b) : (a))
/*! Macro for finding the maximum of two values */
#define SC_MAX(a, b) ((a) < (b) ? (b) : (a))

/* Definition of constants */

/*! Maximum number of characters per line */
#define SC_LENLINE 160
/*! Maximum number of string characters */
#define SC_MAXLEN 1024
/*! Required relative accuracy for data comparisons */
#define SC_TOL 1e-7
/*! Maximum number of space- or tab-separated strings per line
    (see ::sc_basic_readline) */
#define SC_MAXPAR 21
/*! Lower limit for wavelength range dominated by emission lines originating
    in the lower atmosphere in \f${\mu}\f$m */
#define SC_THERMIRLIM 2.31

/* Default column names */

/*! Default wavelength column name */
#define SC_DEFLAMCOL "LAMBDA"
/*! Default flux column name */
#define SC_DEFFLUXCOL "FLUX"
/*! Default flux error column name */
#define SC_DEFDFLUXCOL "DFLUX"
/*! Default mask column name */
#define SC_DEFMASKCOL "MASK"

/* Parameters for damped sinc kernel */

/*! Radius of precomputed damped sinc kernel in pixels. */
#define SC_SINCRAD_PRECOMP 6.0
/*! Number of points per unit for pre-calculated damped sinc kernel */
#define SC_SINCNBIN 10000
/*! Damping constant for damped sinc kernel in pixels */
#define SC_SINCDAMP 3.25

/* Definition of error codes and corresponding standard error messages */

/*! Enumeration structure for sky model related error codes */

typedef enum _sc_error_code_ {
    SC_ERROR_FOPEN = CPL_ERROR_EOL + 11, // "File opening failed"
    SC_ERROR_UFS = CPL_ERROR_EOL + 12, // "Unexpected file structure"
    SC_ERROR_IFE = CPL_ERROR_EOL + 13, // "Invalid file name extension"
    SC_ERROR_NDA = CPL_ERROR_EOL + 20, // "No data"
    SC_ERROR_INSUFF_DATA = CPL_ERROR_EOL + 21, // "Insufficient data points"
    SC_ERROR_IDG = CPL_ERROR_EOL + 22, // "Inconsistent data grids"
    SC_ERROR_IDR = CPL_ERROR_EOL + 23, // "Invalid data range"
    SC_ERROR_IOD = CPL_ERROR_EOL + 24, // "Invalid order of data points"
    SC_ERROR_IIP = CPL_ERROR_EOL + 30, // "Invalid input parameter(s)"
    SC_ERROR_IOV = CPL_ERROR_EOL + 31, // "Invalid object value(s)"
    SC_ERROR_IOS = CPL_ERROR_EOL + 32, // "Invalid object structure"
    SC_ERROR_SUBROUTINE = CPL_ERROR_EOL + 40, // "Error in subroutine"
    // ERRORS RELATED TO access(), mkdir(), chdir()
    SC_ERROR_ACCES       = CPL_ERROR_EOL + 50, // access(), mkdir(), chdir()
    SC_ERROR_LOOP        = CPL_ERROR_EOL + 51, // access(), mkdir(), chdir()
    SC_ERROR_NAMETOOLONG = CPL_ERROR_EOL + 52, // access(), mkdir(), chdir()
    SC_ERROR_NOENT       = CPL_ERROR_EOL + 53, // access(), mkdir(), chdir()
    SC_ERROR_NOTDIR      = CPL_ERROR_EOL + 54, // access(), mkdir(), chdir()
    SC_ERROR_ROFS        = CPL_ERROR_EOL + 55, // access(), mkdir()
    SC_ERROR_FAULT       = CPL_ERROR_EOL + 56, // access(), mkdir(), chdir()
    SC_ERROR_INVAL       = CPL_ERROR_EOL + 57, // access()
    SC_ERROR_IO          = CPL_ERROR_EOL + 58, // access(),          chdir()
    SC_ERROR_NOMEM       = CPL_ERROR_EOL + 59, // access(), mkdir(), chdir()
    SC_ERROR_TXTBSY      = CPL_ERROR_EOL + 60, // access()
    SC_ERROR_EXIST       = CPL_ERROR_EOL + 61, //           mkdir()
    SC_ERROR_NOSPC       = CPL_ERROR_EOL + 62, //           mkdir()
    SC_ERROR_PERM        = CPL_ERROR_EOL + 63, //           mkdir()
    // General errors
    SC_ERROR_BADUSERINPUT = CPL_ERROR_EOL + 70,
    SC_ERROR_LINK        = CPL_ERROR_EOL + 71,
    SC_ERROR_CD          = CPL_ERROR_EOL + 72, // could not change directory
    SC_ERROR_GETCWD      = CPL_ERROR_EOL + 73, // could not get current
                                               // working directory
    SC_ERROR_RFM         = CPL_ERROR_EOL + 81, // general RFM error
    SC_ERROR_UNDEF = CPL_ERROR_EOL + 80
} sc_error_code;

/* Aliases for sky model related error codes */

#define SC_ERROR_FOF SC_ERROR_FOPEN
#define SC_ERROR_ISM SC_ERROR_NOMEM
#define SC_ERROR_EIS SC_ERROR_SUBROUTINE
#define SC_ERROR_ISD SC_ERROR_INSUFF_DATA

/* Standard messages for sky model related errors */

#define SC_ERROR_FOPEN_TXT "File opening failed"
#define SC_ERROR_UFS_TXT "Unexpected file structure"
#define SC_ERROR_IFE_TXT "Invalid file name extension"
#define SC_ERROR_BDR_TXT "Bad directory"
#define SC_ERROR_NDA_TXT "No data"
#define SC_ERROR_ISD_TXT "Insufficient data points"
#define SC_ERROR_IDG_TXT "Inconsistent data grids"
#define SC_ERROR_IDR_TXT "Invalid data range"
#define SC_ERROR_IOD_TXT "Invalid order of data points"
#define SC_ERROR_IIP_TXT "Invalid input parameter(s)"
#define SC_ERROR_IOV_TXT "Invalid object value(s)"
#define SC_ERROR_IOS_TXT "Invalid object structure"
#define SC_ERROR_SUBROUTINE_TXT "Error in subroutine"
// ERRORS RELATED TO access(), mkdir(), chdir()
#define SC_ERROR_ACCES_TXT       "Permission denied"
#define SC_ERROR_LOOP_TXT        "Too many symbolic links"
#define SC_ERROR_NAMETOOLONG_TXT "Pathname too long"
#define SC_ERROR_NOENT_TXT       "File/dir does not exist"
#define SC_ERROR_NOTDIR_TXT      "Component used as directory in pathname " \
                                 "is not a directory"
#define SC_ERROR_ROFS_TXT        "Write permission requested for file/dir " \
                                 "on read-only file system"
#define SC_ERROR_FAULT_TXT       "Pathname points outside accessible " \
                                 "address space"
#define SC_ERROR_INVAL_TXT       "Mode was incorrectly specified"
#define SC_ERROR_IO_TXT          "I/O error occurred"
#define SC_ERROR_NOMEM_TXT       "Insufficient memory"
#define SC_ERROR_TXTBSY_TXT      "Write access requested to executable " \
                                 "which is being executed"
#define SC_ERROR_EXIST_TXT       "File/dir already exists"
#define SC_ERROR_NOSPC_TXT       "No space left on device"
#define SC_ERROR_PERM_TXT        "File system does not support creation of " \
                                 "directories"
#define SC_ERROR_LINK_TXT        "Could not create symbolic link"
#define SC_ERROR_UNDEF_TXT       "Undefined error"
#define SC_ERROR_CD_TXT          "Could not change directory"
#define SC_ERROR_GETCWD_TXT      "Could not get current working directory"

#define SC_ERROR_FOF_TXT SC_ERROR_FOPEN_TXT
#define SC_ERROR_ISM_TXT SC_ERROR_NOMEM_TXT
#define SC_ERROR_EIS_TXT SC_ERROR_SUBROUTINE_TXT

/*****************************************************************************
 *                                TYPEDEF                                    *
 ****************************************************************************/

/* Definition of structures */

/*!
 * Structure for a parameter value in different data types
 *
 * \param c  string (maximum length ::SC_LENLINE)
 * \param i  integer number
 * \param d  float number of double precision
 */

typedef struct _scpar_ {
    char c[SC_LENLINE+2];
    int i;
    double d;
} scpar;

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code sc_basic_readline(FILE *stream, scpar par[], int *npar);
double sc_basic_mjd2fracyear(double mjd);
double sc_basic_fracyear2date(int *year, int *month, int *day,
                              int *hh, int *mm, double *ss,
                              const double *fracyear);
cpl_error_code sc_basic_rebin(cpl_table *outspec, const char *outlam,
                              const char *outflux, const cpl_table *inspec,
                              const char *inlam, const char *influx);
cpl_error_code sc_basic_convolve(cpl_table *spec, const char *colname,
                                 const cpl_array *kernel);
cpl_error_code sc_basic_convolvewindow(cpl_array *convflux,
                                       const cpl_array *flux,
                                       const int range[2],
                                       const cpl_array *kernel);
cpl_error_code sc_basic_filtermedian(cpl_table *spec, const char *colname,
                                     const int npix);
cpl_error_code sc_basic_linecount(long *n_lines, FILE *fp);
cpl_error_code sc_basic_rhum2ppmv_old(const cpl_array *temp,
                                      const cpl_array *pres,
                                      const cpl_array *rhum,
                                      cpl_array *ppmv);
void sc_basic_rhum2ppmv(double *ppmv, const double *tem, const double *p,
                        const double *hum);
cpl_error_code sc_basic_planck(cpl_array *bb, const cpl_array *wavelength,
                               const double temp);
void sc_basic_dirslash(char *dir);
void sc_basic_abspath(char *out, const char *dir, const char *cwd);
cpl_error_code sc_basic_access(const char *pathname, const int mode);
cpl_error_code sc_basic_greg2jd(long *jd, const int year, const int month,
                                const int day);
cpl_error_code sc_basic_jd2greg(int *year, int *month, int *day,
                                const long jd);
cpl_error_code sc_basic_absfile(char *absfilename, const char *filename);
cpl_boolean sc_basic_parameterlists_same(cpl_parameterlist *list1,
                                         cpl_parameterlist *list2);
cpl_error_code sc_basic_status2txt(char *msg, const int status);
cpl_error_code sc_basic_clipmean(double *mean, double *rms,
                                 cpl_array *arr, const cpl_boolean clip);

cpl_error_code sc_basic_col2arr(cpl_array *arr, const cpl_table *tab,
                                const char *colname);
cpl_error_code sc_basic_arr2col(cpl_table *tab, const char *colname,
                                const cpl_array *arr);
int sc_basic_sortarr_double(const void *p1, const void *p2);
cpl_error_code sc_basic_sortarr(cpl_array *arr, cpl_type type);
cpl_error_code sc_basic_copytable_content(cpl_table *outtab,
                                          const cpl_table *intab);
cpl_error_code sc_basic_copytable_full(cpl_table *outtab,
                                       const cpl_table *intab);
void sc_basic_initstring(char *str, const long n);
void sc_basic_terminatestring(char *str);
char *sc_basic_strtrim(char *str);
void sc_basic_strtrim_inplace(char *str);
cpl_boolean sc_basic_isnumber(char *str);
cpl_boolean sc_basic_isinteger(char *str);
cpl_error_code sc_basic_interpollin(const double *x_out, double *y_out,
                                    const long n_out, const double *x_ref,
                                    const double *y_ref, const long n_ref,
                                    const int extrapolate);
extern int sc_basic_gaussfunc(double *fgauss, const double *xgauss,
                              const int n_data, const double *par);
cpl_error_code sc_basic_getfilename(char *dir, char *filename, char *suffix,
                                    const char *path);
cpl_error_code sc_basic_getmaskval_vector(double maskval[2],
                                          const cpl_vector *vector);
cpl_error_code sc_basic_getmaskval_image(double maskval[2],
                                         const cpl_image *image);
cpl_error_code sc_basic_calcsinc(cpl_vector *sinc);

#endif /* SC_BASIC_H */

#ifdef __cplusplus
}
#endif

/**@}*/
