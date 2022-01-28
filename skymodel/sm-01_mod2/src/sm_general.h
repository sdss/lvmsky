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
 * \defgroup sm_general General Sky Model Module
 *
 * This module provides basic definitions and functions for both sky model
 * modules (\em sm_module1 and \em sm_module2).
 */

/*!
 * \ingroup sm_general
 */

/**@{*/

/*!
 * \file sm_general.h
 *
 * Header for basic routines used for the sky model
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  22 Sep 2009
 * \date   06 Oct 2015
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
#include <unistd.h>
#include <sys/stat.h>
#include <errno.h>

/* CPL library */

#include <cpl.h>

/*****************************************************************************
 *                                DEFINES                                    *
 ****************************************************************************/

#ifndef SM_GENERAL_H
#define SM_GENERAL_H

/* Basic definitions */

/*! Macro for finding the minimum of two values */
#define SM_MIN(a, b) ((a) > (b) ? (b) : (a))
/*! Macro for finding the maximum of two values */
#define SM_MAX(a, b) ((a) < (b) ? (b) : (a))

/*! Boolean type (F = false, T = true) */

#ifndef SM_BOOL
#define SM_BOOL
typedef enum SM_BOOL {SM_F, SM_T} smbool;
#endif

/*! Boolean false */
#ifdef F
#undef F
#endif
#define F SM_F

/*! Boolean true */
#ifdef T
#undef T
#endif
#define T SM_T

/* Definitions related to sm_general library */

/* Definition of constants */

/*! Maximum number of error codes */
#define SM_NMAXERR 100
/*! Maximum number of characters per line */
#define SM_LENLINE 160
/*! Maximum number of string characters */
#define SM_MAXLEN 4000
/*! Maximum number of space-, tab-, or '='-separated strings per line
    (see ::sm_param_read) */
#define SM_MAXPAR 21
/*! Required relative accuracy for data comparisons */
#define SM_TOL 1e-7
/*! steradians in \f${\rm arcsec}^2\f$ */
#define SM_SR_IN_ARCSEC2 4.254517e+10
/*! Universal gas constant [J/(molK)] */
#define SM_GAS_CONST 8.31447
/*! Molar mass of dry air [kg/mol] */
#define SM_MOL_AIR_DRY 0.0289644
/*! Earth-surface gravitational acceleration [m/\f${\rm s}^2\f$] */
#define SM_GRAV_ACC 9.80665

/* Definition of error codes and corresponding standard error messages */

/*! Enumeration structure for sky model related error codes */

typedef enum _sm_error_code_ {
    SM_ERROR_FOPEN = CPL_ERROR_EOL + 11, // "File opening failed"
    SM_ERROR_UFS = CPL_ERROR_EOL + 12, // "Unexpected file structure"
    SM_ERROR_IFE = CPL_ERROR_EOL + 13, // "Invalid file name extension"
    SM_ERROR_NDA = CPL_ERROR_EOL + 20, // "No data"
    SM_ERROR_INSUFF_DATA = CPL_ERROR_EOL + 21, // "Insufficient data points"
    SM_ERROR_IDG = CPL_ERROR_EOL + 22, // "Inconsistent data grids"
    SM_ERROR_IDR = CPL_ERROR_EOL + 23, // "Invalid data range"
    SM_ERROR_IOD = CPL_ERROR_EOL + 24, // "Invalid order of data points"
    SM_ERROR_IIP = CPL_ERROR_EOL + 30, // "Invalid input parameter(s)"
    SM_ERROR_IOV = CPL_ERROR_EOL + 31, // "Invalid object value(s)"
    SM_ERROR_IOS = CPL_ERROR_EOL + 32, // "Invalid object structure"
    SM_ERROR_SUBROUTINE = CPL_ERROR_EOL + 40, // "Error in subroutine"
    // ERRORS RELATED TO access(), mkdir(), chdir()
    SM_ERROR_ACCES       = CPL_ERROR_EOL + 50, // access(), mkdir(), chdir()
    SM_ERROR_LOOP        = CPL_ERROR_EOL + 51, // access(), mkdir(), chdir()
    SM_ERROR_NAMETOOLONG = CPL_ERROR_EOL + 52, // access(), mkdir(), chdir()
    SM_ERROR_NOENT       = CPL_ERROR_EOL + 53, // access(), mkdir(), chdir()
    SM_ERROR_NOTDIR      = CPL_ERROR_EOL + 54, // access(), mkdir(), chdir()
    SM_ERROR_ROFS        = CPL_ERROR_EOL + 55, // access(), mkdir()
    SM_ERROR_FAULT       = CPL_ERROR_EOL + 56, // access(), mkdir(), chdir()
    SM_ERROR_INVAL       = CPL_ERROR_EOL + 57, // access()
    SM_ERROR_IO          = CPL_ERROR_EOL + 58, // access(),          chdir()
    SM_ERROR_NOMEM       = CPL_ERROR_EOL + 59, // access(), mkdir(), chdir()
    SM_ERROR_TXTBSY      = CPL_ERROR_EOL + 60, // access()
    SM_ERROR_EXIST       = CPL_ERROR_EOL + 61, //           mkdir()
    SM_ERROR_NOSPC       = CPL_ERROR_EOL + 62, //           mkdir()
    SM_ERROR_PERM        = CPL_ERROR_EOL + 63, //           mkdir()
    // General errors
    SM_ERROR_BADUSERINPUT = CPL_ERROR_EOL + 70,
    SM_ERROR_LINK        = CPL_ERROR_EOL + 71,
    SM_ERROR_RFM         = CPL_ERROR_EOL + 81, // general RFM error
    SM_ERROR_UNDEF = CPL_ERROR_EOL + 80
} sm_error_code;

/* Aliases for sky model related error codes */

#define SM_ERROR_FOF SM_ERROR_FOPEN
#define SM_ERROR_ISM SM_ERROR_NOMEM
#define SM_ERROR_EIS SM_ERROR_SUBROUTINE
#define SM_ERROR_ISD SM_ERROR_INSUFF_DATA

/* Standard messages for sky model related errors */

#define SM_ERROR_FOPEN_TXT "File opening failed"
#define SM_ERROR_UFS_TXT "Unexpected file structure"
#define SM_ERROR_IFE_TXT "Invalid file name extension"
#define SM_ERROR_BDR_TXT "Bad directory"
#define SM_ERROR_NDA_TXT "No data"
#define SM_ERROR_ISD_TXT "Insufficient data points"
#define SM_ERROR_IDG_TXT "Inconsistent data grids"
#define SM_ERROR_IDR_TXT "Invalid data range"
#define SM_ERROR_IOD_TXT "Invalid order of data points"
#define SM_ERROR_IIP_TXT "Invalid input parameter(s)"
#define SM_ERROR_IOV_TXT "Invalid object value(s)"
#define SM_ERROR_IOS_TXT "Invalid object structure"
#define SM_ERROR_SUBROUTINE_TXT "Error in subroutine"
// ERRORS RELATED TO access(), mkdir(), chdir()
#define SM_ERROR_ACCES_TXT       "Permission denied"
#define SM_ERROR_LOOP_TXT        "Too many symbolic links"
#define SM_ERROR_NAMETOOLONG_TXT "Pathname too long"
#define SM_ERROR_NOENT_TXT       "File/dir does not exist"
#define SM_ERROR_NOTDIR_TXT      "Component used as directory in pathname " \
                                 "is not a directory"
#define SM_ERROR_ROFS_TXT        "Write permission requested for file/dir " \
                                 "on read-only file system"
#define SM_ERROR_FAULT_TXT       "Pathname points outside accessible " \
                                 "address space"
#define SM_ERROR_INVAL_TXT       "Mode was incorrectly specified"
#define SM_ERROR_IO_TXT          "I/O error occurred"
#define SM_ERROR_NOMEM_TXT       "Insufficient memory"
#define SM_ERROR_TXTBSY_TXT      "Write access requested to executable " \
                                 "which is being executed"
#define SM_ERROR_EXIST_TXT       "File/dir already exists"
#define SM_ERROR_NOSPC_TXT       "No space left on device"
#define SM_ERROR_PERM_TXT        "File system does not support creation of " \
                                 "directories"
#define SM_ERROR_LINK_TXT        "Could not create symbolic link"
#define SM_ERROR_UNDEF_TXT       "Undefined error"

#define SM_ERROR_FOF_TXT SM_ERROR_FOPEN_TXT
#define SM_ERROR_ISM_TXT SM_ERROR_NOMEM_TXT
#define SM_ERROR_EIS_TXT SM_ERROR_SUBROUTINE_TXT

/*****************************************************************************
 *                                TYPEDEF                                    *
 ****************************************************************************/

/*!
 * Structure for data points
 *
 * \param lam     wavelength
 * \param flux    flux
 * \param dflux1  (lower) flux error
 * \param dflux2  upper flux error
 */

typedef struct _smdat_ {
    double lam;
    double flux;
    double dflux1;
    double dflux2;
} smdat;

/*!
 * Structure for spectra
 *
 * \param type   number of columns (1 = no flux, 2 = no errors,
 *               3 = symmetric error, 4 = lower and upper error)
 * \param n      number of wavelengths
 * \param *dat   vector of data points defined by ::smdat structure
 */

typedef struct _smspec_ {
    int type;
    int n;
    smdat *dat;
} smspec;

/*!
 * Structure for a parameter value in different data types
 *
 * \param c  string (maximum length ::SM_LENLINE)
 * \param i  integer number
 * \param d  float number of double precision
 * \param n  number of value (for arrays; counted backwards)
 */

typedef struct _smparam_ {
    char c[SM_LENLINE+2];
    int i;
    double d;
    int n;
} smparam;

/*!
 * Structure for coordinate grids
 *
 * \param nx     number of \e x coordinate values
 * \param ny     number of \e y coordinate values
 * \param *xpos  vector of \e x coordinate values
 * \param *ypos  vector of \e y coordinate values
 * \param **val  matrix of data values for the coordinate grid
 */

typedef struct _smgrid_ {
    int nx;
    int ny;
    double *xpos;
    double *ypos;
    double **val;
} smgrid;

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code sm_spec_malloc(smspec *spec, const int size);
cpl_error_code sm_spec_create(smspec *outspec, const double limlam[2],
                              const double dlam);
cpl_error_code sm_spec_read(smspec *spec, const char *filename);
cpl_error_code sm_spec_readrange(smspec *spec, const char *filename,
                                 const double limlam[2], const int step);
cpl_error_code sm_spec_readcpl(smspec *spec, const cpl_table *cpltab);
cpl_error_code sm_spec_readfits(smspec *spec, const char *filename);
cpl_error_code sm_spec_readfitsrange(smspec *spec, const char *filename,
                                     const double limlam[2], const int step);
cpl_error_code sm_spec_copy(smspec *outspec, const smspec *inspec);
cpl_error_code sm_spec_compgrids(const smspec *spec1, const smspec *spec2);
cpl_error_code sm_spec_join(smspec *spec, const smspec *errfunc,
                            const int errflag);
cpl_error_code sm_spec_split(smspec *spec, smspec *errfunc,
                             const int errflag);
cpl_error_code sm_spec_changetype(smspec *spec, const int type);
cpl_error_code sm_spec_getval(double *flux, const smspec *spec,
                              const double lam);
cpl_error_code sm_spec_scalerange(smspec *spec, const double limlam[2],
                                  const char op, const double c);
cpl_error_code sm_spec_scale(smspec *spec, const char op, const double c);
cpl_error_code sm_spec_modval(smspec *spec, const double lam, const char op,
                              const double c);
cpl_error_code sm_spec_calc(smspec *spec, const char op,
                            const smspec *opspec);
cpl_error_code sm_spec_funct(smspec *spec, const char *funct,
                             const char baselab);
cpl_error_code sm_spec_functnoerr(smspec *spec, const char *funct);
cpl_error_code sm_spec_convunits(smspec *spec, const double factor,
                                 const int lamexp);
cpl_error_code sm_spec_changegrid(smspec *spec, const double factor,
                                  const char *scale);
cpl_error_code sm_spec_minmax(double limflux[2], const smspec *spec,
                              const double limlam[2]);
cpl_error_code sm_spec_average(double mean[3], const smspec *spec,
                               const double limlam[2]);
cpl_error_code sm_spec_write(const smspec *spec, const char *filename);
cpl_error_code sm_spec_print(const smspec *spec);
cpl_error_code sm_spec_writecpl(cpl_table *cpltab, const smspec *spec);
cpl_error_code sm_spec_writecplcolumn(cpl_table *cpltab, const char *colname,
                                      const smspec *spec, const int datatype);
cpl_error_code sm_spec_writefits(const smspec *spec, const char *filename);
cpl_error_code sm_spec_free(smspec *spec);
cpl_error_code sm_spec_extract(smspec *outspec, const smspec *inspec,
                               const double limlam[2]);
cpl_error_code sm_spec_rebin(smspec *outspec, const smspec *inspec);
cpl_error_code sm_spec_interpol(smspec *outspec, const smspec *inspec);
cpl_error_code sm_spec_convolve(smspec *spec, const int nkpix,
                                const double *kernel);

cpl_error_code sm_param_read(FILE *stream, smparam par[]);
cpl_error_code sm_param_readcheck(FILE *stream, smparam par[],
                                  const char *parname, const int npar);
cpl_error_code sm_param_check(smparam par[], const char *parname,
                              const int npar);

cpl_error_code sm_grid_malloc(smgrid *xy, const int nx, const int ny);
cpl_error_code sm_grid_read(smgrid *xy, const char *filename);
cpl_error_code sm_grid_extract(double *outval, const smgrid *xy,
                               const double x0, const double y0);
cpl_error_code sm_grid_write(const smgrid *xy, const char *filename);
cpl_error_code sm_grid_print(const smgrid *xy);
cpl_error_code sm_grid_free(smgrid *xy);

cpl_error_code sm_basic_chdir(const char *dir);
cpl_error_code sm_basic_access(const char *pathname, const int mode);
cpl_error_code sm_basic_mkdir(const char *dir, const mode_t mode);
cpl_error_code sm_basic_createdir(const char *dir, const mode_t mode);
void sm_basic_initstring(char *str, const long n);
cpl_boolean sm_basic_isnumber(char *str);
char *sm_basic_replacestring(char *instring, char *oldsubstr,
                             char *newsubstr);
char *sm_basic_rmcntrl(char *str);
void sm_basic_rmcntrl_inplace(char *str);
char *sm_basic_strtrim(char *str);
void sm_basic_strtrim_inplace(char *str);
void sm_basic_terminatestring(char *str);
cpl_error_code sm_basic_interpollin(const double *x_out, double *y_out,
                                    const long n_out, const double *x_ref,
                                    const double *y_ref, const long n_ref);

#endif /* SM_GENERAL_H */

#ifdef __cplusplus
}
#endif

/**@}*/
