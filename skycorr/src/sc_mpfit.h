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
 * \file sc_mpfit.h
 *
 * Header for routines related to the handling of CMPFIT
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  17 Feb 2011
 * \date   24 Feb 2014
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

/* Sky correction headers */

#include <sc_basic.h>
#include <sc_modsky.h>

/*****************************************************************************
 *                                DEFINES                                    *
 ****************************************************************************/

#ifndef SC_MPFIT_H
#define SC_MPFIT_H

/* Definition of constants */

/*! Default error value (taken if no specific error is available) */
#define SC_DEFERRVAL 9.99
/*! Minimum correction factor for line groups */
#define SC_CORRFAC_MIN 0.01
/*! Maximum correction factor for line groups */
#define SC_CORRFAC_MAX 100.
/*! Minimum number of valid lines for a save line group */
#define SC_MINNLIN 3
/*! Maximum relative deviation between correction factor of a system A group
    with less than ::SC_MINNLIN lines and the system mean */
#define SC_MAXRELFAC 2.
/*! Maximum relative uncertainty of fit parameters for line flux correction
    (otherwise the initial parameter value is taken) */
#define SC_MAXPARERR 0.5

/*****************************************************************************
 *                                GLOBALS                                    *
 ****************************************************************************/

/* Declaration of global variables */

/*! Number of fitting function calls */
extern int nfev;
/*! Last modsky call? */
extern cpl_boolean lastcall;

/*****************************************************************************
 *                                TYPEDEF                                    *
 ****************************************************************************/

/* Definition of structures */

/*!
 * Structure for holding the spectral data and the driver file parameters
 *
 * \param scispec  CPL table with science spectrum
 * \param skyspec  CPL table with sky spectrum and line group weights
 * \param fitpar   CPL table with fit parameters
 * \param sinc     CPL vector with damped sinc kernel
 * \param parlist  general CPL parameter list
 */

typedef struct _scvars_ {
    cpl_table *scispec;
    cpl_table *skyspec;
    cpl_table *fitpar;
    cpl_vector *sinc;
    const cpl_parameterlist *parlist;
} scvars;

/*!
 * Structure for holding the fit parameters and their constraints
 *
 * \param n     number of parameters
 * \param p     parameter vector for fitting
 * \param pars  structure for parameter constraints
 */

typedef struct _scpars_ {
    int n;
    double *p;
    mp_par *pars;
} scpars;

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code sc_mpfit(mp_result *result, cpl_table *scispec,
                        cpl_table *skyspec, cpl_table *fitpar,
                        const cpl_parameterlist *parlist);
int sc_mpfit_calcdev(int m, int n, double *p, double *dy, double **dvec,
                     void *vars);
cpl_error_code sc_mpfit_modinitpar(cpl_table *fitpar,
                                   cpl_table *scispec,
                                   cpl_table *skyspec,
                                   const cpl_parameterlist *parlist);
cpl_error_code sc_mpfit_setpar(scpars *fitpars,
                               const cpl_table *fitpar, const char fittype);
cpl_error_code sc_mpfit_allocmempar(scpars *fitpars, const int npar);
cpl_error_code sc_mpfit_freemempar(scpars *fitpars);
cpl_error_code sc_mpfit_allocmemresult(mp_result *result, const int m,
                                       const int n);
cpl_error_code sc_mpfit_initresult(mp_result *result, const int m,
                                   const int n);
cpl_error_code sc_mpfit_copyresult(mp_result *outresult,
                                   const mp_result *inresult);
cpl_error_code sc_mpfit_freememresult(mp_result *result);
cpl_error_code sc_mpfit_substbadfitpar(cpl_table *fitpar,
                                       const cpl_table *initfitpar,
                                       const cpl_parameterlist *parlist);
cpl_error_code sc_mpfit_writeresults(const mp_result *result,
                                     const cpl_table *scispec,
                                     const cpl_table *fitpar,
                                     const cpl_parameterlist *parlist,
                                     const int calls, const double runtime);

#endif /* SC_MPFIT_H */

#ifdef __cplusplus
}
#endif

/**@}*/
