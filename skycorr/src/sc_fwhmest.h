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
 * \file sc_fwhmest.h
 *
 * Header for continuum interpolation and subtraction routine
 *
 * \author Marco Barden, Stefan Noll, & ESO In-Kind Team Innsbruck
 *
 * \date   12 May 2013
 *
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

/*****************************************************************************
 *                                DEFINES                                    *
 ****************************************************************************/

#ifndef SC_FWHMEST_H
#define SC_FWHMEST_H

/*****************************************************************************
 *                                TYPEDEF                                    *
 ****************************************************************************/

/* Definition of structures */

/*!
 * Structure for y(x) fitting
 *
 * \param x      independent variable of model
 * \param y      measured "y" values
 * \param y_err  measurement uncertainty in y
 */

typedef struct _scxydata_ {
    cpl_array *x;
    cpl_array *y;
    cpl_array *y_err;
} scxydata;

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code sc_fwhmest(double *fwhm, double *rms, cpl_table *spec,
                          cpl_table *lines, const cpl_parameterlist *parlist);
int sc_fwhmest_gaussfunc(int n_data, int n_pars, double *p, double *deviates,
                         double **derivs, void *pdata);

#endif // SC_FWHMEST_H

#ifdef __cplusplus
}
#endif

/**@}*/
