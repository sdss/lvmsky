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
 * \file sc_modsky.h
 *
 * Header for routines related to the modification of the sky spectrum
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  20 Feb 2011
 * \date   24 Feb 2014
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

#ifndef SC_MODSKY_H
#define SC_MODSKY_H

/* Definition of constants */

/*! This factor times the maximum error in the sky spectrum results in a
    reference error for the substitution of weight = 0 */
#define SC_RELMAXERR 10.

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code sc_modsky(cpl_table *scispec, cpl_table *skyspec,
                         cpl_table *fitpar, const cpl_vector *sinc,
                         const cpl_parameterlist *parlist);
cpl_error_code sc_modsky_modlines(cpl_table *skyspec, cpl_table *fitpar,
                                  const char grouptype);
cpl_error_code sc_modsky_moderrors(cpl_table *skyspec);
cpl_error_code sc_modsky_modwavegrid(cpl_table *skyspec, cpl_table *fitpar);
cpl_error_code sc_modsky_rebin(cpl_table *scispec, cpl_table *skyspec,
                               const cpl_vector *sinc,
                               const cpl_parameterlist *parlist);
cpl_error_code sc_modsky_getpixelshifts(cpl_table *skyspec,
                                        const cpl_table *scispec);
cpl_error_code sc_modsky_sincrebin(cpl_table *scispec, const char *scicol,
                                   const cpl_table *skyspec,
                                   const char *skycol,
                                   const cpl_vector *sinc);
cpl_error_code sc_modsky_calcsinckernel(cpl_array *kernel,
                                        const double shift);
cpl_error_code sc_modsky_calcdev(cpl_table *scispec);

#endif /* SC_MODSKY_H */

#ifdef __cplusplus
}
#endif

/**@}*/
