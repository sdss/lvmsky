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
 * \file sc_contsub.h
 *
 * Header for routines related to continuum identification, interpolation, and
 * subtraction
 *
 * \author Marco Barden, Stefan Noll, & ESO In-Kind Team Innsbruck
 *
 * \date   12 Jun 2014
 */

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************************
 *                                INCLUDES                                   *
 ****************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* Sky correction headers */

#include <sc_basic.h>

/*****************************************************************************
 *                                DEFINES                                    *
 ****************************************************************************/

#ifndef SC_CONTSUB_H
#define SC_CONTSUB_H

/* Definition of constants */

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code sc_contsub(cpl_table *spec);
cpl_error_code sc_contsub_identcont(cpl_table *spec, const cpl_table *groups,
                                    const cpl_table *linetab,
                                    const double fluxlim,
                                    const cpl_parameterlist *parlist);
int sc_contsub_check(cpl_table *spec);
cpl_error_code sc_contsub_skipthermir(cpl_table *spec);

#endif /* SC_CONTSUB_H */

#ifdef __cplusplus
}
#endif

/**@}*/
