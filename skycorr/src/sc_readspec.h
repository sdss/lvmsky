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
 * \file sc_readspec.h
 *
 * Header for routines related to the reading of observing data
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  16 Feb 2011
 * \date   12 Jul 2013
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
#include <sc_par.h>
#include <sc_conv.h>

/*****************************************************************************
 *                                DEFINES                                    *
 ****************************************************************************/

#ifndef SC_READSPEC_H
#define SC_READSPEC_H

/* Definition of constants */

/*! Number of FITS header keywords defined in sc_readspec_header */
#define SC_NKEY 3

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code sc_readspec(cpl_table *spec, cpl_parameterlist *parlist);
cpl_error_code sc_readspec_fits(cpl_table *spec,
                                const cpl_parameterlist *parlist,
                                const char *filename);
cpl_error_code sc_readspec_header(cpl_parameterlist *parlist,
                                  const char *filename);
cpl_error_code sc_readspec_setweight(cpl_table *spec,
                                     const cpl_parameterlist *parlist);

#endif /* SC_READSPEC_H */

#ifdef __cplusplus
}
#endif

/**@}*/
