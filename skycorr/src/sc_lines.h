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
 * \file sc_lines.h
 *
 * Header for routines handling the airglow line list
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  17 Jul 2013
 * \date   17 Jul 2013
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

#ifndef SC_LINES_H
#define SC_LINES_H

/* Definition of constants */

/*! Average Earth radius in km */
#define SC_ERAD 6371.

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code sc_lines(cpl_table *groups, cpl_parameterlist *parlist);
cpl_error_code sc_lines_readlist(cpl_table *groups,
                                 const cpl_parameterlist *parlist);
cpl_error_code sc_lines_getmodelbins(cpl_parameterlist *parlist);
cpl_error_code sc_lines_readsolflux(cpl_parameterlist *parlist);
cpl_error_code sc_lines_readvarpar(cpl_array *varpar,
                                   cpl_parameterlist *parlist);
cpl_error_code sc_lines_scalelines(cpl_table *groups, const cpl_array *varpar,
                                   const cpl_parameterlist *parlist);
cpl_error_code sc_lines_vactoair(cpl_table *groups,
                                 const cpl_parameterlist *parlist);

#endif /* SC_LINES_H */

#ifdef __cplusplus
}
#endif

/**@}*/
