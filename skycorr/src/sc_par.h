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
 * \file sc_par.h
 *
 * \brief Header for routines for handling parameter files
 *
 * \author Wolfgang Kausch, Stefan Noll, & ESO In-Kind Team Innsbruck
 * \since  02 Feb 2011
 * \date   17 Apr 2013
 */

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

/* Config header */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* Sky correction headers */

#include <sc_basic.h>

/*****************************************************************************
 *                                 DEFINES                                   *
 ****************************************************************************/

#ifndef SC_PAR_H
#define SC_PAR_H

/* Definition of constants */

/*****************************************************************************
 *                                 PROTOTYPES                                *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code sc_par_readfile(cpl_parameterlist *parlist,
                               const char *filename);
cpl_error_code sc_par_writefile(const cpl_parameterlist *parlist,
                                const char *filename);
cpl_error_code sc_par_addkeywords(cpl_propertylist *header,
                                  const cpl_parameterlist *parlist);

#endif /* SC_PAR_H */

#ifdef __cplusplus
}
#endif

/**@}*/
