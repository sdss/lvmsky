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
 * \file sc_skycorr.h
 *
 * Header for top-level sky correction routines
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  18 Feb 2011
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
#include <sc_par.h>
#include <sc_readspec.h>
#include <sc_lines.h>
#include <sc_specdiss.h>
#include <sc_contsub.h>
#include <sc_fwhmest.h>
#include <sc_weights.h>
#include <sc_mpfit.h>

/*****************************************************************************
 *                                DEFINES                                    *
 ****************************************************************************/

#ifndef SC_SKYCORR_H
#define SC_SKYCORR_H

/* Definition of constants */

/*! Required relative accuracy for comparison of the wavelength grids of
    science and sky spectrum */
#define SC_GRIDTOL 1e-2

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code sc_skycorr(const char *parfile_in);
cpl_error_code sc_skycorr_makeoutputdir(const cpl_parameterlist *parlist);
cpl_error_code sc_skycorr_readspec(cpl_table *scispec, cpl_table *skyspec,
                                   cpl_parameterlist *parlist);
cpl_error_code sc_skycorr_subcont(cpl_table *spec, cpl_table *linetab,
                                  cpl_parameterlist *parlist,
                                  const cpl_table *groups);
cpl_error_code sc_skycorr_fit(cpl_table *scispec, cpl_table *skyspec,
                              cpl_table *groups, cpl_parameterlist *parlist);
cpl_error_code sc_skycorr_subsky(cpl_table *scispec);
cpl_error_code sc_skycorr_writespec(cpl_table *scispec,
                                    cpl_parameterlist *parlist);
cpl_error_code sc_skycorr_plot(cpl_table *scispec,
                               const cpl_parameterlist *parlist);

#endif /* SC_SKYCORR_H */

#ifdef __cplusplus
}
#endif

/**@}*/
