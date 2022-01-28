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

/*!
 * \file sc_specdiss.h
 *
 * \brief Header for spectral analysis library
 *
 * \author Wolfgang Kausch, Stefan Noll, & ESO In-Kind Team Innsbruck
 * \since  01 Feb 2011
 * \date   08 Sep 2013
 *
 */

/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

/* Config header */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* Sky correction headers */

#include <sc_basic.h>
#include <sc_contsub.h>

/*****************************************************************************
 *                                 DEFINES                                   *
 ****************************************************************************/

#ifndef SC_SPECDISS_H
#define SC_SPECDISS_H

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************************
 *                                 TYPEDEF                                   *
 ****************************************************************************/

/*****************************************************************************
 *                                 GLOBALS                                   *
 ****************************************************************************/

/*****************************************************************************
 *                                 PROTOTYPES                                *
 ****************************************************************************/

/* Main routines */
cpl_error_code sc_specdiss_find_emissionlines(cpl_table *input_spectrum,
                           cpl_table *linetab,
                           cpl_parameterlist *parlist,
                           const cpl_table *groups);
cpl_error_code sc_specdiss_find_isolatedlines(cpl_table *linetab,
                           cpl_table *input_spectrum,
                           cpl_parameterlist *parlist);
cpl_error_code sc_specdiss_merge_speclinelist(cpl_table *spec,
                           const cpl_table *groups,
                           const cpl_parameterlist *parlist);
cpl_error_code sc_specdiss_identify_airglowlines(cpl_table *spec);


/* Internal routines */
cpl_error_code sc_specdiss_create_hist(cpl_table *hist,
                                      const cpl_table *input_table,
                                      char *col_name, const int n_bins);
cpl_error_code sc_specdiss_find_valuesoutside(cpl_table *input_table,
                                              char *newcolname,
                                              char *colname,
                                              const double lowlim,
                                              const double uplim);
cpl_error_code sc_specdiss_find_valuesabove(cpl_table *input_table,
                                            char *newcolname,
                                            char *colname,
                                            const double lowlim);
cpl_error_code sc_specdiss_find_localmaxima(cpl_table *input_table,
                                            char *newcolname,
                                            char *colname,
                                            const double lower_threshold,
                                            const double upper_threshold);
int sc_specdiss_count_lines(const cpl_table *intable);
void sc_specdiss_init_linetab(cpl_table *linetab);

#endif /* SC_SPECDISS_H */

#ifdef __cplusplus
}
#endif

