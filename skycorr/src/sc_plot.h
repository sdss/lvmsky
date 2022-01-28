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
 * \file sc_plot.h
 *
 * \brief Header for plotting library
 *
 * \author Wolfgang Kausch & ESO In-Kind Team Innsbruck
 * \since   02 Feb 2011
 * \date    28 Aug 2013
 */

/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <sc_basic.h>

/*****************************************************************************
 *                                 DEFINES                                   *
 ****************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef SC_PLOT_H
#define SC_PLOT_H

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************************
 *                                 TYPEDEF                                   *
 ****************************************************************************/

typedef struct _scplottags_ {
    char title[SC_MAXLEN];
    char xlabel[SC_MAXLEN];
    char ylabel[SC_MAXLEN];
} scplottags;

/*****************************************************************************
 *                                 GLOBALS                                   *
 ****************************************************************************/

/*****************************************************************************
 *                                 PROTOTYPES                                *
 ****************************************************************************/

cpl_error_code sc_plot_single_spec(const cpl_table *spec,
                            cpl_parameterlist *plottags);

/*cpl_error_code sc_plot_double_spec(const cpl_table *spec1,
                            const cpl_table *spec2,
                            cpl_parameterlist *parlist);*/
cpl_error_code sc_plot_double_spec(const cpl_table *spec1,
                            const cpl_table *spec2,
                            cpl_parameterlist *plottags);

cpl_error_code sc_overplot_spec(const cpl_table *spec,
                            cpl_parameterlist *plottags);

cpl_error_code sc_plot_hist(const cpl_table *histdat,
                            cpl_parameterlist *plottags);
cpl_error_code sc_plot_single_spec_with_lines(const cpl_table *spec,
                            cpl_parameterlist *plottags);

/* Plot tag routines */
void sc_setplottags_single_spec(cpl_parameterlist *plottags,
                            const char *x_column,
                            const char *y_column,
                            const char *title,
                            const char *x_label,
                            const char *y_label,
                            const cpl_parameterlist *parlist);
void sc_setplottags_double_spec(cpl_parameterlist *plottags,
                            const char *x_column1,
                            const char *y_column1,
                            const char *title1,
                            const char *x_label1,
                            const char *y_label1,
                            const char *x_column2,
                            const char *y_column2,
                            const char *title2,
                            const char *x_label2,
                            const char *y_label2,
                            const cpl_parameterlist *parlist);
void sc_setplottags_overplot_spec(cpl_parameterlist *plottags,
                            const char *x_column,
                            const char *y_column1,
                            const char *y_column2,
                            const char *specname1,
                            const char *specname2,
                            const char *title,
                            const char *x_label,
                            const char *y_label,
                            const cpl_parameterlist *parlist);
void sc_setplottags_hist(cpl_parameterlist *plottags, const char *title,
                            const char *x_label, const char *y_label,
                            const cpl_parameterlist *parlist);
void sc_setplottags_single_spec_lines(cpl_parameterlist *plottags,
                            const char *title,
                            const char *x_label,
                            const char *y_label,
                            const cpl_parameterlist *parlist);

#ifdef __cplusplus
}
#endif

#endif
