/*
 *  This file is part of the Sky Background software package.
 *  Copyright (C) 2009-2018 European Southern Observatory
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

/*!
 * \defgroup sm_module2_run Sky Model Module 2 Running
 *
 * This module provides functions for running the sky model module 2
 * (\em sm_module2).
 */

/*!
 * \ingroup sm_module2_run
 */

/**@{*/

/*!
 * \file sm_run.h
 *
 * Header for routines related to the ETC sky model stand-alone programme
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  27 May 2010
 * \date   06 Oct 2015
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

/* Sky model headers */

#include <sm_general.h>
#include <sm_skyemcomp.h>

/*****************************************************************************
 *                                DEFINES                                    *
 ****************************************************************************/

#ifndef SM_RUN_H
#define SM_RUN_H

/* Definition of constants */

/*! Threshold for recalculation of kernel depending on relative change of
    wavelength */
#define SM_LIMRELLAMVAR 0.01
/*! for integral of Voigt profile */
#define SM_BINS_PER_FWHM 50

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code sm_run_readmodelpar(cpl_parameterlist *params,
                                   const char *parfile);
cpl_error_code sm_run_readinstrupar(cpl_parameterlist *lamgrid,
                                    cpl_parameterlist *params,
                                    const char *parfile);
cpl_error_code sm_run_createtable(cpl_table *skytable,
                                  const cpl_parameterlist *lamgrid);
cpl_error_code sm_run_calcmodel(cpl_table *skytable,
                                cpl_parameterlist *params);
cpl_error_code sm_run_convolvespec(cpl_table *skytable,
                                   cpl_parameterlist *params);
cpl_error_code sm_run_calckernel(cpl_array *kernel,
                                 const cpl_parameterlist *params);
cpl_error_code sm_run_calcvoigtkernel(cpl_array *kernel,
                                      const cpl_parameterlist *params);
cpl_error_code sm_run_calcboxkernel(cpl_array *kernel,
                                    const cpl_parameterlist *params);
cpl_error_code sm_run_convolvearray(cpl_array *array,
                                    const cpl_array *kernel);
cpl_error_code sm_run_writekernel(const cpl_array *kernel,
                                  const cpl_parameterlist *params);
cpl_error_code sm_run_convolvewindow(cpl_table *outtable,
                                     const cpl_table *intable,
                                     const int range[2],
                                     const cpl_array *kernel);
cpl_error_code sm_run_copytablewindow(cpl_table *outtable,
                                      const cpl_table *intable,
                                      const int range[2]);
cpl_error_code sm_run_writefits(const cpl_table *skytable,
                                const char *radfile, const char *transfile);

#endif /* SM_RUN_H */

#ifdef __cplusplus
}
#endif

/**@}*/
