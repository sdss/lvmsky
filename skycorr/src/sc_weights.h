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
 * \file sc_weights.h
 *
 * Header for library for preparing the line group weights of each pixel of
 * the input sky spectrum
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  14 Feb 2011
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

#ifndef SC_WEIGHTS_H
#define SC_WEIGHTS_H

/* Definition of constants */

/*! Initial depth of array columns for line group weights */
#define SC_COLDEPTH 10
/*! Extra pixels at both sides of the airglow model spectrum
    (pixel scale of input sky spectrum) */
#define SC_EXTRACOVER 5
/*! Oversampling factor for airglow model spectrum */
#define SC_SAMPFAC 5.
/*! Threshold for recalculation of kernel depending on relative change of
    wavelength */
#define SC_LIMRELLAMVAR 0.01
/*! Bins per FWHM for integral of Gaussian */
#define SC_BINS_PER_FWHM 200.
/*! Radius of Gaussian kernel in FWHM/2 */
#define SC_KERNFAC 3.

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code sc_weights(cpl_table *skyspec, cpl_table *fitpar,
                          cpl_table *groups, cpl_parameterlist *parlist);
cpl_error_code sc_weights_initfitpar(cpl_table *fitpar,
                                    cpl_parameterlist *parlist,
                                    const cpl_table *groups);
cpl_error_code sc_weights_getpixcontrib(cpl_table *skyspec, cpl_table *fitpar,
                                        cpl_table *groups,
                                        const cpl_parameterlist *parlist);
cpl_error_code sc_weights_createwavegrid(cpl_table *groupspec,
                                         const cpl_table *skyspec);
cpl_error_code sc_weights_getgrouplines(cpl_table *groupspec,
                                        cpl_table *fitpar, cpl_table *groups,
                                        const cpl_parameterlist *parlist,
                                        const char grouptype,
                                        const int group);
cpl_error_code sc_weights_convolve(cpl_table *groupspec, const double fwhm,
                                   const double reflam,
                                   const double speedpar);
cpl_error_code sc_weights_calckernel(cpl_array *kernel, const double fwhm);
cpl_error_code sc_weights_rebinspec(cpl_table *skyspec,
                                    const cpl_table *groupspec);
cpl_error_code sc_weights_fillgrouparrays(cpl_table *skyspec,
                                           const char grouptype,
                                           const int group);
cpl_error_code sc_weights_normpixcontrib(cpl_table *skyspec,
                                         const char grouptype);
cpl_error_code sc_weights_getdomgroups(cpl_table *skyspec,
                                       const cpl_parameterlist *parlist,
                                       const char grouptype);

#endif /* SC_WEIGHTS_H */

#ifdef __cplusplus
}
#endif

/**@}*/
