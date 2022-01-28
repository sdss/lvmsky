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
 * \file sc_conv.h
 *
 * Header for routines related to file conversion
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  16 Apr 2013
 * \date   31 Jan 2014
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

/*****************************************************************************
 *                                DEFINES                                    *
 ****************************************************************************/

#ifndef SC_CONV_H
#define SC_CONV_H

/* Definition of constants */

/*! Fraction of left or right edge pixels */
#define SC_EDGEFRAC 0.025
/*! Factor (to be multiplied by extreme fluxes calculated excluding the edges)
    for rejecting unreliable edge fluxes */
#define SC_MAXFLUXFAC 10.

/*****************************************************************************
 *                                TYPEDEF                                    *
 ****************************************************************************/

/*!
 * Structure for data of FITS tables
 *
 * \param next    number of FITS extensions
 * \param **head  array of CPL property lists of size \e next + 1
 * \param **tab   array of CPL tables of size \e next + 1
 */

typedef struct _sctarr_ {
    int next;
    cpl_propertylist **head;
    cpl_table **tab;
} sctarr;

/*!
 * Structure for data of 1D FITS images
 *
 * \param next    number of FITS extensions
 * \param **head  array of CPL property lists of size \e next + 1
 * \param **vec   array of CPL vectors of size \e next + 1
 */

typedef struct _scvarr_ {
    int next;
    cpl_propertylist **head;
    cpl_vector **vec;
} scvarr;

/*!
 * Structure for data of 2D FITS images
 *
 * \param next    number of FITS extensions
 * \param **head  array of CPL property lists of size \e next + 1
 * \param **ima   array of CPL images of size \e next + 1
 */

typedef struct _sciarr_ {
    int next;
    cpl_propertylist **head;
    cpl_image **ima;
} sciarr;

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code sc_conv_preptable(const char *parfile_in);
cpl_error_code sc_conv_readfile(sctarr *tabdat,
                                const cpl_parameterlist *parlist);
cpl_error_code sc_conv_checkfitsformat(int *fitsformat, const char *filename);
cpl_error_code sc_conv_setcolnames(cpl_array *colnames,
                                   const cpl_parameterlist *parlist);
cpl_error_code sc_conv_setextnames_varr(cpl_table *extnames,
                                        const scvarr *vecdat,
                                        const cpl_parameterlist *parlist);
cpl_error_code sc_conv_setextnames_iarr(cpl_table *extnames,
                                        const sciarr *imadat,
                                        const cpl_parameterlist *parlist);
cpl_error_code sc_conv_varr2tarr(sctarr *tabdat, const scvarr *vecdat,
                                 const cpl_table *extnames);
double sc_conv_getwcskey(const cpl_propertylist *plist, const char *key);
cpl_error_code sc_conv_modtable(sctarr *tabdat,
                                const cpl_parameterlist *parlist);
cpl_error_code sc_conv_writetable(const sctarr *tabdat,
                                  const cpl_parameterlist *parlist);
cpl_error_code sc_conv_writeresults(const char *parfile_in);
cpl_error_code sc_conv_readresults(cpl_table *results,
                                   const cpl_parameterlist *parlist);
cpl_error_code sc_conv_writefile(cpl_table *results,
                                 const cpl_parameterlist *parlist);
cpl_error_code sc_conv_readprepfits(sctarr *tabdat,
                                    const cpl_parameterlist *parlist);
cpl_error_code sc_conv_results2tarr(sctarr *tabdat, cpl_table *results);
cpl_error_code sc_conv_erasemaskcol(sctarr *tabdat,
                                    const cpl_parameterlist *parlist);
cpl_error_code sc_conv_resultstarr2varr(scvarr *vecdat,
                                        const cpl_table *extnames,
                                        const sctarr *tabdat);
cpl_error_code sc_conv_getmaskval(double maskval[2], const sctarr *tabdat,
                                  const cpl_table *extnames);
cpl_error_code sc_conv_extract1d(const char *parfile_in,
                                 const char *scifilename,
                                 const char *skyfilename,
                                 const double minrelprofflux,
                                 const double lowpixfrac);
cpl_error_code sc_conv_iarr2varr_sci(scvarr *vecdat, cpl_vector *selpix,
                                     const sciarr *imadat,
                                     const cpl_table *extnames,
                                     const double minrelprofflux);
cpl_error_code sc_conv_iarr2varr_sky(scvarr *vecdat, const sciarr *imadat,
                                     const cpl_table *extnames,
                                     const double lowpixfrac,
                                     const cpl_vector *selpix);
cpl_error_code sc_conv_ascii_read(sctarr *tabdat, const char *filename,
                                  const cpl_array *colnames);
cpl_error_code sc_conv_ascii_write(const char *filename,
                                   const sctarr *tabdat);
cpl_error_code sc_conv_tarr_init(sctarr *tabdat, const int next);
cpl_error_code sc_conv_tarr_read(sctarr *tabdat, const char *filename);
cpl_error_code sc_conv_tarr_write(const char *filename, const sctarr *tabdat);
cpl_error_code sc_conv_tarr_delete(sctarr *tabdat);
cpl_error_code sc_conv_varr_init(scvarr *vecdat, const int next);
cpl_error_code sc_conv_varr_read(scvarr *vecdat, const char *filename);
cpl_error_code sc_conv_varr_write(const char *filename, const scvarr *vecdat);
cpl_error_code sc_conv_varr_delete(scvarr *vecdat);
cpl_error_code sc_conv_iarr_init(sciarr *imadat, const int next);
cpl_error_code sc_conv_iarr_read(sciarr *imadat, const char *filename);
cpl_error_code sc_conv_iarr_write(const char *filename, const sciarr *imadat);
cpl_error_code sc_conv_iarr_delete(sciarr *imadat);

#endif /* SC_CONV_H */

#ifdef __cplusplus
}
#endif

/**@}*/
