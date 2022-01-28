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
 * \defgroup sm_module2 Sky Model Module 2
 *
 * This module provides functions for the computation of the sky model for the
 * ESO ETC.
 */

/*!
 * \ingroup sm_module2
 */

/**@{*/

/*!
 * \file sm_skyemcomp.h
 *
 * Header for routines concerning the ETC sky model
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  02 Oct 2009
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

/*****************************************************************************
 *                                DEFINES                                    *
 ****************************************************************************/

#ifndef SM_SKYEMCOMP_H
#define SM_SKYEMCOMP_H

/* Definition of constants */

/*! List of data paths and file names  */
#define SM_FILENAMELIST "sm_filenames.dat"
/*! Line jumps for wavelength searches in LBLRTM/RFM library spectra */
#define SM_RRSTEP 1000
/*! \f$\mu{\rm m}\f$ */
#define SM_LAM_UNIT 1e-6
/*! Minimum wavelength for extraction of radiance \f$\to\f$ faster code */
#define SM_RADMINLAM 1.3 * 1e-6 / SM_LAM_UNIT
/*! average earth radius in km */
#define SM_ERAD 6371.
/*! Extension of Gaussian for airglow lines in \f$\sigma\f$ */
#define SM_SIGMAX 4.
/*! Minimum number of bins for \f$\sigma\f$ width of airglow lines */
#define SM_NSIGBIN 20

/*****************************************************************************
 *                                TYPEDEF                                    *
 ****************************************************************************/

/* Definition of structures */

/*!
 * Structure for sky emission parameters
 *
 * \param sm_h          observatory height in km
 * \param sm_hmin       Lower height limit in km
 * \param alt           altitude of object above horizon [0,90]
 * \param airmass       line-of-sight airmass
 * \param alpha         separation of sun and moon as seen from earth [0,180]
 * \param rho           separation of moon and object [0,180]
 * \param altmoon       altitude of moon above horizon [-90,90]
 * \param moondist      distance to moon (mean distance = 1; [0.945,1.055])
 * \param pres          pressure at observer altitude in hPa
 * \param ssa           single scattering albedo for aerosols [0,1]
 *                      (default: 0.97)
 * \param calcds        calculation of double scattering of moonlight
 *                      ('Y' or 'N')
 * \param o3column      relative UV/optical ozone column density (1 -> 258 DU)
 * \param moonscal      scaling factor for scattered moonlight (default: 1.0)
 * \param lon_ecl       heliocentric ecliptic longitude of object [-180,180]
 * \param lat_ecl       ecliptic latitude of object [-90,90]
 * \param ncomp         number of emissivity/temperature pairs
 *                      (\f$\le\f$ ::SM_MAXPAR)
 * \param eps[]         emissivity (factor for conversion black body \f$\to\f$
 *                      grey body)
 * \param temp[]        temperature in K
 * \param msolflux      monthly-averaged solar radio flux [sfu]
 * \param season        period of the year (1: Dec/Jan, ..., 6: Oct/Nov.;
 *                      0: entire year)
 * \param time          period of the night (x/3 of night, x = 1,2,3;
 *                      0: entire night)
 * \param vac_air       vac[uum] or air wavelengths
 * \param pwv           precipitable water vapour in mm (-1: bimonthly mean)
 * \param rtcode        radiative transfer code L(BLRTM) or R(FM) for
 *                      molecular spectra
 * \param resol         resolution of molecular spectra in library (crucial
 *                      for run time)
 * \param libpath       directory path for library of molecular spectra
 * \param libstruct     name of file containing the structure of the selected
                        LBLRTM/RFM library
 * \param libstruct1    name of structure file for time-dependent library
 * \param libstruct2    name of structure file for PWV-dependent library
 * \param datapath      directory path for input data files
 * \param solspecname   solar spectrum
 * \param mieextname    aerosol extinction table (optical depths;
 *                      NONE = parametrisation)
 * \param lunirrname    file for lunar irradiance model parameters
 * \param miephasename  file for Mie scattering phase functions
 * \param sscatcorname  file for multiple scattering correction of single
 *                      scattering calculations
 * \param dscatcorname  file for multiple scattering correction of double
 *                      scattering calculations
 * \param o3transname   file for UV/optical ozone transmission
 * \param starspecname  mean spectrum of scattered starlight
 * \param zodtabname    V-brightness of zodiacal light
 *            [\f$10^{-8}\,{\rm W\,m}^{-2}\,\mu{\rm m}^{-1}\,{\rm sr}^{-1}\f$]
 * \param linetabname   sky line table
 * \param vardatname    file for airglow lines scaling parameters
 * \param acontname     file for airglow continuum (scaling parameters
 *                      included)
 * \param incl          rules inclusion of sky emission components
 *        - format: "xxxxxxx" where x = "Y" (yes) or x = "N" (no)
 *        - position:
 *               - 1: scattered moonlight
 *               - 2: scattered starlight
 *               - 3: zodiacal light
 *               - 4: thermal emission by telescope/instrument
 *               - 5: molecular emission of lower atmosphere
 *               - 6: sky emission lines of upper atmosphere
 *               - 7: airglow continuum (residual continuum)
 */

typedef struct _smparmodel_ {
    double sm_h;
    double sm_hmin;
    double alt;
    double airmass;
    double alpha;
    double rho;
    double altmoon;
    double moondist;
    double pres;
    double ssa;
    char calcds[SM_LENLINE+1];
    double o3column;
    double moonscal;
    double lon_ecl;
    double lat_ecl;
    int ncomp;
    double eps[SM_MAXPAR];
    double temp[SM_MAXPAR];
    double msolflux;
    int season;
    int time;
    char vac_air[SM_LENLINE+1];
    double pwv;
    char rtcode[SM_LENLINE+1];
    double resol;
    char libpath[SM_MAXLEN+1];
    char libstruct[SM_LENLINE+1];
    char libstruct1[SM_LENLINE+1];
    char libstruct2[SM_LENLINE+1];
    char datapath[SM_MAXLEN+1];
    char solspecname[SM_LENLINE+1];
    char mieextname[SM_LENLINE+1];
    char lunirrname[SM_LENLINE+1];
    char miephasename[SM_LENLINE+1];
    char sscatcorname[SM_LENLINE+1];
    char dscatcorname[SM_LENLINE+1];
    char o3transname[SM_LENLINE+1];
    char starspecname[SM_LENLINE+1];
    char zodtabname[SM_LENLINE+1];
    char linetabname[SM_LENLINE+1];
    char vardatname[SM_LENLINE+1];
    char acontname[SM_LENLINE+1];
    char incl[SM_LENLINE+1];
} smparmodel;

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code sm_etc_calcmodel(cpl_table *skytable,
                                const cpl_parameterlist *params);
cpl_error_code sm_etc_getparams(smparmodel *modelpar,
                                const cpl_parameterlist *params);
const cpl_parameter
*sm_etc_parameterlist_find_const(const cpl_parameterlist *self,
                                 const char *name, cpl_type type);
cpl_error_code sm_etc_splitstring(cpl_array *outval, char *instring);
cpl_error_code sm_etc_readfilenames(smparmodel *modelpar,
                                    const cpl_parameterlist *params);
cpl_error_code sm_etc_writetable(cpl_table *skytable, const smspec *radiance,
                                 const smspec *transmission,
                                 const cpl_table *radcomp,
                                 const cpl_table *transcomp);
cpl_error_code sm_comp_skyemcomp(smspec *radiance, smspec *transmission,
                                 cpl_table *radcomp, cpl_table *transcomp,
                                 const smparmodel modelpar);
cpl_error_code sm_comp_extrapoltrans(smspec *spec, const smparmodel modelpar,
                                     const double lim[2]);
double sm_comp_alttoairmass(const double alt);
cpl_error_code sm_comp_getmolspec(smspec *spec,
                                  const cpl_parameterlist *libfilepar,
                                  const double lim[2]);
cpl_error_code sm_comp_readlibstruct(cpl_parameterlist *libfilepar,
                                     const smparmodel modelpar,
                                     const char *spectype);
double sm_comp_vactoair_single(const double lam);
cpl_error_code sm_comp_vactoair_spec(smspec *spec);
double sm_comp_airtovac_single(const double lam);
cpl_error_code sm_comp_airtovac_spec(smspec *spec);
cpl_error_code sm_comp_scaleo3(smspec *molabstrans, smspec *abstrans,
                               smspec *o3trans, const smspec *libtrans,
                               const smparmodel modelpar);
double sm_comp_calcrayleighscat1(const double lam, const double pres,
                                const double sm_h);
cpl_error_code sm_comp_calcrayleighscat(smspec *rscattrans,
                                        const smparmodel modelpar);
double sm_comp_calcmiescat1(const double lam);
cpl_error_code sm_comp_calcmiescat(smspec *mscattrans,
                                    const smparmodel modelpar);
cpl_error_code sm_comp_scaletranscurv(smspec *trans,
                                      const smparmodel modelpar);
cpl_error_code sm_comp_writetranscomp(cpl_table *transcomp,
                                      const char *colname,
                                      const smspec *trans,
                                      const smparmodel modelpar);
cpl_error_code sm_comp_lunskybright_o(smspec *spec, const smparmodel modelpar,
                                      const smspec *solspec,
                                      const smspec *molabstrans,
                                      const smspec *rscattrans,
                                      const smspec *mscattrans);
cpl_error_code sm_comp_lunskybright(smspec *spec, const smparmodel modelpar,
                                    const smspec *solspec,
                                    const smspec *abstrans,
                                    const smspec *o3trans,
                                    const smspec *rscattrans,
                                    const smspec *mscattrans);
cpl_error_code sm_comp_getmoonalbedo(smspec *albedo,
                                     const smparmodel modelpar);
cpl_error_code sm_comp_scatstarlight(smspec *spec, const smparmodel modelpar,
                                     const smspec *molabstrans);
cpl_error_code sm_comp_zodskybright(smspec *spec, const smparmodel modelpar,
                                    const smspec *solspec,
                                    const smspec *molabstrans,
                                    const smspec *rscattrans,
                                    const smspec *mscattrans);
cpl_error_code sm_comp_telem(smspec *spec, const smparmodel modelpar);
cpl_error_code sm_comp_extrapolrad(smspec *spec, const smspec *trans,
                                   const smparmodel modelpar,
                                   const double lim[2]);
cpl_error_code sm_comp_getlinespec(smspec *spec, const smparmodel modelpar,
                                   const smspec *rscattrans,
                                   const smspec *mscattrans);
cpl_error_code sm_comp_readvarpar(cpl_table *varpar,
                                  const smparmodel modelpar);
cpl_error_code sm_comp_scalelinetab(cpl_table *linedat, const smspec *linetab,
                                    const cpl_table *varpar);
cpl_error_code sm_comp_convertlinetab(smspec *spec, const cpl_table *linedat);
cpl_error_code sm_comp_convertlinetabo(smspec *spec,
                                       const cpl_table *linedat);
cpl_error_code sm_comp_airglowcont(smspec *spec, const smparmodel modelpar,
                                   const smspec *molabstrans,
                                   const smspec *rscattrans,
                                   const smspec *mscattrans);

#endif /* SM_SKYEMCOMP_H */

#ifdef __cplusplus
}
#endif

/**@}*/
