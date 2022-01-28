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
 * \ingroup sm_module2
 */

/**@{*/

/*!
 * \file sm_scatmoonlight.h
 *
 * Header for routines for calculation of scattered moonlight
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  05 Jun 2012
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

#ifndef SM_SCATMOONLIGHT_H
#define SM_SCATMOONLIGHT_H

/* Definition of constants */

/*! Mode for calculation of scattered moonlight (0 = most accurate,
    3 = fastest) */
#define SM_SCATMODE 0
/*! Mode for calculation of multiple scattering correction for moonlight
    (0 = most accurate, 3 = fastest) */
#define SM_MSCATMODE 0
/*! Number of zenith distances for double scattering coordinate grid
    (large grid) */
#define SM_NZ 40
/*! Number of zenith distances for double scattering coordinate grid
    (small grid) */
#define SM_NZD 22
/*! Number of path elements for integration (single scattering) */
#define SM_NS 20
/*! Number of path elements for integration (double scattering) */
#define SM_NSD 12
/*! Scaling factor in km for unitless logarithmic bins */
#define SM_LNSCALE 0.5
/*! Radius of Earth in km */
#define SM_R 6371.
/*! Top of atmosphere in km */
#define SM_HMAX 200. /* Staude (1975) */
/*! Density of atmosphere at surface in \f${\rm cm}^{-3}\f$ */
#define SM_N0_R 2.67e19 /* Staude (1975) */
/*! Scale height of atmosphere in km */
#define SM_H0_R 7.99 /* Staude (1975) */
/*! Density of aerosols at surface in \f${\rm cm}^{-3}\f$ */
#define SM_N0_M 1.11e4 /* Staude (1975) */
/*! Scale height of aerosol number density in km */
#define SM_H0_M 1.2 /* Staude (1975) */
/*! Reference wavelength for effective absorption airmass calculation in
    \f$\mu{\rm m}\f$ */
#define SM_REFLAM0 1.0
/*! Maximum zenith distance for the approximation \f$X = cot(z)\f$ in deg */
#define SM_ZMAX 75.
/*! Radius in deg for refinement of coordinate grid */
#define SM_DMAX 10.

/*****************************************************************************
 *                                TYPEDEF                                    *
 ****************************************************************************/

/* Definition of structures */

/*!
 * Structure for results of scattering calculations
 *
 * \param iscatr   intensity scattered into line of sight by single Rayleigh
 *                 scattering
 * \param iscatm   intensity scattered into line of sight by single Mie
 *                 scattering
 * \param iscatrr  intensity scattered into line of sight by double Rayleigh
 *                 scattering
 * \param iscatrm  intensity scattered into line of sight by Mie and then
 *                 Rayleigh scattering
 * \param iscatrm  intensity scattered into line of sight by Rayleigh and
 *                 then Mie scattering
 * \param iscatmm  intensity scattered into line of sight by double Mie
 *                 scattering
 * \param iscatrg  intensity scattered into line of sight by ground
 *                 and then Rayleigh scattering
 * \param iscatmg  intensity scattered into line of sight by ground and
 *                 then Mie scattering
 * \param iscat1s  intensity scattered into line of sight by single scattering
 * \param iscat2s  intensity scattered into line of sight by double scattering
 * \param iscat    intensity scattered into line of sight
 * \param it1r     as \e iscatr but with zero extinction
 * \param it1m     as \e iscatm but with zero extinction
 * \param it1rr    as \e iscatrr but with zero extinction
 * \param it1rm    as \e iscatrm but with zero extinction
 * \param it1mr    as \e iscatmr but with zero extinction
 * \param it1mm    as \e iscatmm but with zero extinction
 * \param it1rg    as \e iscatrg but with zero extinction
 * \param it1mg    as \e iscatmg but with zero extinction
 * \param it1_1s   as \e iscat1s but with zero extinction
 * \param it1_2s   as \e iscat2s but with zero extinction
 * \param it1      as \e iscat but with zero extinction
 * \param xeff     effective airmass for single + double scattering
 */

typedef struct _smscatres_ {
    double iscatr;
    double iscatm;
    double iscatrr;
    double iscatrm;
    double iscatmr;
    double iscatmm;
    double iscatrg;
    double iscatmg;
    double iscat1s;
    double iscat2s;
    double iscat;
    double it1r;
    double it1m;
    double it1rr;
    double it1rm;
    double it1mr;
    double it1mm;
    double it1rg;
    double it1mg;
    double it1_1s;
    double it1_2s;
    double it1;
    double xeff;
} smscatres;

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */
cpl_error_code sm_scat_moon(smspec *spec,
		            const double sm_h, const double sm_hmin,
		            const double z, const double zmoon,
                            const double rho, const double ssa,
                            const smspec *molabstrans,
                            const smspec *rscattrans,
                            const smspec *mscattrans, const smgrid *miephase,
                            const smgrid *sscatcor, const smspec *dscatcor,
                            const char *calcds);
cpl_error_code sm_scat_createlambdagrid(smspec *lamgrid, const smspec *spec);
cpl_error_code sm_scat_createtaugrid(cpl_table *taugrid,
                                     const smspec *molabstrans);
cpl_error_code sm_scat_calcscat(smscatres *scatres,
		                        const double sm_h, const double sm_hmin,
                                const double z0,
                                const double zeta0, const double theta0,
                                const double cextr, const double cextm,
                                const double cexta, const cpl_table *grid,
                                const smgrid *miephase,
                                const double lam0, const int mode);
cpl_error_code sm_scat_calcdoublescat(smscatres *scatres,
                                      const int ns,
                                      const double z, const double zeta,
                                      const double alphaa0, const double hx,
                                      const double bxr, const double bxm,
                                      const double b1r, const double b1m,
                                      const double cextr, const double cextm,
                                      const double cexta,
                                      const cpl_table *mgrid,
                                      const smgrid *miephase,
									  const double lam0, const double sm_hmin);
cpl_error_code sm_scat_initscatres(smscatres *scatres);
cpl_error_code sm_scat_procscatres(smscatres *scatres, const double flux0,
                                   const double lam0, const double tau0,
                                   const double cextr, const double cextm,
                                   const double ssa);
cpl_error_code sm_scat_initscatgrid(cpl_table *grid, const int nz);
cpl_error_code sm_scat_modscatgrid(cpl_table *grid, const double z0,
                                   const double a0);
double sm_scat_geteffcoldens(const double z, const double sig, const int ns,
                             const char scattype, const double sm_hmin);
double sm_scat_calceffcoldens(const double z, const double sig, const int ns,
                              const char scattype);
double sm_scat_esteffcoldens(const double z, const double hmin,
                             const double hmax, const char scattype,
							 const double sm_hmin);
double sm_scat_getmolecdens(const double h);
double sm_scat_getaerosoldens(const double h);
double sm_scat_getgrefl(const double lam);
cpl_error_code sm_scat_corrabs(smspec *spec, const cpl_table *xefftab,
                               const smspec *molabstrans);
cpl_error_code sm_scat_estmultiscat(const char *sscatcorfile,
                                    const char *dscatcorfile,
                                    const smparmodel modelpar);
cpl_error_code sm_scat_createztab(cpl_table *ztab);
cpl_error_code sm_scat_createzmtab(cpl_table *zmtab);
cpl_error_code sm_scat_createrhotab(cpl_table *rhotab);
cpl_error_code sm_scat_createdatatab(cpl_table *datatab,
                                     const cpl_table *ztab,
                                     const cpl_table *zmtab,
                                     const cpl_table *rhotab,
                                     const smspec *lamgrid);
cpl_error_code sm_scat_copyscatres(cpl_table *datatab, const int row,
                                   const smscatres *scatres);
cpl_error_code sm_scat_createzavtab(cpl_table *zavtab,
                                    const cpl_table *rhotab,
                                    const smspec *lamgrid);
cpl_error_code sm_scat_calczavdoublescat(cpl_table *zavtab,
                                         cpl_table *datatab,
                                         const cpl_table *ztab,
                                         const cpl_table *zmtab);
cpl_error_code sm_scat_createrhoavtab(cpl_table *rhoavtab,
                                      const smspec *lamgrid);
cpl_error_code sm_scat_calcrhoavdoublescat(cpl_table *rhoavtab,
                                           cpl_table *zavtab,
                                           const cpl_table *rhotab);
cpl_error_code sm_scat_calcscatcorr(cpl_table *zavtab, cpl_table *rhoavtab);
cpl_error_code sm_scat_writesinglescatcorr(const char *sscatcorfile,
                                           const cpl_table *zavtab,
                                           const cpl_table *rhotab,
                                           const smspec *lamgrid);
cpl_error_code sm_scat_writedoublescatcorr(const char *dscatcorfile,
                                           const cpl_table *rhoavtab);

#endif /* SCATMOONLIGHT_H */

#ifdef __cplusplus
}
#endif

/**@}*/
