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
 * \callgraph
 *
 * \file sm_scatmoonlight.c
 *
 * Routines for calculation of scattered moonlight
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  05 Jun 2012
 * \date   06 Oct 2015
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <sm_scatmoonlight.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

cpl_error_code sm_scat_moon(smspec *spec, 
                            const double sm_h, const double sm_hmin,
                            const double z, const double zmoon,
                            const double rho, const double ssa,
                            const smspec *molabstrans,
                            const smspec *rscattrans,
                            const smspec *mscattrans, const smgrid *miephase,
                            const smgrid *sscatcor, const smspec *dscatcor,
                            const char *calcds)
{
    /*!
     * \callgraph
     * Atmospheric scattering calculations for moonlight. While single
     * scattering calculations are always performed, double scattering is
     * only calculated if requested. The latter considers mixed scattering at
     * air molecules, aerosols, and the surface. Required input parameters are
     * line-of-sight zenith distance, zenith distance of Moon, angle between
     * line of sight and Moon, and single scattering albedo for aerosols.
     * Except for the latter, all parameters have to be provided in deg. Be
     * aware that there are impossible combinations of the three angles, which
     * cause a programme termination without scattering results. In the case
     * of a successful completion, the procedure returns a spectrum of
     * scattering intensities per \f${\rm arcsec}^{2}\f$ for an input total
     * Moon flux of 1. The required conversion from wavelength to \f$\tau\f$
     * is carried out by means of the input zenithal transmission curves for
     * Rayleigh and Mie scattering. For estimating the contribution of
     * multiple scattering, corrections are independently derived from single
     * and double scattering calculations for a large set of realistic
     * parameter sets (see ::sm_scat_estmultiscat and subroutines). They are
     * provided as data tables to this routine. Before multiplication with the
     * wavelength-dependent scattering intensities, the transmission curve for
     * molecular absorption is adapted to the effective airmasses related to
     * the scattering paths of both scattering processes. This procedure works
     * best for Rayleigh scattering and absorption by \f${\rm O}_2\f$.
     *
     * \b INPUT:
     * \param spec         desired wavelength grid
     * \param sm_h         observatory height in km 
     * \param sm_hmin      Lower height limit in km
     * \param z            line-of-sight zenith distance in deg
     * \param zmoon        zenith distance of Moon in deg
     * \param rho          angle between line of sight and Moon in deg
     * \param ssa          single scattering albedo for aerosols
     * \param molabstrans  transmission curve from radiative transfer code
     * \param rscattrans   transmission curve for Rayleigh scattering
     * \param mscattrans   transmission curve for Mie scattering
     * \param miephase     grid structure for Mie scattering phase functions
     * \param sscatcor     grid structure for multiple scattering correction
     *                     factors for single scattering
     * \param dscatcor     spectrum of multiple scattering intensities to
     *                     correct double scattering
     * \param multiscat    grid structure for multiple scattering correction
     *                     factors
     * \param calcds       calculation of double scattering ('Y' or 'N')
     *
     * \b OUTPUT:
     * \param spec         spectrum of scattered moonlight (input flux = 1)
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     * - IDG: Inconsistent data grids
     */

    cpl_table *grid, *xefftab;
    smspec iscatspec;
    smscatres scatres;
    char errtxt[SM_MAXLEN+1];
    int mode = SM_SCATMODE; // calculation mode (see header)
    int nz = 0, i = 0, n = 0, j = 0;
    double zrad = 0., zmrad = 0., rhorad = 0., sig0 = 0., b0r = 0., b0m = 0.;
    double lam0 = 0., rtrans = 0., mtrans = 0., tau0r = 0., tau0m = 0.;
    double tau0 = 0., cextr = 0., cextm = 0., cexta = 0., flux0 = 1.;
    double ims = 0., iscat = 0., mscorr = 0., reftau0r = 0., reftau0m = 0.;
    double reftau0 = 0., tau0a = 0., tau0sum = 0.;

    /* Check for impossible parameter combinations */
    if (fabs(z - zmoon) - SM_TOL > rho || z + zmoon + SM_TOL < rho) {
        sprintf(errtxt, "%s: invalid combination of z, zmoon, and rho",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Conversion from deg to rad */
    zrad = CPL_MATH_RAD_DEG * z;
    zmrad = CPL_MATH_RAD_DEG * zmoon;
    rhorad = CPL_MATH_RAD_DEG * rho;

    /* Test agreement of wavelength grids in input spectra */
    if ((int) sm_spec_compgrids(molabstrans, spec) == (int) SM_ERROR_IDG ||
        (int) sm_spec_compgrids(rscattrans, spec) == (int) SM_ERROR_IDG ||
        (int) sm_spec_compgrids(mscattrans, spec) == (int) SM_ERROR_IDG) {
        sprintf(errtxt, "%s: *molabstrans, *rscattrans, or *mscattrans != "
                "smspec *spec", SM_ERROR_IDG_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IDG, "%s", errtxt);
    }

    /* Height of observer in km relative to centre of Earth */
    sig0 = SM_R + sm_h;

    /* Get effective column density in cm^-2 for zenith */
    b0r = sm_scat_geteffcoldens(0., sig0, SM_NS, 'r', sm_hmin);
    b0m = sm_scat_geteffcoldens(0., sig0, SM_NS, 'm', sm_hmin);

    /* Create coordinate grid for double scattering calculations if
       requested */
    grid = cpl_table_new(0);
    if (calcds[0] == 'Y' || calcds[0] == 'y') {
        if (mode == 0) {
            nz = SM_NZ;
        } else {
            nz = SM_NZD;
        }
        sm_scat_initscatgrid(grid, nz);
    }

    /* Create smspec structure for scattering intensities */
    sm_scat_createlambdagrid(&iscatspec, spec);

    /* Calculations for scattering */

    for (i = 0; i < iscatspec.n; i++) {

        /* Get transmission values by Rayleigh and Mie scattering for listed
           wavelength */
        lam0 = iscatspec.dat[i].lam;
        sm_spec_getval(&rtrans, rscattrans, lam0);
        sm_spec_getval(&mtrans, mscattrans, lam0);

        /* Calculate optical depths for listed wavelength */
        tau0r = -log(rtrans);
        tau0m = -log(mtrans);
        tau0 = tau0r + tau0m;

        /* Correct for numerical uncertainties */
        if (lam0 < miephase->xpos[0]) {
            lam0 = miephase->xpos[0];
        }
        if (lam0 > miephase->xpos[(miephase->nx)-1]) {
            lam0 = miephase->xpos[(miephase->nx)-1];
        }

        /* Calculate Rayleigh and aerosol extinction cross section in cm^2 */
        cextr = tau0r / b0r;
        cextm = tau0m / b0m;

        /* No molecular absorption */
        cexta = 0;

        /* Default values of output parameters */
        sm_scat_initscatres(&scatres);

        /* Calculate intensity of radiation scattered into line of sight */
        sm_scat_calcscat(&scatres,
        		 sm_h, sm_hmin,
        		 zrad, zmrad, rhorad, cextr, cextm, cexta,
                         grid, miephase, lam0, mode);

        /* Process results of scattering calculations (retrieve final
           intensities and the effective airmass) */
        sm_scat_procscatres(&scatres, flux0, lam0, tau0, cextr, cextm, ssa);

        /* Correction for multiple scattering based on single + double
           scattering calculations */
        if (calcds[0] == 'Y' || calcds[0] == 'y') {
            /* Double scattering calculations: add scattering intensity for
               3 and more scattering events */
            sm_spec_getval(&ims, dscatcor, lam0);
            iscat = scatres.iscat + ims;
        } else {
            /* Single scattering calculations: correction for 2 and more
               scattering events */
            sm_grid_extract(&mscorr, sscatcor, lam0, rho);
            //sm_grid_extract(&mscorr, sscatcor, tau0, rho);
            iscat = scatres.iscat1s * mscorr;
        }

        /* Write scattering intensity into spectrum */
        iscatspec.dat[i].flux = iscat;

    }

    /* Interpolate scattering intensities to match wavelength grid of output
       spectrum */
    sm_spec_interpol(spec, &iscatspec);

    /* Free memory */
    sm_spec_free(&iscatspec);

    /* Create CPL table for effective airmasses */
    xefftab = cpl_table_new(0);
    sm_scat_createtaugrid(xefftab, molabstrans);
    cpl_table_new_column(xefftab, "X_eff", CPL_TYPE_DOUBLE);
    n = cpl_table_get_nrow(xefftab);
    cpl_table_fill_column_window(xefftab, "X_eff", 0, n, 0.);

    /* Get transmission values by Rayleigh and Mie scattering for reference
       wavelength */
    sm_spec_getval(&rtrans, rscattrans, SM_REFLAM0);
    sm_spec_getval(&mtrans, mscattrans, SM_REFLAM0);

    /* Calculate optical depths for reference wavelength */
    reftau0r = -log(rtrans);
    reftau0m = -log(mtrans);
    reftau0 = reftau0r + reftau0m;

    /* Calculate Rayleigh and aerosol extinction cross section in cm^2 for
       reference wavelength */
    cextr = reftau0r / b0r;
    cextm = reftau0m / b0m;

    /* Calculations for absorption */

    for (j = 0; j < n; j++) {

        /* Get absorption optical depth from CPL table */
        tau0a = cpl_table_get(xefftab, "tau0", j, NULL);

        /* Calculate molecular absorption cross section in cm^2 assuming same
           particle distribution as for Rayleigh scattering */
        cexta = tau0a / b0r;

        /* Calculate total optical depth */
        tau0sum = reftau0 + tau0a;

        /* Default values of output parameters */
        sm_scat_initscatres(&scatres);

        /* Calculate effective airmass for different absorption optical
           depths */
        sm_scat_calcscat(&scatres,
                         sm_h, sm_hmin,
                         zrad, zmrad, rhorad, cextr, cextm, cexta,
                         grid, miephase, SM_REFLAM0, mode);

        /* Process results of scattering calculations (retrieve final
           intensities and the effective airmass) */
        sm_scat_procscatres(&scatres, flux0, lam0, tau0sum, cextr, cextm,
                            ssa);

        /* Write effective airmass to CPL table */
        cpl_table_set(xefftab, "X_eff", j, scatres.xeff);

    }

    /* Correct scattering intensities for molecular absorption */
    sm_scat_corrabs(spec, xefftab, molabstrans);

    /* Free memory */
    cpl_table_delete(grid);
    cpl_table_delete(xefftab);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_scat_createlambdagrid(smspec *lamgrid, const smspec *spec)
{
    /*!
     * Creates grid of wavelengths for scattering calculations. The
     * wavelengths are fixed. However, only those are considered which are
     * relevant for the wavelength range of the given spectrum. If the input
     * spectrum is NULL, no wavelengths are excluded.
     *
     * \b INPUT:
     * \param spec     wavelength grid of output spectrum
     *
     * \b OUTPUT:
     * \param lamgrid  wavelength grid for scattering calculations
     *
     * \b ERRORS:
     * - ISM: Insufficient memory
     */

    int i = 0, imin = 0, imax = 0, nsel = 0;
    double minlam = 0., maxlam = 0.;

    /* List of wavelengths */
    int nlam = 40;
    double lam[] = { 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38,
                     0.39, 0.40, 0.41, 0.42, 0.43, 0.45, 0.47, 0.50, 0.55,
                     0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00,
                     1.05, 1.10, 1.15, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70,
                     1.80, 2.00, 2.50,30.00};

    /* Get wavelength range from spectrum if not NULL*/
    if (spec == NULL) {
        minlam = 0.;
        maxlam = HUGE_VAL;
    } else {
        minlam = spec->dat[0].lam;
        maxlam = spec->dat[(spec->n)-1].lam;
    }

    /* Find first relevant wavelength from list */
    for (imin = nlam-1, i = 0; i < nlam; i++) {
        if (lam[i] > minlam) {
            if (i == 0) {
                imin = i;
            } else {
                imin = i-1;
            }
            break;
        }
    }

    /* Find last relevant wavelength from list */
    for (imax = 0, i = nlam-1; i >= 0; i--) {
        if (lam[i] < maxlam) {
            if (i == nlam-1) {
                imax = i;
            } else {
                imax = i+1;
            }
            break;
        }
    }

    /* Get number of wavelengths in output grid */
    nsel = imax - imin + 1;

    /* Allocate memory (set spec->type = 3) */
    if ((int) sm_spec_malloc(lamgrid, nsel) == (int) SM_ERROR_ISM) {
        return SM_ERROR_ISM;
    }

    /* Reset spec->type = 2 */
    lamgrid->type = 2;

    /* Transfer wavelengths */
    for (i = imin; i <= imax; i++) {
        lamgrid->dat[i-imin].lam = lam[i];
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_scat_createtaugrid(cpl_table *taugrid,
                                     const smspec *molabstrans)
{
    /*!
     * Writes CPL table with grid of optical depths for scattering
     * calculations with absorption. The optical depths are fixed. However,
     * only those are considered which are relevant for the input transmission
     * curve for molecular absorption.
     *
     * \b INPUT:
     * \param taugrid      empty CPL table
     * \param molabstrans  transmission curve from radiative transfer code
     *
     * \b OUTPUT:
     * \param taugrid      optical depth grid for scattering calculations
     *
     * \b ERRORS:
     * - none
     */

    int i = 0, imin = 0, imax = 0, nsel = 0;
    double limtrans[] = {0., 0.}, limlam[] = {0., HUGE_VAL};
    double mintau = 0., maxtau = 0.;

    /* List of optical depths */
    /*
    int ntau = 27;
    double tau[] = { 0.00, 0.01, 0.02, 0.03, 0.05, 0.10, 0.15, 0.20, 0.25,
                     0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.20,
                     1.40, 1.60, 1.80, 2.00, 2.50, 3.00, 3.50, 4.00, 5.00};
    */
    int ntau = 11;
    double tau[] = { 0.00, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.00, 2.00,
                     5.00,10.00};

    /* Create wavelength column in CPL table */
    cpl_table_new_column(taugrid, "tau0", CPL_TYPE_DOUBLE);

    /* Get transmission range from transmission curve */
    sm_spec_minmax(limtrans, molabstrans, limlam);

    /* Calculate minimum and maximum optical depth */
    mintau = -log(limtrans[1]);
    maxtau = -log(limtrans[0]);

    /* Find first relevant optical depth from list */
    for (imin = ntau-1, i = 0; i < ntau; i++) {
        if (tau[i] > mintau) {
            if (i == 0) {
                imin = i;
            } else {
                imin = i-1;
            }
            break;
        }
    }

    /* Find last relevant optical depth from list */
    for (imax = 0, i = ntau-1; i >= 0; i--) {
        if (tau[i] < maxtau) {
            if (i == ntau-1) {
                imax = i;
            } else {
                imax = i+1;
            }
            break;
        }
    }

    /* Set table size */
    nsel = imax - imin + 1;
    cpl_table_set_size(taugrid, nsel);

    /* Write wavelengths into CPL table */
    for (i = imin; i <= imax; i++) {
        cpl_table_set(taugrid, "tau0", i-imin, tau[i]);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_scat_calcscat(smscatres *scatres, 
                                const double sm_h,
                                const double sm_hmin,
                                const double z0,
                                const double zeta0, const double theta0,
                                const double cextr, const double cextm,
                                const double cexta, const cpl_table *grid,
                                const smgrid *miephase,const double lam0,
                                const int mode)
{
    /*!
     * \callgraph
     * Calculation of intensity of scattered moonlight, which is scattered
     * into the given line of sight through the Earth's atmosphere, and
     * effective airmass relative to a direct path from zenith to observer.
     * Rayleigh scattering, aerosol extinction, and molecular absorption are
     * considered. For the latter, the same particle distribution as for
     * Rayleigh scattering is assumed (optimal for molecular oxygen). The
     * wavelength-dependent phase functions for Mie scattering were calculated
     * based on optimised log-normal aerosol distributions (Warneck \&
     * Williams 2012) and the Bohren-Huffman Mie scattering algorithm (IDL
     * code BHMIE). The single scattering formalism in an spherical atmosphere
     * used is described by Wolstencroft \& van Breda (1967), Staude (1975),
     * and Bernstein et al. (2002). As in the latter reference, polarisation
     * is neglected. Unlike the references, the subroutine
     * ::sm_scat_calcdoublescat also estimates the effect of double scattering
     * and ground reflection if the input coordinate grid has not zero
     * entries. For retrieving the final scattering intensities, the routine
     * ::sm_scat_procscatres has to be started afterwards.
     *
     * \b INPUT:
     * \param sm_h      observatory height in km
     * \param sm_hmin   Lower height limit in km
     * \param z0        line-of-sight zenith distance in rad
     * \param zeta0     zenith distance of radiation source in rad
     * \param theta0    angle between line of sight and radiation source
     * \param cextr     Rayleigh extinction cross section in
     *                  \f${\rm cm}^{2}\f$
     * \param cextm     Mie extinction cross section in \f${\rm cm}^{2}\f$
     * \param cexta     molecular absorption cross section in
     *                  \f${\rm cm}^{2}\f$
     * \param grid      CPL table with coordinate grid
     * \param miephase  grid structure for Mie scattering phase functions
     * \param lam0      wavelength in \f$\mu{\rm m}\f$
     * \param mode      calculation mode (0 = most accurate, 3 = fastest)
     *
     * \b OUTPUT:
     * \param scatres   ::smscatres structure containing the intensities
     *                  scattered into line of sight for different scattering
     *                  processes and the effective airmass
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     */

    cpl_table *mgrid = NULL;
    char errtxt[SM_MAXLEN+1];
    int ns = 0, nsc = 0, nsd = 0, nrow = 0, i = 0;
    double sigmax = 0., sigmin = 0., sig0 = 0., sinz0 = 0., cosz0 = 0.;
    double smax0 = 0., dlns = 0., b0r = 0., b0m = 0., costheta0 = 0.;
    double p0r = 0., p0m = 0., sinzeta0 = 0., coszeta0 = 0., cosalpha0a0 = 0.;
    double alpha0a0 = 0., spmin = 0., s = 0., spmax = 0., ds = 0., sig = 0.;
    double sinpsi = 0., psi = 0., cospsi = 0., coszeta = 0., zeta = 0.;
    double z = 0., sinzeta = 0., cosalphaa0 = 0., alphaa0 = 0., zetamax = 0.;
    double nsigr = 0., nsigm = 0., b1r = 0., b1m = 0., b2r = 0., b2m = 0.;
    double bxr = 0., bxm = 0., taur = 0., taum = 0., tau = 0., trans = 0.;
    double res[] = {0., 0.};

    /* Calculation mode (higher numbers -> faster but less accurate) */
    if (mode <= 0) {
        ns = SM_NS;
        nsc = SM_NS;
        nsd = SM_NS;
    } else if (mode == 1) {
        ns = SM_NS;
        nsc = SM_NS;
        nsd = SM_NSD;
    } else if (mode == 2) {
        ns = SM_NS;
        nsc = -1;
        nsd = SM_NSD;
    } else if (mode >= 3) {
        ns = SM_NSD;
        nsc = -1;
        nsd = SM_NSD;
    }

    /* Top of atmosphere in km relative to centre of Earth */
    sigmax = SM_R + SM_HMAX;
    /* Lower height limit in km relative to centre of Earth */
    sigmin = SM_R + sm_hmin;
    /* Height of observer in km relative to centre of Earth */
    sig0 = SM_R + sm_h;
    /* Sine and cosine of line-of-sight zenith distance */
    sinz0 = sin(z0);
    cosz0 = cos(z0);
    /* Maximum distance of top of atmosphere from observer */
    smax0 = sqrt(sigmax * sigmax - sig0 * sig0 * sinz0 * sinz0)
            - sig0 * cosz0;
    /* Logarithmic path length bin size */
    dlns = log(smax0 / SM_LNSCALE) / (double) ns;

    /* Effective column density in cm^-2 for line-of-sight zenith distance */
    b0r = sm_scat_geteffcoldens(z0, sig0, nsc, 'r', sm_hmin);
    b0m = sm_scat_geteffcoldens(z0, sig0, nsc, 'm', sm_hmin);

    /* Cosine of scattering angle at X */
    costheta0 = cos(theta0);

    /* Get phase function value for scattering angle, wavelength (Mie only),
       and type of scattering */
    /* Theoretical Rayleigh phase function */
    p0r = 0.75 * (1 + costheta0 * costheta0);
    /* Mie phase function from grid structure */
    sm_grid_extract(&p0m, miephase, lam0, theta0 / CPL_MATH_RAD_DEG);

    /* Sine and cosine of angle between zenith and radiation source at
       position of observer */
    sinzeta0 = sin(zeta0);
    coszeta0 = cos(zeta0);

    /* Angle between line-of-sight and radiation source azimuth at position of
       observer */
    if (sinz0 == 0. || sinzeta0 == 0.) {
        cosalpha0a0 = 1.;
        /* Error if combination of input parameters is not possible */
        if (fabs(fabs(z0 - zeta0) - theta0) > SM_TOL) {
            sprintf(errtxt, "%s: invalid combination of z0, zeta, and theta",
                    SM_ERROR_IIP_TXT);
            scatres->iscat = -1;
            return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s",
                                         errtxt);
        }
    } else {
        cosalpha0a0 = (costheta0 - cosz0 * coszeta0) / (sinz0 * sinzeta0);
        if (cosalpha0a0 > -1 - SM_TOL && cosalpha0a0 < -1) {
            cosalpha0a0 = -1.;
        } else if (cosalpha0a0 > 1 && cosalpha0a0 < 1 + SM_TOL) {
            cosalpha0a0 = 1.;
        }
        /* Error if combination of input parameters is not possible */
        if (cosalpha0a0 < -1 || cosalpha0a0 > 1) {
            sprintf(errtxt, "%s: invalid combination of z0, zeta0, and "
                    "theta0", SM_ERROR_IIP_TXT);
            scatres->iscat = -1;
            return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s",
                                         errtxt);
        }
    }
    alpha0a0 = acos(cosalpha0a0);

    /* Refine coordinate grid for Mie forward scattering if double scattering
       is requested */
    mgrid = cpl_table_duplicate(grid);
    nrow = cpl_table_get_nrow(mgrid);
    if (nrow > 0) {
        sm_scat_modscatgrid(mgrid, zeta0, alpha0a0);
        sm_scat_modscatgrid(mgrid, z0, 0.);
    }

    /* Initialise lower limit of path element */
    spmin = 0.;

    /* Compute scattered intensity for different scattering positions in the
       atmosphere (representing volume elements) along the line of sight */

    for (i = -1; i < ns; i++) {

        /* Distance of scattering point X from observer in km (exponential
           increase of distances) */
        if (i == -1) {
            s = SM_LNSCALE / 2.;
        } else {
            s = SM_LNSCALE * exp(((double) i + 0.5) * dlns);
        }
        /* Upper limit of path element */
        spmax = SM_LNSCALE * exp((double) (i + 1) * dlns);
        /* Length of path element */
        ds = spmax - spmin;
        /* Upper limit becomes lower limit in next iteration */
        spmin = spmax;
        /* Height of scattering point in km relative to centre of Earth */
        sig = sqrt(sig0 * sig0 + s * s + 2 * sig0 * s * cosz0);

        /* Consider 3D effects? */
        if (mode == 0) {
            /* Angle between observer and scattering point at centre of
               Earth */
            sinpsi = s * sinz0 / sig;
            psi = asin(sinpsi);
            cospsi = sqrt(1 - sinpsi * sinpsi);
            /* Zenith distance of radiation source at scattering point X */
            coszeta = coszeta0 * cospsi + sinzeta0 * sinpsi * cosalpha0a0;
            zeta = acos(coszeta);
            /* Only transform the following parameters if double scattering is
               requested */
            if (nrow > 0) {
                /* Zenith distance of line of sight at scattering point X */
                z = z0 - psi;
                /* Angle between radiation source and line-of-sight azimuth at
                   scattering point X */
                sinzeta = sqrt(1 - coszeta * coszeta);
                if (sinzeta == 0. || sinpsi == 0.) {
                    cosalphaa0 = 1.;
                } else {
                    cosalphaa0 = (cospsi * coszeta - coszeta0) /
                                 (sinpsi * sinzeta);
                    if (cosalphaa0 < -1) {
                        cosalphaa0 = -1.;
                    } else if (cosalphaa0 > 1) {
                        cosalphaa0 = 1.;
                    }
                }
                alphaa0 = acos(cosalphaa0);
            } else {
                z = z0;
                alphaa0 = alpha0a0;
            }
        } else {
            psi = 0.;
            zeta = zeta0;
            z = z0;
            alphaa0 = alpha0a0;
        }

        /* Maximum zenith distance of radiation source for scattering point X
           (limited by lower height limit) */
        zetamax = CPL_MATH_PI_2 + acos(sigmin / sig);
        /* Check that selected zenith distance is below maximum value */
        if (zeta > zetamax) {
            /* Centre of radiation source below horizon -> Skip */
            continue;
        }

        /* Get particle density at scattering point in cm^-3 for Rayleigh and
           Mie (aerosol) scattering */
        /* Rayleigh scattering */
        nsigr = sm_scat_getmolecdens(sig - SM_R);
        /* Mie (aerosol) scattering */
        nsigm = sm_scat_getaerosoldens(sig - SM_R);

        /* Effective column density in cm^-2 for path from scattering point X
           to observer */
        b1r = b0r - sm_scat_geteffcoldens(z0 - psi, sig, nsc, 'r', sm_hmin);
        b1m = b0m - sm_scat_geteffcoldens(z0 - psi, sig, nsc, 'm', sm_hmin);
        /* Effective column density in cm^-2 for path from top of atmosphere
           to scattering point X */
        b2r = sm_scat_geteffcoldens(zeta, sig, nsc, 'r', sm_hmin);
        b2m = sm_scat_geteffcoldens(zeta, sig, nsc, 'm', sm_hmin);

        /* Column density in cm^-2 of path element at X */
        bxr = nsigr * 1e5 * ds;
        bxm = nsigm * 1e5 * ds;

        /* Transmission along full light path */
        taur = (cextr + cexta) * (b1r + b2r);
        taum = cextm * (b1m + b2m);
        tau = taur + taum;
        trans = exp(-tau);

        /* Unextinguished, scattered light (for effective airmass) */
        res[0] = bxr * p0r;
        scatres->it1r += res[0];
        res[1] = bxm * p0m;
        scatres->it1m += res[1];

        /* Scattered light for line-of-sight path element */
        scatres->iscatr += res[0] * trans;
        scatres->iscatm += res[1] * trans;

        /* Calculate double scattering if requested */
        if (nrow > 0) {
            sm_scat_calcdoublescat(scatres, nsd, z, zeta, alphaa0, sig - SM_R,
                                   bxr, bxm, b1r, b1m, cextr, cextm, cexta,
                                   mgrid, miephase, lam0, sm_hmin);
        }

    }

    /* Free memory */
    cpl_table_delete(mgrid);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_scat_calcdoublescat(smscatres *scatres, const int ns,
                                      const double z, const double zeta,
                                      const double alphaa0, const double hx,
                                      const double bxr, const double bxm,
                                      const double b1r, const double b1m,
                                      const double cextr, const double cextm,
                                      const double cexta,
                                      const cpl_table *mgrid,
                                      const smgrid *miephase,
                                      const double lam0, const double sm_hmin)
{
    /*!
     * \callgraph
     * Calculation of intensity of scattered moonlight by double scattering.
     * For a scattering point coming from the single scattering calculations
     * by ::sm_scat_calcscat, another scattering event is considered, which
     * can be caused by either a molecule, an aerosol particle, or the ground.
     * The routine returns for each possible combination of two scattering
     * events the corresponding scattering intensities (neglecting constant
     * factors as calculated by ::sm_scat_calcscat) as ::smscatres structure.
     * In order to reduce the run time, the angular positions of the Moon and
     * the line of sight and the effective column densities are approximated.
     *
     * \b INPUT:
     * \param scatres   ::smscatres structure containing the results of
     *                  previous scattering calculations
     * \param ns        reference number of path elements
     * \param z         line-of-sight zenith distance in rad at scattering
     *                  point X
     * \param zeta      zenith distance of radiation source at X in rad
     * \param alphaa0   azimuth angle of radiation source relative to line of
     *                  sight in rad at X
     * \param hx        height of X in km
     * \param bxr       molecular column density in \f${\rm cm}^{-2}\f$ of
     *                  path element at X
     * \param bxm       aerosol column density in \f${\rm cm}^{-2}\f$ of path
     *                  element at X
     * \param b1r       molecular column density in \f${\rm cm}^{-2}\f$ for
     *                  path from X to observer
     * \param b1m       aerosol column density in \f${\rm cm}^{-2}\f$ for path
     *                  from X to observer
     * \param cextr     Rayleigh extinction cross section in
     *                  \f${\rm cm}^{2}\f$
     * \param cextm     Mie extinction cross section in \f${\rm cm}^{2}\f$
     * \param cexta     molecular absorption cross section in
     *                  \f${\rm cm}^{2}\f$
     * \param mgrid     CPL table with modified coordinate grid
     * \param miephase  grid structure for Mie scattering phase functions
     * \param lam0      wavelength for Mie scattering in \f$\mu{\rm m}\f$
     *
     * \b OUTPUT:
     * \param scatres   ::smscatres structure containing the updated
     *                  intensities for double scattering
     *
     * \b ERRORS:
     * - none
     */

    cpl_boolean newzp = CPL_TRUE;
    int nrow = 0, j = 0, ghit = 0, nl = 0, k = 0;
    double sigmax = 0., sigmin = 0., sig = 0., sinz = 0., cosz = 0.;
    double sinzeta = 0., coszeta = 0., dlns0 = 0., zpo = 0., zp = 0.;
    double apa = 0., domega = 0., sinzp = 0., coszp = 0., cosapa = 0.;
    double alphaap = 0., cosalphaap = 0., costheta = 0., theta = 0., pr = 0.;
    double pm = 0., term = 0., lmax = 0., dlnl = 0., lpmin = 0., l = 0.;
    double lpmax = 0., dl = 0., sigp = 0., zetapmax = 0., phi = 0., zpp = 0.;
    double costhetap = 0., thetap = 0., ppr = 0., ppm = 0., nsigpr = 0.;
    double nsigpm = 0., b0pr = 0., b0pm = 0., b2pr = 0., b2pm = 0., b3pr = 0.;
    double b3pm = 0., taupr = 0., taupm = 0., taup = 0., transp0 = 0.;
    double res[] = {0., 0., 0., 0., 0., 0.};
    double bxpr[ns+1], bxpm[ns+1], transp[ns+1];

    /* Top of atmosphere in km relative to centre of Earth */
    sigmax = SM_R + SM_HMAX;
    /* Lower height limit in km relative to centre of Earth */
    sigmin = SM_R + sm_hmin;
    /* Height of scattering point in km relative to centre of Earth */
    sig = SM_R + hx;
    /* Sine and cosine of line-of-sight zenith distance at scatering point
       X */
    sinz = sin(z);
    cosz = cos(z);
    /* Sine and cosine of angle between zenith and radiation source at
       scattering point X */
    sinzeta = sin(zeta);
    coszeta = cos(zeta);
    /* Reference logarithmic path length bin size */
    dlns0 = log(SM_HMAX / SM_LNSCALE) / (double) ns;

    /* Get number of points of spherical coordinate grid */
    nrow = cpl_table_get_nrow(mgrid);

    /* Compute contribution to line-of-sight intensity by secondary scattering
       from all directions */

    for (zpo = -1, j = 0; j < nrow; j++) {

        /* Read data grid table: zenith distance of scattering point X' at X,
           azimuth of X' at X, and solid angle */
        zp = cpl_table_get(mgrid, "z", j, NULL);
        apa = cpl_table_get(mgrid, "A", j, NULL); // = apa0
        domega = cpl_table_get(mgrid, "domega", j, NULL);

        /* Check change of zenith distance of X' at X */
        if (fabs(zp - zpo) > SM_TOL) {
            newzp = CPL_TRUE;
        } else {
            newzp = CPL_FALSE;
        }
        zpo = zp;

        /* Sine and cosine for read parameters */
        if (newzp == CPL_TRUE) {
            sinzp = sin(zp);
            coszp = cos(zp);
        }
        cosapa = cos(apa);

        /* Azimuth difference between X' and radiation source at scattering
           point X */
        alphaap = alphaa0 - apa;
        cosalphaap = cos(alphaap);
        /* Scattering angle at scattering point X */
        costheta = cosz * coszp + sinz * sinzp * cosapa;
        if (costheta < -1.) {
            costheta = -1.;
        } else if (costheta > 1.) {
            costheta = 1.;
        }
        theta = acos(costheta);

        /* Phase function value for scattering angle at scattering point X for
           Rayleigh and Mie (aerosol) scattering */
        /* Rayleigh scattering */
        pr = 0.75 * (1 + costheta * costheta);
        /* Mie (aerosol) scattering */
        sm_grid_extract(&pm, miephase, lam0, theta / CPL_MATH_RAD_DEG);

        /* Scattering angle at scattering point X' */
        costhetap = coszp * coszeta + sinzp * sinzeta * cosalphaap;
        if (costhetap < -1.) {
            costhetap = -1.;
        } else if (costhetap > 1.) {
            costhetap = 1.;
        }
        thetap = acos(costhetap);

        /* Phase function value for scattering angle at scattering point X'
           for Rayleigh and Mie (aerosol) scattering */
        /* Rayleigh scattering */
        ppr = 0.75 * (1 + costhetap * costhetap);
        /* Mie (aerosol) scattering */
        sm_grid_extract(&ppm, miephase, lam0, thetap / CPL_MATH_RAD_DEG);

        /* Calculate the following parameters only if zp changes */

        if (newzp == CPL_TRUE) {

            /* Effective column density in cm^-2 for path from top of
               atmosphere to X along direction XX' */
            if (zp > CPL_MATH_PI - zp) {
                b0pr = sm_scat_geteffcoldens(CPL_MATH_PI - zp, sig, -1, 'r', sm_hmin);
                b0pm = sm_scat_geteffcoldens(CPL_MATH_PI - zp, sig, -1, 'm', sm_hmin);
            } else {
                b0pr = sm_scat_geteffcoldens(zp, sig, -1, 'r', sm_hmin);
                b0pm = sm_scat_geteffcoldens(zp, sig, -1, 'm', sm_hmin);
            }

            /* Distance of scattering point X from ground along path */
            term = sigmin * sigmin - sig * sig * sinzp * sinzp;
            if (zp <= CPL_MATH_PI_2 || term < 0) {
                /* No hit of Earth */
                lmax = -1;
                ghit = 0;
            } else {
                lmax = - sqrt(term) - sig * coszp;
                ghit = 1;
            }
            /* Maximum distance of top of atmosphere from scattering point
               X */
            if (ghit == 0) {
                lmax = sqrt(sigmax * sigmax - sig * sig * sinzp * sinzp)
                       - sig * coszp;
            }

            /* Reduce number of path elements for path down to the ground */
            if (ghit == 1) {
                nl = (int) ceil(log(lmax / SM_LNSCALE) / dlns0);
                if (nl > ns) {
                    nl = ns;
                }
            } else {
                nl = ns;
            }
            /* Logarithmic path length bin size */
            dlnl = log(lmax / SM_LNSCALE) / (double) nl;
            /* Initialise lower limit of path element */
            lpmin = 0.;

        }

        /* Compute scattered intensity for different scattering positions in
           the atmosphere (representing volume elements) along the line
           between scattering points X' and X */

        for (k = -1; k < nl; k++) {

            /* Calculate the following parameters only if zp changes.
               Otherwise take already calculated bxpr, bxpm, and transp. */

            if (newzp == CPL_TRUE) {

                /* Initialise array values */
                bxpr[k+1] = -1;
                bxpm[k+1] = -1;
                transp[k+1] = -1;

                /* Skip if scattering point is below minimum height */
                if (lmax >= 0 && lpmin >= lmax) {
                    break;
                }

                /* Distance of X' from X in km (exponential increase of
                   distances) */
                if (k == -1) {
                    l = SM_LNSCALE / 2.;
                } else {
                    l = SM_LNSCALE * exp(((double) k + 0.5) * dlnl);
                }
                /* Upper limit of path element */
                lpmax = SM_LNSCALE * exp((double) (k + 1) * dlnl);
                /* Cut path element if ground is intersected */
                if (lmax >= 0 && lpmax > lmax) {
                    lpmax = lmax;
                }
                if (lmax >= 0 && l > lmax) {
                    l = lmax;
                }
                /* Length of path element */
                dl = lpmax - lpmin;
                /* Upper limit becomes lower limit in next iteration */
                lpmin = lpmax;

                /* Height of scattering point X' in km relative to centre of
                   Earth */
                sigp = sqrt(sig * sig + l * l + 2 * sig * l * coszp);

                /* Maximum zenith distance of radiation source for scattering
                   point X' (limited by lower height limit) */
                zetapmax = CPL_MATH_PI_2 + acos(sigmin / sigp);
                /* Check that radiation source is visible at X' */
                if (zeta > zetapmax) {
                    /* Centre of radiation source below horizon -> Skip */
                    continue;
                }

                /* Angle between scattering points X' and X at centre of Earth
                   (approximation for small angles) */
                phi = l * sinzp / sigp;
                /* Zenith distance z' at scattering point X' */
                zpp = zp - phi;

                /* Particle density at scattering point X' in cm^-3 for
                   Rayleigh and Mie (aerosol) scattering */
                /* Rayleigh scattering */
                nsigpr = sm_scat_getmolecdens(sigp - SM_R);
                /* Mie (aerosol) scattering */
                nsigpm = sm_scat_getaerosoldens(sigp - SM_R);

                /* Effective column density in cm^-2 for path from X' to X */
                if (zpp > CPL_MATH_PI - zpp) {
                    b2pr = fabs(b0pr - sm_scat_geteffcoldens(CPL_MATH_PI -zpp,
                                                             sigp, -1, 'r', sm_hmin));
                    b2pm = fabs(b0pm - sm_scat_geteffcoldens(CPL_MATH_PI -zpp,
                                                             sigp, -1, 'm', sm_hmin));
                } else {
                    b2pr = fabs(b0pr - sm_scat_geteffcoldens(zpp, sigp, -1,
                                                             'r', sm_hmin));
                    b2pm = fabs(b0pm - sm_scat_geteffcoldens(zpp, sigp, -1,
                                                             'm', sm_hmin));
                }
                /* Effective column density in cm^-2 for path from top of
                   atmosphere to scattering point X' */
                b3pr = sm_scat_geteffcoldens(zeta, sigp, -1, 'r', sm_hmin);
                b3pm = sm_scat_geteffcoldens(zeta, sigp, -1, 'm', sm_hmin);

                /* Column density in cm^-2 of path element at X' */
                bxpr[k+1] = nsigpr * 1e5 * dl;
                bxpm[k+1] = nsigpm * 1e5 * dl;

                /* Transmission along full light path */
                taupr = (cextr + cexta) * (b1r + b2pr + b3pr);
                taupm = cextm * (b1m + b2pm + b3pm);
                taup = taupr + taupm;
                transp[k+1] = exp(-taup);

            }

            /* Skip impossible scattering point (no calculation of
               transmission) */
            if (transp[k+1] == -1) {
                continue;
            }

            /* Unextinguished, scattered light (for effective airmass) */
            res[0] = bxr * pr * bxpr[k+1] * ppr * domega;
            scatres->it1rr += res[0];
            res[1] = bxr * pr * bxpm[k+1] * ppm * domega;
            scatres->it1rm += res[1];
            res[2] = bxm * pm * bxpr[k+1] * ppr * domega;
            scatres->it1mr += res[2];
            res[3] = bxm * pm * bxpm[k+1] * ppm * domega;
            scatres->it1mm += res[3];

            /* Light by scattering at X' and X */
            scatres->iscatrr += res[0] * transp[k+1];
            scatres->iscatrm += res[1] * transp[k+1];
            scatres->iscatmr += res[2] * transp[k+1];
            scatres->iscatmm += res[3] * transp[k+1];

        }

        /* Compute scattered intensity for possible scattering at the ground
           G */

        if (ghit == 1 && zeta <= CPL_MATH_PI_2) {

            /* Calculate the following parameters only if zp changes.
               Otherwise take already calculated transmission. */

            if (newzp == CPL_TRUE) {

                /* Effective column density in cm^-2 for path from X' to X */
                b2pr = fabs(b0pr - sm_scat_geteffcoldens(CPL_MATH_PI - zp,
                                                         sigmin, -1, 'r', sm_hmin));
                b2pm = fabs(b0pm - sm_scat_geteffcoldens(CPL_MATH_PI - zp,
                                                         sigmin, -1, 'm', sm_hmin));
                /* Effective column density in cm^-2 for path from top of
                   atmosphere to scattering point G */
                b3pr = sm_scat_geteffcoldens(zeta, sigmin, -1, 'r', sm_hmin);
                b3pm = sm_scat_geteffcoldens(zeta, sigmin, -1, 'm', sm_hmin);

                /* Transmission along full light path */
                taupr = (cextr + cexta) * (b1r + b2pr + b3pr);
                taupm = cextm * (b1m + b2pm + b3pm);
                taup = taupr + taupm;
                transp0 = exp(-taup);

            }

            /* Unextinguished, scattered light (for effective airmass) */
            res[4] = bxr * pr * domega;
            scatres->it1rg += res[4];
            res[5] = bxm * pm * domega;
            scatres->it1mg += res[5];

            /* Light by scattering at G and X */
            scatres->iscatrg += res[4] * transp0;
            scatres->iscatmg += res[5] * transp0;

        }

    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_scat_initscatres(smscatres *scatres)
{
    /*!
     * Initialises ::smscatres structure for results of scattering
     * calculations
     *
     * \b INPUT:
     * \param scatres   empty ::smscatres structure
     *
     * \b OUPUT:
     * \param scatres   initialised ::smscatres structure
     *
     * \b ERRORS:
     * - none
     */

    /* Scattering results considering extinction */
    scatres->iscatr = 0.;
    scatres->iscatm = 0.;
    scatres->iscatrr = 0.;
    scatres->iscatrm = 0.;
    scatres->iscatmr = 0.;
    scatres->iscatmm = 0.;
    scatres->iscatrg = 0.;
    scatres->iscatmg = 0.;
    scatres->iscat1s = 0.;
    scatres->iscat2s = 0.;
    scatres->iscat = 0.;

    /* Scattering results for no extinction */
    scatres->it1r = 0.;
    scatres->it1m = 0.;
    scatres->it1rr = 0.;
    scatres->it1rm = 0.;
    scatres->it1mr = 0.;
    scatres->it1mm = 0.;
    scatres->it1rg = 0.;
    scatres->it1mg = 0.;
    scatres->it1_1s = 0.;
    scatres->it1_2s = 0.;
    scatres->it1 = 0.;

    /* Effective airmass */
    scatres->xeff = 100.;

    return CPL_ERROR_NONE;
}


cpl_error_code sm_scat_procscatres(smscatres *scatres, const double flux0,
                                   const double lam0, const double tau0,
                                   const double cextr, const double cextm,
                                   const double ssa)
{
    /*!
     * Processes results for scattering intensities provided by a ::smscatres
     * structure. The individual scattering intensities are derived by
     * multiplying constant factors. Moreover, the total scattering intensity
     * is calculated. Finally, the effective airmass for the correction of
     * molecular absorption is estimated. If the input total intensity is not
     * zero, the routine terminates without changes in the ::smscatres
     * structure.
     *
     * \b INPUT:
     * \param scatres  ::smscatres structure with results of scattering
     *                 calculations
     * \param flux0    flux factor to be multiplied by the scattering
     *                 intensities
     * \param lam0     wavelength in \f$\mu{\rm m}\f$
     * \param tau0     optical depth at zenith
     * \param cextr    Rayleigh extinction cross section in \f${\rm cm}^{2}\f$
     * \param cextm    Mie extinction cross section in \f${\rm cm}^{2}\f$
     * \param ssa      single scattering albedo for aerosols
     *
     * \b OUPUT:
     * \param scatres  ::smscatres structure with final scattering intensities
     *
     * \b ERRORS:
     * - none
     */

    double i0 = 0., grefl = 0., fac_1s = 0., fac_2s = 0., fac_sg = 0.;
    double cscatr = 0., cscatm = 0.;

    /* Return without changes if the scattering results has already been
       processed */
    if (scatres->iscat != 0.) {
        return CPL_ERROR_NONE;
    }

    /* Conversion from total flux to flux per arcsec^2 */
    i0 = flux0 / (4 * CPL_MATH_PI * SM_SR_IN_ARCSEC2);

    /* Get ground reflectance for given wavelength */
    grefl = sm_scat_getgrefl(lam0);

    /* Get general factors for conversion into intensities */
    fac_1s = i0;
    fac_2s = i0 / (4 * CPL_MATH_PI);
    fac_sg = grefl * i0 / (2 * CPL_MATH_PI);

    /* Get scattering cross section in cm^2 depending on extinction cross
       section and type of scattering (Rayleigh or Mie) */
    /* No absorption for Rayleigh scattering */
    cscatr = cextr;
    /* Consider single scattering albedo for Mie scattering by aerosols
       (assumption: no wavelength dependence) */
    cscatm = ssa * cextm;

    /* Multiplication of general factors for final scattered light
       contribution to light of sight */

    /* Scattering results considering extinction */

    /* Single scattering */
    scatres->iscatr *= cscatr * fac_1s;
    scatres->iscatm *= cscatm * fac_1s;
    /* Double scattering (only one particle type) */
    scatres->iscatrr *= cscatr * cscatr * fac_2s;
    scatres->iscatrm *= cscatr * cscatm * fac_2s;
    scatres->iscatmr *= cscatm * cscatr * fac_2s;
    scatres->iscatmm *= cscatm * cscatm * fac_2s;
    /* Scattering at ground */
    scatres->iscatrg *= cscatr * fac_sg;
    scatres->iscatmg *= cscatm * fac_sg;
    /* Summed intensities for single and double scattering */
    scatres->iscat1s = scatres->iscatr + scatres->iscatm;
    scatres->iscat2s = scatres->iscatrr + scatres->iscatrm +
                       scatres->iscatmr + scatres->iscatmm +
                       scatres->iscatrg + scatres->iscatmg;
    scatres->iscat = scatres->iscat1s + scatres->iscat2s;

    /* Scattering results for no extinction */

    /* Single scattering */
    scatres->it1r *= cscatr * fac_1s;
    scatres->it1m *= cscatm * fac_1s;
    /* Double scattering (only one particle type) */
    scatres->it1rr *= cscatr * cscatr * fac_2s;
    scatres->it1rm *= cscatr * cscatm * fac_2s;
    scatres->it1mr *= cscatm * cscatr * fac_2s;
    scatres->it1mm *= cscatm * cscatm * fac_2s;
    /* Scattering at ground */
    scatres->it1rg *= cscatr * fac_sg;
    scatres->it1mg *= cscatm * fac_sg;
    /* Summed intensities for single and double scattering */
    scatres->it1_1s = scatres->it1r + scatres->it1m;
    scatres->it1_2s = scatres->it1rr + scatres->it1rm +
                      scatres->it1mr + scatres->it1mm +
                      scatres->it1rg + scatres->it1mg;
    scatres->it1 = scatres->it1_1s + scatres->it1_2s;

    /* Estimate effective airmass */
    if (scatres->it1 > 0 && scatres->iscat > 0) {
        scatres->xeff = log(scatres->it1 / scatres->iscat) / tau0;
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_scat_initscatgrid(cpl_table *grid, const int nz)
{
    /*!
     * Initialises coordinate grid for scattering. Zenith distance (rad),
     * azimuth (rad), and solid angle (sr) are written to a CPL table. Zenith
     * and nadir are included as grid points. Around the equator, the zenith
     * distance bins are made more narrow.
     *
     * \b INPUT:
     * \param grid  empty CPL table
     * \param nz    number of zenith distance values
     *
     * \b OUTPUT:
     * \param grid  CPL table with coordinate grid
     *
     * \b ERRORS:
     * - Invalid input parameters
     */

    char errtxt[SM_MAXLEN+1];
    int nzp = 0, nzpp = 0, iz1 = 0, iz2 = 0, i = 0, na = 0, n = 0, j = 0;
    double dz = 0., zp = 0., da = 0., zmin = 0., zmax = 0., dzl = 0.;
    double dzu = 0., domega = 0.;
    double z[nz], zg[3] = {0.6, 0.3, 0.1};

    /* Create table columns */
    cpl_table_set_size(grid, 0);
    cpl_table_new_column(grid, "z", CPL_TYPE_DOUBLE);
    cpl_table_new_column(grid, "A", CPL_TYPE_DOUBLE);
    cpl_table_new_column(grid, "dzl", CPL_TYPE_DOUBLE);
    cpl_table_new_column(grid, "dzu", CPL_TYPE_DOUBLE);
    cpl_table_new_column(grid, "dA", CPL_TYPE_DOUBLE);
    cpl_table_new_column(grid, "domega", CPL_TYPE_DOUBLE);

    /* Check nunber of zenith distances */
    if (nz < 6) {
        sprintf(errtxt, "%s: nz < 6", SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Get number of zenith distances in main grid (narrow bins excluded) */
    nzp = nz - 4;

    /* Get maximum difference between adjacent zenith distances and
       index numbers related to bins with smaller bin sizes */
    if (pow(-1, nzp) == -1) {
        nzpp = nzp - 1;
        iz1 = (nzp - 1) / 2;
        iz2 = iz1 + 3;
    } else {
        nzpp = nzp;
        iz1 = nzp / 2;
        iz2 = iz1 + 2;
    }

    /* Get maximum difference between adjacent zenith distances */
    dz =  CPL_MATH_PI / (double) nzpp;

    /* Set zenith distances of coordinate grid */

    for (i = 0; i < nz; i++) {

        if (i == 0) {
            zp = 0;
        } else if (i == nz - 1) {
            zp = CPL_MATH_PI;
        } else if (i == iz1 || i == iz2 + 2) {
            zp += dz * zg[0];
        } else if (i == iz1 + 1 || i == iz2 + 1) {
            zp += dz * zg[1];
        } else if (i == iz1 + 2 && iz1 + 3 == iz2) {
            zp += dz * zg[2];
        } else if (i == iz2 && iz1 + 3 == iz2) {
            zp += dz * zg[2];
        } else if (i == iz2 && iz1 + 2 == iz2) {
            zp += dz * zg[2] * 2;
        } else {
            zp += dz;
        }

        z[i] = zp;

    }

    /* Set azimuths and solid angles of coordinate grid */

    for (i = 0; i < nz; i++) {

        /* Get number of azimuths and corresponding bin size for given
           zenith distance */
        if (i == 0 || i == nz - 1) {
            na = 1;
        } else {
            na = (int) (4 * floor(2. * nzpp * sin(z[i]) / 4. + 0.5));
            if (na == 0) na = 4;
        }
        da = CPL_MATH_2PI / (double) na;

        /* Add rows to table */
        cpl_table_set_size(grid, n + na);

        /* Get minimum and maximum zenith distances */
        if (i == 0) {
            zmin = z[i];
            zmax = (z[i] + z[i+1]) / 2;
        } else if (i == nz - 1) {
            zmin = (z[i-1] + z[i]) / 2;
            zmax = z[i];
        } else {
            zmin = (z[i-1] + z[i]) / 2;
            zmax = (z[i] + z[i+1]) / 2;
        }

        /* Get zenith distance extension of cell */
        dzl = z[i] - zmin;
        dzu = zmax - z[i];

        /* Get solid angle for each grid point */
        domega = da * fabs(cos(zmin) - cos(zmax));

        /* Write coordinates to CPL table */
        for (j = 0; j < na; j++) {
            cpl_table_set(grid, "z", n + j, z[i]);
            cpl_table_set(grid, "A", n + j, j * da);
            cpl_table_set(grid, "dzl", n + j, dzl);
            cpl_table_set(grid, "dzu", n + j, dzu);
            cpl_table_set(grid, "dA", n + j, da);
            cpl_table_set(grid, "domega", n + j, domega);
        }

        /* Count number of grid points (= number of table rows) */
        n += na;

    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_scat_modscatgrid(cpl_table *grid, const double z0,
                                   const double a0)
{
    /*!
     * Produces a finer coordinate grid for scattering towards the direction
     * given by the input zenith distance (rad) and azimuth (rad). The
     * limiting distance from the requested position for the changes is given
     * by SM_DMAX. Each selected grid point is substituted by nine new points
     * representing a ninth of the initial solid angle.
     *
     * \b INPUT:
     * \param grid  CPL table with coordinate grid
     * \param z0    zenith distance in rad
     * \param a0    azimuth in rad
     *
     * \b OUTPUT:
     * \param grid  CPL table with modified coordinate grid
     *
     * \b ERRORS:
     * - none
     */

    int nrow = 0, i = 0, n = 0, k = 0, j = 0;
    double deltamax = 0., z = 0., a = 0., delta = 0., dzl = 0., dzu = 0.;
    double zmin = 0., dz = 0., da = 0., domega = 0., dzp = 0., dzhp = 0.;
    double zp = 0., dap = 0., ap = 0., domegap = 0., ap0 = 0.;

    /* Search radius in rad */
    deltamax = SM_DMAX * CPL_MATH_RAD_DEG;

    /* Get number of input data points */
    nrow = cpl_table_get_nrow(grid);

    /* Find cells close to given direction and divide them into smaller
       elements  */

    while (i < nrow) {

        /* Get grid coordinates */
        z = cpl_table_get(grid, "z", i, NULL);
        a = cpl_table_get(grid, "A", i, NULL);

        /* Calculate angular distance from requested direction */
        delta = acos(cos(z0) * cos(z) + sin(z0) * sin(z) * cos(a0 - a));

        /* Skip if the cell is outside the search radius */
        if (delta > deltamax) {
            i++;
            continue;
        }

        /* Get cell size */
        dzl = cpl_table_get(grid, "dzl", i, NULL);
        dzu = cpl_table_get(grid, "dzu", i, NULL);
        zmin = z - dzl;
        dz = dzu + dzl;
        da = cpl_table_get(grid, "dA", i, NULL);
        domega = cpl_table_get(grid, "domega", i, NULL);

        /* Set new cell sizes */
        if (z < SM_TOL || z > CPL_MATH_PI - SM_TOL) {
            dzp = 2 * dz / 3.;
            dzhp = dzp / 2.;
            if (z < SM_TOL) {
                zp = zmin + dzp;
            } else {
                zp = zmin + dzhp;
            }
            dap = da / 4.;
            ap = a;
            domegap = domega / 5.;
        } else {
            dzp = dz / 3.;
            dzhp = dzp / 2.;
            zp = zmin + dzhp;
            dap = da / 3.;
            ap0 = a - dap;
            ap = ap0;
            domegap = domega / 9.;
        }

        /* Insert new table rows */
        if (z < SM_TOL || z > CPL_MATH_PI - SM_TOL) {
            cpl_table_insert_window(grid, i, 4);
            nrow += 4;
        } else {
            cpl_table_insert_window(grid, i, 8);
            nrow += 8;
        }

        /* Write coordinates to CPL table */
        if (z < SM_TOL) {
            n = i;
            cpl_table_set(grid, "z", n, z);
            cpl_table_set(grid, "A", n, a);
            cpl_table_set(grid, "dzl", n, 0.);
            cpl_table_set(grid, "dzu", n, dzhp);
            cpl_table_set(grid, "dA", n, da);
            cpl_table_set(grid, "domega", n, domegap);
            for (k = 0; k < 4; k++) {
                n = i + k + 1;
                cpl_table_set(grid, "z", n, zp);
                cpl_table_set(grid, "A", n, ap);
                cpl_table_set(grid, "dzl", n, dzhp);
                cpl_table_set(grid, "dzu", n, dzhp);
                cpl_table_set(grid, "dA", n, dap);
                cpl_table_set(grid, "domega", n, domegap);
                ap += dap;
            }
        } else if (z > CPL_MATH_PI - SM_TOL) {
            for (k = 0; k < 4; k++) {
                n = i + k;
                cpl_table_set(grid, "z", n, zp);
                cpl_table_set(grid, "A", n, ap);
                cpl_table_set(grid, "dzl", n, dzhp);
                cpl_table_set(grid, "dzu", n, dzhp);
                cpl_table_set(grid, "dA", n, dap);
                cpl_table_set(grid, "domega", n, domegap);
                ap += dap;
            }
            n = i + 4;
            cpl_table_set(grid, "z", n, z);
            cpl_table_set(grid, "A", n, a);
            cpl_table_set(grid, "dzl", n, dzhp);
            cpl_table_set(grid, "dzu", n, 0.);
            cpl_table_set(grid, "dA", n, da);
            cpl_table_set(grid, "domega", n, domegap);
        } else {
            for (j = 0; j < 3; j++) {
                for (k = 0; k < 3; k++) {
                    n = i + 3 * j + k;
                    cpl_table_set(grid, "z", n, zp);
                    cpl_table_set(grid, "A", n, ap);
                    cpl_table_set(grid, "dzl", n, dzhp);
                    cpl_table_set(grid, "dzu", n, dzhp);
                    cpl_table_set(grid, "dA", n, dap);
                    cpl_table_set(grid, "domega", n, domegap);
                    ap += dap;
                }
                zp += dzp;
                ap = ap0;
            }
        }

        /* Next coordinates of initial grid */
        i = n + 1;

    }

    return CPL_ERROR_NONE;
}


double sm_scat_geteffcoldens(const double z, const double sig,
                             const int ns, const char scattype,
                             const double sm_hmin)
{
    /*!
     * \callgraph
     * Gets effective column density in \f${\rm cm}^{-2}\f$ for a path from
     * top of atmosphere to point of reference by considering the given
     * zenith distance. The output is either for molecules (scattype = 'r') or
     * aerosol particles (scattype = 'm'). There are two calculation modes
     * that can be selected: a numerical integration with ns + 1 path elements
     * in maximum or a faster approximation avoiding the integration (ns < 0).
     *
     * \b INPUT:
     * \param z         zenith distance in rad at point of reference
     * \param sig       height in km of point of reference relative to centre
     *                  of Earth
     * \param ns        number of path elements (-1 \f$\to\f$ approximation)
     * \param scattype  'r' = Rayleigh scattering,
     *                  'm' = Mie (aerosol) scattering
     *
     * \b RETURN:
     * - effective column density in \f${\rm cm}^{-2}\f$ for given path
     *
     * \b ERRORS:
     * - none
     */

    double b = 0.;

    /* Select calculation mode */
    if (ns > 0) {
        b = sm_scat_calceffcoldens(z, sig, ns, scattype);
    } else {
        b = sm_scat_esteffcoldens(z, sig - SM_R, SM_HMAX, scattype, sm_hmin);
    }

    /* Return effective column density */
    return b;
}


double sm_scat_calceffcoldens(const double z, const double sig, const int ns,
                              const char scattype)
{
    /*!
     * \callgraph
     * Calculates effective column density in \f${\rm cm}^{-2}\f$ for a path
     * from top of atmosphere to point of reference by considering the given
     * zenith distance. The size of the volume elements for the integration
     * increases with height in order to consider the decrease of density with
     * height.
     *
     * \b INPUT:
     * \param z         zenith distance in rad at point of reference
     * \param sig       height in km of point of reference relative to centre
     *                  of Earth
     * \param ns        maximum number of data points (- 1) along a path in
     *                  the atmosphere
     * \param scattype  'r' = Rayleigh scattering,
     *                  'm' = Mie (aerosol) scattering
     *
     * \b RETURN:
     * - effective column density in \f${\rm cm}^{-2}\f$ for given path
     *   (-1 in the case of errors)
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     */

    char errtxt[SM_MAXLEN+1];
    int j = 0, nsp = 0;
    double sigmax = 0., sinz = 0., cosz = 0., smax = 0., lnd = 0., dlns = 0.;
    double sp = 0., spmax = 0., spmin = 0., ds = 0., sigp = 0., nsig = 0.;
    double b = 0.;

    /* Top of atmosphere in km relative to centre of Earth */
    sigmax = SM_R + SM_HMAX;
    /* Sine and cosine of given zenith distance */
    sinz = sin(z);
    cosz = cos(z);
    /* Maximum distance of top of atmosphere from point of reference */
    smax = sqrt(sigmax * sigmax - sig * sig * sinz * sinz) - sig * cosz;
    /* Logarithmic, scaled path length to zenith */
    lnd = log((sigmax - sig) / SM_LNSCALE);
    /* Number of volume elements for integration (<= ns depending on height
       of point of reference */
    nsp = floor((double) ns * lnd / log(SM_HMAX / SM_LNSCALE) + 0.5);
    if (nsp == 0) nsp = 1;
    /* Length of logarithmic path element */
    dlns = log(smax / SM_LNSCALE) / (double) nsp;

    /* Compute effective column density */

    for (j = -1; j < nsp; j++) {
        /* Position along path */
        if (j == -1) {
            sp = SM_LNSCALE / 2.;
        } else {
            sp = SM_LNSCALE * exp(((double) j + 0.5) * dlns);
        }
        /* Upper limit of path element */
        spmax = SM_LNSCALE * exp((double) (j + 1) * dlns);
        /* Length of path element */
        ds = spmax - spmin;
        /* Upper limit becomes lower limit in next iteration */
        spmin = spmax;
        /* Height of path/volume element in km relative to centre of Earth */
        sigp = sqrt(sig * sig + sp * sp + 2 * sig * sp * cosz);
        /* Particle density for volume element in cm^-3 for Rayleigh or Mie
           (aerosol) scattering */
        if (scattype == 'r') {
            /* Rayleigh scattering (barometric formula) */
            nsig = sm_scat_getmolecdens(sigp - SM_R);
        } else if (scattype == 'm') {
            /* Mie (aerosol) scattering (rough exponential model) */
            nsig = sm_scat_getaerosoldens(sigp - SM_R);
        } else {
            sprintf(errtxt, "%s: scattype != 'r' and 'm'", SM_ERROR_IIP_TXT);
            cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
            return -1;
        }
        /* Add column density for given volume element */
        b += nsig * ds;
    }

    /* km in cm */
    b *= 1e5;

    /* Return effective column density */
    return b;
}


double sm_scat_esteffcoldens(const double z, const double hmin,
                             const double hmax, const char scattype,
                             const double sm_hmin)
{
    /*!
     * Estimates effective column density in \f${\rm cm}^{-2}\f$ for a path
     * between the two given heights. The input zenith distance is related to
     * the upper point as seen from the lower point. The approximation formula
     * estimates the effect of the curvature of the Earth and requires an
     * exponential decrease of the molecular or aerosol density with
     * increasing height. Up to the horizon (z = 90 deg), an accuracy better
     * than 12% is achieved. Reasonable results are also provided for zenith
     * distances greater than 90 deg.
     *
     * \b INPUT:
     * \param z         zenith distance in rad
     * \param hmin      lower height in km
     * \param hmax      upper height in km
     * \param scattype  'r' = Rayleigh scattering,
     *                  'm' = Mie (aerosol) scattering
     *
     * \b RETURN:
     * - effective column density in \f${\rm cm}^{-2}\f$ for given path
     *   (-1 in the case of errors)
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     */

    char errtxt[SM_MAXLEN+1];
    double sigmin = 0., n0 = 0., h0 = 0., hminc = 0., zcorr = 0., sigmax = 0.;
    double smax = 0., s = 0., sigminc = 0., sig = 0., xeff = 0., b = 0.;

    /* Lower height relative to centre of Earth */
    sigmin = hmin + SM_R;

    /* Check input parameters */
    if (hmin >= hmax * (1 - SM_TOL)) {
        sprintf(errtxt, "%s: hmin >= hmax", SM_ERROR_IIP_TXT);
        cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
        return -1;
    } else if (z > CPL_MATH_PI_2) {
        if (z > CPL_MATH_PI_2 + acos((sm_hmin + SM_R) / sigmin)) {
            sprintf(errtxt, "%s: z > z_max (ground)", SM_ERROR_IIP_TXT);
            cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
            return -1;
        }
    }

    /* Rayleigh or Mie (aerosol) scattering? */
    if (scattype == 'r') {
        /* Rayleigh scattering */
        n0 = SM_N0_R;
        h0 = SM_H0_R;
    } else if (scattype == 'm') {
        /* Mie (aerosol) scattering */
        n0 = SM_N0_M;
        h0 = SM_H0_M;
    } else {
        sprintf(errtxt, "%s: scattype != 'r' and 'm'", SM_ERROR_IIP_TXT);
        cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
        return -1;
    }

    /* Correct zenith distance for curvature of Earth */
    if (z <= SM_ZMAX * CPL_MATH_RAD_DEG) {
        hminc = hmin;
        zcorr = 0.;
    } else {
        sigmax = sigmin + h0;
        smax = sqrt(sigmax * sigmax - sigmin * sigmin);
        if (z == CPL_MATH_PI_2) {
            hminc = hmin;
            s = smax;
        } else if (z > CPL_MATH_PI_2) {
            sigminc = sigmin * sin(CPL_MATH_PI - z);
            hminc = sigminc - SM_R;
            s = smax;
        } else {
            hminc = hmin;
            s = h0 * tan(z);
            if (s > smax) {
                s = smax;
            }
        }
        sig = sqrt(sigmin * sigmin + s * s);
        //zcorr = atan((sig - sigmin) / s);
        zcorr = (sig - sigmin) / s;

    }

    /* Estimate effective airmass */
    if (z > CPL_MATH_PI_2) {
        xeff = 1 / cos(CPL_MATH_PI_2 - zcorr);
    } else {
        xeff = 1 / cos(z - zcorr);
    }

    /* Calculate effective column density */
    if (hmax >= SM_HMAX) {
        b = 1e5 * n0 * h0 * exp(-hminc / h0) * xeff;
    } else {
        b = 1e5 * n0 * h0 * (exp(-hminc / h0) - exp(-hmax / h0)) * xeff;
    }

    /* Return effective column density */
    return b;
}


double sm_scat_getmolecdens(const double h)
{
    /*!
     * Provides molecular density in \f${\rm cm}^{-3}\f$ for given height in
     * km. Uses barometric formula. The constants are from Staude (1975).
     *
     * \b INPUT:
     * \param h  height above ground in km
     *
     * \b RETURN:
     * - molecular density in \f${\rm cm}^{-3}\f$
     *
     * \b ERRORS:
     * - none
     */

    double n = 0.;

    n = SM_N0_R * exp(-h / SM_H0_R);

    return n;
}


double sm_scat_getaerosoldens(const double h)
{
    /*!
     * Provides aerosol density in \f${\rm cm}^{-3}\f$ for given height in km.
     * Uses exponential formula. The constants are from Staude (1975).
     *
     * \b INPUT:
     * \param h  height above ground in km
     *
     * \b RETURN:
     * - aerosol density in \f${\rm cm}^{-3}\f$
     *
     * \b ERRORS:
     * - none
     */

    double n = 0.;

    n = SM_N0_M * exp(-h / SM_H0_M);

    return n;
}


double sm_scat_getgrefl(const double lam)
{
    /*!
     * Returns ground reflectance for given wavelength. The reflectance data
     * are from Sutter et al. (2007), who studied Atacama desert soils. The
     * spectrum for sample A1 is scaled to match reflectances of 0.09, 0.10,
     * and 0.12 at 0.354, 0.388, and 0.5 \f$\mu{\rm m}\f$, respectively, as
     * estimated from OMI data.
     *
     * \b INPUT:
     * \param lam0  wavelength in \f$\mu{\rm m}\f$
     *
     * \b RETURN:
     * - ground reflectance
     *
     * \b ERRORS:
     * - none
     */

    int i = 0;
    double grefl = 0.;

    /* Number of data points */
    int n = 26;
    /* Wavelengths */
    double l[] =  {0.350,0.354,0.381,0.388,0.412,0.442,0.473,0.500,0.504,
                   0.535,0.565,0.596,0.658,0.719,0.781,0.842,0.904,0.965,
                   1.027,1.088,1.150,1.212,1.273,1.335,1.581,1.735};
    /* Ground reflectance from Sutter et al. */
    double gr[] = {0.066,0.066,0.066,0.068,0.075,0.086,0.096,0.107,0.109,
                   0.124,0.154,0.184,0.208,0.226,0.233,0.235,0.238,0.240,
                   0.244,0.246,0.248,0.250,0.253,0.255,0.261,0.265};
    /* Correction factors to match OMI data at 354, 388, and 500 nm */
    double c[] =  {1.364,1.364,1.449,1.471,1.396,1.302,1.205,1.121,1.121,
                   1.121,1.121,1.121,1.121,1.121,1.121,1.121,1.121,1.121,
                   1.121,1.121,1.121,1.121,1.121,1.121,1.121,1.121};

    /* Get interpolated ground reflectance from list (constant extrapolation
       if input wavelength is outside covered range) */
    if (lam <= l[0]) {
        grefl = gr[0] * c[0];
    } else if (lam > l[n-1]) {
        grefl = gr[n-1] * c[n-1];
    } else {
        while (lam > l[i]) {
            i++;
        }
        grefl = (gr[i] * c[i] - gr[i-1] * c[i-1]) * (lam - l[i-1]) /
                (l[i] - l[i-1]) + gr[i-1] * c[i-1];
    }

    /* Return ground reflectance */
    return grefl;
}


cpl_error_code sm_scat_corrabs(smspec *spec, const cpl_table *xefftab,
                               const smspec *molabstrans)
{
    /*!
     * Corrects spectrum of scattering intensities for molecular absorption.
     * The scattering intensities are roughly corrected by means of effective
     * airmasses corresponding to the extinction optical depths based on
     * scattering and absorption.
     *
     * \b INPUT:
     * \param spec         spectrum of scattering intensities
     * \param xefftab      effective airmasses for different optical depths
     * \param molabstrans  transmission curve from radiative transfer code
     *
     * \b OUTPUT:
     * \param spec         corrected spectrum of scattering intensities
     *
     * \b ERRORS:
     * - IDG: Inconsistent data grids
     */

    char errtxt[SM_MAXLEN+1];
    int n = 0, i = 0, k0 = 0, k = 0;
    const double *tau, *xeff;
    double trans = 0., tau0 = 0., frac = 0., xeff0 = 0., corr = 0.;

    /* Test agreement of wavelength grids in input spectra */
    if ((int) sm_spec_compgrids(molabstrans, spec) == (int) SM_ERROR_IDG) {
        sprintf(errtxt, "%s: smspec *molabstrans != smspec *spec",
                SM_ERROR_IDG_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IDG, "%s", errtxt);
    }

    /* Number of rows in CPL table */
    n = cpl_table_get_nrow(xefftab);

    /* Get pointers to table columns */
    tau = cpl_table_get_data_double_const(xefftab, "tau0");
    xeff = cpl_table_get_data_double_const(xefftab, "X_eff");

    /* Calculate scattering intensity for each wavelength of the input grid */

    for (i = 0; i < spec->n; i++) {

        /* Get transmission for both scattering processes and molecular
           absorption */
        trans = molabstrans->dat[i].flux;

        /* Handle zero transmission */
        if (trans <= 0.) {
            spec->dat[i].flux = 0.;
            continue;
        }

        /* Convert transmission into optical depth */
        tau0 = -log(trans);

        /* Avoid negative tau values */
        if (tau0 < 0) {
            tau0 = 0;
        }

        /* Find index for first higher optical depth in table */

        for (k0 = -1, k = 0; k < n; k++) {
            if (tau[k] > tau0) {
                k0 = k;
                break;
            }
        }

        if (k0 == -1) {
            k0 = n;
        }

        /* Get pixel fraction */
        if (k0 != 0 && k0 != n) {
            frac = (tau0 - tau[k0]) / (tau[k0-1] - tau[k0]);
        }

        /* Interpolate effective airmass for given optical depth */
        if (k0 == 0) {
            xeff0 = xeff[k0];
        } else if (k0 == n) {
            xeff0 = xeff[k0-1];
        } else {
            xeff0 = frac * (xeff[k0-1] - xeff[k0]) + xeff[k0];
        }

        /* Avoid effective airmasses below 1 */
        if (xeff0 < 1.) {
            xeff0 = 1.;
        }

        /* Correction for full extinction along light path (rough
           approximation) */
        corr = pow(trans, xeff0);

        /* Write corrected scattering intensity into output spectrum */
        spec->dat[i].flux = spec->dat[i].flux * corr;

    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_scat_estmultiscat(const char *sscatcorfile,
                                    const char *dscatcorfile,
                                    const smparmodel modelpar)
{
    /*!
     * \callgraph
     * Performs single and double scattering calculations of moonlight for a
     * list of parameter sets consisting of realistic combinations of
     * line-of-sight and Moon zenith distances, angles between line of sight
     * and Moon \f$\rho\f$, and wavelength. The results are used to estimate
     * multiple scattering corrections for single and double scattering
     * corrections. For single scattering corrections, a
     * wavelength-\f$\rho\f$ grid with correction factors is written to an
     * ASCII file. Another file is written for triple and higher order
     * scattering. In the latter case, weighted multiple scattering
     * intensities depending on wavelength are provided.
     *
     * \b INPUT:
     * \param sscatcorfile  path and name of file with multiple scattering
     *                      correction factors for single scattering
     *                      calculations
     * \param dscatcorfile  path and name of file with multiple scattering
     *                      intensities for correction of double scattering
     *                      calculations
     * \param modelpar      sky emission parameters (see typedef of
     *                      ::smparmodel)
     *
     * \b OUTPUT:
     * - none
     *
     * \b ERRORS:
     * - none
     */

    cpl_table *dsgrid, *ztab, *zmtab, *rhotab, *datatab, *zavtab, *rhoavtab;
    smgrid miephase;
    smspec lamgrid, rscattrans, mscattrans;
    smscatres scatres;
    char miephasefile[FILENAME_MAX];
    int mode = SM_MSCATMODE; // calculation mode (see header)
    int nz = 0, n = 0, i = 0;
    double sig0 = 0., b0r = 0., b0m = 0., lam0 = 0., rtrans = 0., mtrans = 0.;
    double tau0r = 0., tau0m = 0., tau0 = 0., cextr = 0., cextm = 0.;
    double cexta = 0., z = 0., zmoon = 0., rho = 0., zrad = 0., zmrad = 0.;
    double rhorad = 0., flux0 = 1.;

    /* Height of observer in km relative to centre of Earth */
    sig0 = SM_R + modelpar.sm_h;

    /* Get effective column density in cm^-2 for zenith */
    b0r = sm_scat_geteffcoldens(0., sig0, SM_NS, 'r', modelpar.sm_hmin);
    b0m = sm_scat_geteffcoldens(0., sig0, SM_NS, 'm', modelpar.sm_hmin);

    /* Create coordinate grid for double scattering calculations if
       requested */
    dsgrid = cpl_table_new(0);
    if (mode == 0) {
        nz = SM_NZ;
    } else {
        nz = SM_NZD;
    }
    sm_scat_initscatgrid(dsgrid, nz);

    /* Create smspec structure for wavelength grid */
    sm_scat_createlambdagrid(&lamgrid, NULL);

    /* Calculate Rayleigh and Mie (aerosol) extinction curves */
    sm_spec_copy(&rscattrans, &lamgrid);
    sm_comp_calcrayleighscat(&rscattrans, modelpar);
    sm_spec_copy(&mscattrans, &lamgrid);
    sm_comp_calcmiescat(&mscattrans, modelpar);

    /* Read table for Mie scattering phase functions */
    sprintf(miephasefile, "%s/%s", modelpar.datapath, modelpar.miephasename);
    sm_grid_read(&miephase, miephasefile);

    /* Create table for line-of-sight zenith distances and their weights */
    ztab = cpl_table_new(0);
    sm_scat_createztab(ztab);

    /* Create table for Moon zenith distances and their weights */
    zmtab = cpl_table_new(0);
    sm_scat_createzmtab(zmtab);

    /* Create table for angles between line of sight and Moon and their
       weights */
    rhotab = cpl_table_new(0);
    sm_scat_createrhotab(rhotab);

    /* Create columns and rows for input parameters and results table */
    datatab = cpl_table_new(0);
    sm_scat_createdatatab(datatab, ztab, zmtab, rhotab, &lamgrid);

    /* Get number of parameter combinations */
    n = cpl_table_get_nrow(datatab);

    /* Calculations for scattering */

    for (i = 0; i < n; i++) {

        /* Counter for stdout */
        printf("%d/%d\r", i, n);

        /* Get transmission values by Rayleigh and Mie scattering for listed
           wavelength */
        lam0 = cpl_table_get(datatab, "lam", i, NULL);
        sm_spec_getval(&rtrans, &rscattrans, lam0);
        sm_spec_getval(&mtrans, &mscattrans, lam0);

        /* Calculate optical depths for listed wavelength */
        tau0r = -log(rtrans);
        tau0m = -log(mtrans);
        tau0 = tau0r + tau0m;

        /* Correct for numerical uncertainties */
        if (lam0 < miephase.xpos[0]) {
            lam0 = miephase.xpos[0];
        }
        if (lam0 > miephase.xpos[(miephase.nx)-1]) {
            lam0 = miephase.xpos[(miephase.nx)-1];
        }

        /* Calculate Rayleigh and aerosol extinction cross section in cm^2 */
        cextr = tau0r / b0r;
        cextm = tau0m / b0m;

        /* No molecular absorption */
        cexta = 0;

        /* Get line-of-sight and Moon zenith distance and angle between line
           of sight and Moon */
        z = cpl_table_get(datatab, "z", i, NULL);
        zmoon = cpl_table_get(datatab, "zmoon", i, NULL);
        rho = cpl_table_get(datatab, "rho", i, NULL);

        /* Conversion deg to rad */
        zrad = CPL_MATH_RAD_DEG * z;
        zmrad = CPL_MATH_RAD_DEG * zmoon;
        rhorad = CPL_MATH_RAD_DEG * rho;

        /* Default values of output parameters */
        sm_scat_initscatres(&scatres);

        /* Calculate intensity of radiation scattered into line of sight */
        sm_scat_calcscat(&scatres, modelpar.sm_h, modelpar.sm_hmin, zrad, zmrad, rhorad, cextr, cextm, cexta,
                         dsgrid, &miephase, lam0, mode);

        /* Process results of scattering calculations (retrieve final
           intensities and the effective airmass) */
        sm_scat_procscatres(&scatres, flux0, lam0, tau0, cextr, cextm,
                            modelpar.ssa);

        /* Copy scattering results in smscatres structure into data table */
        sm_scat_copyscatres(datatab, i, &scatres);

    }

    /* Average scattering intensities for different zentith distances of line
       of sight and Moon */
    zavtab = cpl_table_new(0);
    sm_scat_createzavtab(zavtab, rhotab, &lamgrid);
    sm_scat_calczavdoublescat(zavtab, datatab, ztab, zmtab);

    /* Average scattering intensities for different zentith distances of line
       of sight and Moon */
    rhoavtab = cpl_table_new(0);
    sm_scat_createrhoavtab(rhoavtab, &lamgrid);
    sm_scat_calcrhoavdoublescat(rhoavtab, zavtab, rhotab);

    /* Calculate multiple scattering corrections */
    sm_scat_calcscatcorr(zavtab, rhoavtab);

    /* Write ASCII files for multiple scattering corrections of single and
       double scattering calculations */
    sm_scat_writesinglescatcorr(sscatcorfile, zavtab, rhotab, &lamgrid);
    sm_scat_writedoublescatcorr(dscatcorfile, rhoavtab);

    /* Write CPL tables into FITS files for analysis */
    cpl_table_save(datatab, NULL, NULL, "output/scattab_all.fits",
                   CPL_IO_CREATE);
    cpl_table_save(zavtab, NULL, NULL, "output/scattab_zav.fits",
                   CPL_IO_CREATE);
    cpl_table_save(rhoavtab, NULL, NULL, "output/scattab_rhoav.fits",
                   CPL_IO_CREATE);

    /* Free memory */
    sm_spec_free(&lamgrid);
    sm_spec_free(&rscattrans);
    sm_spec_free(&mscattrans);
    cpl_table_delete(dsgrid);
    cpl_table_delete(ztab);
    cpl_table_delete(zmtab);
    cpl_table_delete(rhotab);
    cpl_table_delete(datatab);
    cpl_table_delete(zavtab);
    cpl_table_delete(rhoavtab);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_scat_createztab(cpl_table *ztab)
 {
    /*!
     * Writes CPL table with list of line-of-sight zenith distances and their
     * weights
     *
     * \b INPUT:
     * \param ztab  empty CPL table
     *
     * \b OUTPUT:
     * \param ztab  CPL table with line-of-sight zenith distances and their
     *              weights
     *
     * \b ERRORS:
     * - none
     */

    int i = 0;

    /* Values to be set */
    int nz = 8;
    double z[]  = {   0.,  10.,  20.,  30.,  40.,  50.,  60.,  70.};
    double pz[] = {0.006,0.097,0.255,0.272,0.222,0.106,0.038,0.004};

    /* Create table columns */
    cpl_table_new_column(ztab, "z", CPL_TYPE_DOUBLE);
    cpl_table_new_column(ztab, "pz", CPL_TYPE_DOUBLE);

    /* Set table size */
    cpl_table_set_size(ztab, nz);

    /* Write values into CPL table */
    for (i = 0; i < nz; i++) {
        cpl_table_set(ztab, "z", i, z[i]);
        cpl_table_set(ztab, "pz", i, pz[i]);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_scat_createzmtab(cpl_table *zmtab)
 {
    /*!
     * Writes CPL table with list of Moon zenith distances and their weights
     *
     * \b INPUT:
     * \param zmtab  empty CPL table
     *
     * \b OUTPUT:
     * \param zmtab  CPL table with Moon zenith distances and their weights
     *
     * \b ERRORS:
     * - none
     */

    int i = 0;

    /* Values to be set */
    int nzm = 9;
    double zm[]  = {   0.,  10.,  20.,  30.,  40.,  50.,  60.,  70.,  80.};
    double pzm[] = {0.009,0.046,0.072,0.097,0.131,0.176,0.142,0.133,0.194};

    /* Create table columns */
    cpl_table_new_column(zmtab, "zm", CPL_TYPE_DOUBLE);
    cpl_table_new_column(zmtab, "pzm", CPL_TYPE_DOUBLE);

    /* Set table size */
    cpl_table_set_size(zmtab, nzm);

    /* Write values into CPL table */
    for (i = 0; i < nzm; i++) {
        cpl_table_set(zmtab, "zm", i, zm[i]);
        cpl_table_set(zmtab, "pzm", i, pzm[i]);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_scat_createrhotab(cpl_table *rhotab)
 {
    /*!
     * Writes CPL table with list of angles between line of sight and Moon
     * and their weights
     *
     * \b INPUT:
     * \param rhotab  empty CPL table
     *
     * \b OUTPUT:
     * \param rhotab  CPL table with angles between line of sight and Moon and
     *                their weights
     *
     * \b ERRORS:
     * - none
     */

    int i = 0;

    /* Values to be set */
    int nrho = 14;
    double rho[]  = {  10.,  20.,  30.,  40.,  50.,  60.,  70.,  80.,  90.,
                      100., 110., 120., 130., 140.};
    //double prho[] = {0.042,0.075,0.102,0.121,0.135,0.131,0.121,0.099,0.073,
    //                 0.053,0.032,0.013,0.004,0.001};
    double prho[] = {0.000,0.000,0.074,0.143,0.160,0.156,0.143,0.117,0.087,
                     0.063,0.037,0.015,0.005,0.001};

    /* Create table columns */
    cpl_table_new_column(rhotab, "rho", CPL_TYPE_DOUBLE);
    cpl_table_new_column(rhotab, "prho", CPL_TYPE_DOUBLE);

    /* Set table size */
    cpl_table_set_size(rhotab, nrho);

    /* Write values into CPL table */
    for (i = 0; i < nrho; i++) {
        cpl_table_set(rhotab, "rho", i, rho[i]);
        cpl_table_set(rhotab, "prho", i, prho[i]);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_scat_createdatatab(cpl_table *datatab,
                                     const cpl_table *ztab,
                                     const cpl_table *zmtab,
                                     const cpl_table *rhotab,
                                     const smspec *lamgrid)
 {
    /*!
     * Writes CPL table with list of line-of-sight and Moon zenith distances,
     * angles between line of sight and Moon, and wavelengths in
     * \f$\mu{\rm m}\f$ for atmospheric scattering calculations for moonlight.
     * The first three parameters are read from CPL tables, whereas the set
     * of wavelengths is taken from an ::smspec structure.
     *
     * \b INPUT:
     * \param datatab  empty CPL table
     * \param ztab     grid of line-of-sight zenith distances
     * \param zmtab    grid of Moon zenith distances
     * \param rhotab   grid of angles between line of sight and Moon
     * \param lamgrid  grid of wavelengths
     *
     * \b OUTPUT:
     * \param datatab  list of parameter sets characterising the intensity
     *                 distribution of scattered moonlight on the sky
     *
     * \b ERRORS:
     * - none
     */

    int nz = 0, nzm = 0, nr = 0;
    int h = 0, i = 0, j = 0, k = 0, n = 0;
    const double *z = NULL, *zmoon = NULL, *rho = NULL;

    /* Create table columns */
    cpl_table_new_column(datatab, "z", CPL_TYPE_DOUBLE);
    cpl_table_new_column(datatab, "zmoon", CPL_TYPE_DOUBLE);
    cpl_table_new_column(datatab, "rho", CPL_TYPE_DOUBLE);
    cpl_table_new_column(datatab, "lam", CPL_TYPE_DOUBLE);
    cpl_table_new_column(datatab, "I_1s", CPL_TYPE_DOUBLE);
    cpl_table_new_column(datatab, "I_2s", CPL_TYPE_DOUBLE);
    cpl_table_new_column(datatab, "I_u2s", CPL_TYPE_DOUBLE);

    /* Get size of parameter tables */
    nz = cpl_table_get_nrow(ztab);
    nzm = cpl_table_get_nrow(zmtab);
    nr = cpl_table_get_nrow(rhotab);

    /* Set maximum table size */
    cpl_table_set_size(datatab, nz * nzm * nr * lamgrid->n);

    /* Get pointers to parameter values in the different input tables */
    z = cpl_table_get_data_double_const(ztab, "z");
    zmoon = cpl_table_get_data_double_const(zmtab, "zm");
    rho = cpl_table_get_data_double_const(rhotab, "rho");

    /* Write parameters to CPL table */

    for (h = 0; h < nz; h++) {

        for (i = 0; i < nzm; i++) {

            for (j = 0; j < nr; j++) {

                /* Skip impossible parameter combinations */
                if (fabs(z[h] - zmoon[i]) > rho[j] ||
                    z[h] + zmoon[i] < rho[j]) {
                    continue;
                }

                for (k = 0; k < lamgrid->n; k++) {

                    /* Write parameter values into CPL table */
                    cpl_table_set(datatab, "z", n, z[h]);
                    cpl_table_set(datatab, "zmoon", n, zmoon[i]);
                    cpl_table_set(datatab, "rho", n, rho[j]);
                    cpl_table_set(datatab, "lam", n, lamgrid->dat[k].lam);

                    /* Count number of valid parameter sets */
                    n++;

                }

            }

        }

    }

    /* Resize table to resulting number of parameter sets */
    cpl_table_set_size(datatab, n);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_scat_copyscatres(cpl_table *datatab, const int row,
                                   const smscatres *scatres)
{
    /*!
     * Copies scattering results from ::smscatres structure to given row of
     * CPL table
     *
     * \b INPUT:
     * \param scatres   empty ::smscatres structure
     *
     * \b OUPUT:
     * \param scatres   initialised ::smscatres structure
     *
     * \b ERRORS:
     * - none
     */

    char errtxt[SM_MAXLEN+1];

    /* Check existence of row */
    if (row < 0 || row >= cpl_table_get_nrow(datatab)) {
        sprintf(errtxt, "%s: invalid row %d for cpl_table *datatab",
                SM_ERROR_IIP_TXT, row);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Copy data */
    cpl_table_set(datatab, "I_1s", row, scatres->iscat1s);
    cpl_table_set(datatab, "I_2s", row, scatres->iscat2s);
    cpl_table_set(datatab, "I_u2s", row, scatres->iscat);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_scat_createzavtab(cpl_table *zavtab,
                                    const cpl_table *rhotab,
                                    const smspec *lamgrid)
{
    /*!
     * Writes CPL table with list of angles between line of sight and Moon
     * \f$\rho\f$ and wavelengths in \f$\mu{\rm m}\f$ for atmospheric
     * scattering calculations for moonlight. The former parameter is read
     * from a CPL table, whereas the set of wavelengths is taken from an
     * ::smspec structure.
     *
     * \b INPUT:
     * \param zavtab   empty CPL table
     * \param rhotab   grid of angles between line of sight and Moon
     * \param lamgrid  grid of wavelengths
     *
     * \b OUTPUT:
     * \param zavtab   list of parameter sets depending on angle between line
     *                 of sight and Moon and wavelength
     *
     * \b ERRORS:
     * - none
     */

    int nr = 0, j = 0, k = 0, n = 0;
    const double *rho = NULL;

    /* Create table columns */
    cpl_table_new_column(zavtab, "rho", CPL_TYPE_DOUBLE);
    cpl_table_new_column(zavtab, "lam", CPL_TYPE_DOUBLE);
    cpl_table_new_column(zavtab, "I_1s", CPL_TYPE_DOUBLE);
    cpl_table_new_column(zavtab, "I_2s", CPL_TYPE_DOUBLE);
    cpl_table_new_column(zavtab, "I_u2s", CPL_TYPE_DOUBLE);

    /* Get size of rho table */
    nr = cpl_table_get_nrow(rhotab);

    /* Set table size */
    cpl_table_set_size(zavtab, nr * lamgrid->n);

    /* Get pointer to rho values in the rho table */
    rho = cpl_table_get_data_double_const(rhotab, "rho");

    /* Write parameters to CPL table */

    for (j = 0; j < nr; j++) {

        for (k = 0; k < lamgrid->n; k++) {

            /* Write parameter values into CPL table */
            cpl_table_set(zavtab, "rho", n, rho[j]);
            cpl_table_set(zavtab, "lam", n, lamgrid->dat[k].lam);

            /* Count parameter sets */
            n++;

        }

    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_scat_calczavdoublescat(cpl_table *zavtab,
                                         cpl_table *datatab,
                                         const cpl_table *ztab,
                                         const cpl_table *zmtab)
{
    /*!
     * Calculates scattering intensities weighted for line-of-sight and Moon
     * zenith distances. The data to be averaged has to be provided by a
     * table which contains scattering intensities for the two zenith
     * distances, the angle between line of sight and Moon \f$\rho\f$, and the
     * wavelength. The weights for the different zenith distances are taken
     * from two dedicated tables. The ouput table provides averaged scattering
     * intensities for different \f$\rho\f$ and wavelengths.
     *
     * \b INPUT:
     * \param zavtab   CPL table with output parameter set
     * \param datatab  scattering intensities from single and double
     *                 scattering calculations
     * \param ztab     grid and weights of line-of-sight zenith distances
     * \param zmtab    grid and weights of Moon zenith distances
     *
     * \b OUTPUT:
     * \param zavtab   CPL table with scattering intensities weighted for
     *                 line-of-sight and Moon zenith distances
     *
     * \b ERRORS:
     * - none
     */

    cpl_table *seltab;
    int n = 0, nd = 0, nz = 0, nzm = 0, i = 0, ns = 0, h = 0, j = 0, k = 0;
    int l = 0;
    double *rho, *lam, *zs, *zms, *ps;
    const double *z, *pz, *zm, *pzm;
    double psum = 0., av = 0.;

    /* Columns in data table to be averaged */
    int ncol = 3;
    char col[][SM_LENLINE+1] = {"I_1s", "I_2s", "I_u2s"};

    /* Get number of rows in output table */
    n = cpl_table_get_nrow(zavtab);

    /* Get pointers to angles between line of sight and Moon rho and
       wavelengths in output table */
    rho = cpl_table_get_data_double(zavtab, "rho");
    lam = cpl_table_get_data_double(zavtab, "lam");

    /* Get number of rows in data table */
    nd = cpl_table_get_nrow(datatab);

    /* Create temorary columns for weights and weighted averages in data
       table */
    cpl_table_new_column(datatab, "p", CPL_TYPE_DOUBLE);
    cpl_table_new_column(datatab, "av", CPL_TYPE_DOUBLE);

    /* Initialise weight column in data table */
    cpl_table_fill_column_window(datatab, "p", 0, nd, 0.);

    /* Get number of line-of-sight and Moon zenith distances */
    nz = cpl_table_get_nrow(ztab);
    nzm = cpl_table_get_nrow(zmtab);

    /* Get pointers to parameter values in the tables listing the different
       zenith distances */
    z = cpl_table_get_data_double_const(ztab, "z");
    pz = cpl_table_get_data_double_const(ztab, "pz");
    zm = cpl_table_get_data_double_const(zmtab, "zm");
    pzm = cpl_table_get_data_double_const(zmtab, "pzm");

    /* Average scattering intensities for each combination of angle rho and
       wavelength */

    for (i = 0; i < n; i++) {

        /* Select rows in data table depending on angle rho and wavelength */
        cpl_table_unselect_all(datatab);
        cpl_table_or_selected_double(datatab, "rho", CPL_EQUAL_TO, rho[i]);
        cpl_table_and_selected_double(datatab, "lam", CPL_EQUAL_TO, lam[i]);
        seltab = cpl_table_extract_selected(datatab);

        /* Check existence of parameter sets for the given selection
           criteria */
        ns = cpl_table_get_nrow(seltab);
        if (ns == 0) {
            cpl_table_delete(seltab);
            continue;
        }

        /* Get pointers to columns in selected table */
        zs = cpl_table_get_data_double(seltab, "z");
        zms = cpl_table_get_data_double(seltab, "zmoon");
        ps = cpl_table_get_data_double(seltab, "p");

        /* Find weight for each entry in selected table */
        for (h = 0, j = 0; j < nz; j++) {
            for (k = 0; k < nzm; k++) {
                if (z[j] != zs[h] || zm[k] != zms[h]) {
                    /* Skip non-existing parameter sets */
                    continue;
                }
                ps[h] = pz[j] * pzm[k];
                h++;
                if (h == ns) break;
            }
        }

        /* Make sure that summed weight for selected table is 1 */
        psum = ns * cpl_table_get_column_mean(seltab, "p");
        cpl_table_divide_scalar(seltab, "p", psum);

        /* Calculate weighted average for different scattering parameters */
        for (l = 0; l < ncol; l++) {
            cpl_table_fill_column_window(seltab, "av", 0, ns, 0.);
            cpl_table_add_columns(seltab, "av", col[l]);
            cpl_table_multiply_columns(seltab, "av", "p");
            av = ns * cpl_table_get_column_mean(seltab, "av");
            cpl_table_set(zavtab, col[l], i, av);
        }

        /* Delete temporary table */
        cpl_table_delete(seltab);

    }

    /* Select all entries in data table */
    cpl_table_select_all(datatab);

    /* Delete temporary columns in data table */
    cpl_table_erase_column(datatab, "p");
    cpl_table_erase_column(datatab, "av");

    return CPL_ERROR_NONE;
}


cpl_error_code sm_scat_createrhoavtab(cpl_table *rhoavtab,
                                      const smspec *lamgrid)
{
    /*!
     * Writes CPL table with list of wavelengths in \f$\mu{\rm m}\f$ for
     * atmospheric scattering calculations for moonlight. The set of
     * wavelengths is taken from an ::smspec structure.
     *
     * \b INPUT:
     * \param rhoavtab  empty CPL table
     * \param lamgrid   grid of wavelengths
     *
     * \b OUTPUT:
     * \param rhoavtab  table of wavelengths plus empty columns for scattering
     *                  intensities
     *
     * \b ERRORS:
     * - none
     */

    int k = 0;

    /* Create table columns */
    cpl_table_new_column(rhoavtab, "lam", CPL_TYPE_DOUBLE);
    cpl_table_new_column(rhoavtab, "I_1s", CPL_TYPE_DOUBLE);
    cpl_table_new_column(rhoavtab, "I_2s", CPL_TYPE_DOUBLE);

    /* Set table size */
    cpl_table_set_size(rhoavtab, lamgrid->n);

    /* Write wavelengths to CPL table */
    for (k = 0; k < lamgrid->n; k++) {
        cpl_table_set(rhoavtab, "lam", k, lamgrid->dat[k].lam);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_scat_calcrhoavdoublescat(cpl_table *rhoavtab,
                                           cpl_table *zavtab,
                                           const cpl_table *rhotab)
{
    /*!
     * Calculates weighted scattering intensities depending on wavelength and
     * provides multiple scattering intensities and correction factors. The
     * data to be averaged has to be provided by a table which contains
     * scattering intensities for the angle between line of sight and Moon
     * \f$\rho\f$ and wavelength. The weights for the different \f$\rho\f$
     * are taken from a dedicated table. The \f$\rho\f$-averaged ouput table
     * provides averaged single and double scattering intensities for the
     * given wavelength grid.
     *
     * \b INPUT:
     * \param rhoavtab  CPL table with wavelength grid and required columns
     *                  for calculations
     * \param zavtab    CPL table with scattering intensities weighted for
     *                  line-of-sight and Moon zenith distances
     * \param rhotab    grid and weights of angles between line of sight and
     *                  Moon
     *
     * \b OUTPUT:
     * \param rhoavtab  CPL table with weighted scattering intensities
     *                  depending on wavelength
     *
     * \b ERRORS:
     * - none
     */

    cpl_table *seltab;
    int n = 0, nd = 0, nrho = 0, i = 0, ns = 0, h = 0, j = 0, l = 0;
    double *lam, *rhos, *ps;
    const double *rho, *prho;
    double psum = 0., av = 0.;

    /* Columns in data table to be averaged */
    int ncol = 2;
    char col[][SM_LENLINE+1] = {"I_1s", "I_2s"};

    /* Get number of rows in output table */
    n = cpl_table_get_nrow(rhoavtab);

    /* Get pointer to wavelengths in output table */
    lam = cpl_table_get_data_double(rhoavtab, "lam");

    /* Get number of rows in data table */
    nd = cpl_table_get_nrow(zavtab);

    /* Create temorary columns for weights and weighted averages in data
       table */
    cpl_table_new_column(zavtab, "p", CPL_TYPE_DOUBLE);
    cpl_table_new_column(zavtab, "av", CPL_TYPE_DOUBLE);

    /* Initialise weight column in data table */
    cpl_table_fill_column_window(zavtab, "p", 0, nd, 0.);

    /* Get number of angles between line of sight and Moon rho */
    nrho = cpl_table_get_nrow(rhotab);

    /* Get pointers to parameter values in rho table */
    rho = cpl_table_get_data_double_const(rhotab, "rho");
    prho = cpl_table_get_data_double_const(rhotab, "prho");

    /* Average scattering intensities for each combination of angle rho and
       wavelength */

    for (i = 0; i < n; i++) {

        /* Select rows in data table depending on wavelength */
        cpl_table_unselect_all(zavtab);
        cpl_table_or_selected_double(zavtab, "lam", CPL_EQUAL_TO, lam[i]);
        seltab = cpl_table_extract_selected(zavtab);

        /* Check existence of parameter sets for the given selection
           criteria */
        ns = cpl_table_get_nrow(seltab);
        if (ns == 0) {
            cpl_table_delete(seltab);
            continue;
        }

        /* Get pointers to rho and weight column in selected table */
        rhos = cpl_table_get_data_double(seltab, "rho");
        ps = cpl_table_get_data_double(seltab, "p");

        /* Find weight for each entry in selected table */
        for (h = 0, j = 0; j < nrho; j++) {
            if (rho[j] != rhos[h]) {
                /* Skip non-existing parameter sets */
                continue;
            }
            ps[h] = prho[j];
            h++;
            if (h == ns) break;
        }

        /* Make sure that summed weight for selected table is 1 */
        psum = ns * cpl_table_get_column_mean(seltab, "p");
        cpl_table_divide_scalar(seltab, "p", psum);

        /* Calculate weighted average for different scattering parameters */
        for (l = 0; l < ncol; l++) {
            cpl_table_fill_column_window(seltab, "av", 0, ns, 0.);
            cpl_table_add_columns(seltab, "av", col[l]);
            cpl_table_multiply_columns(seltab, "av", "p");
            av = ns * cpl_table_get_column_mean(seltab, "av");
            cpl_table_set(rhoavtab, col[l], i, av);
        }

        /* Delete temporary table */
        cpl_table_delete(seltab);

    }

    /* Select all entries in data table */
    cpl_table_select_all(zavtab);

    /* Delete temporary columns in data table */
    cpl_table_erase_column(zavtab, "p");
    cpl_table_erase_column(zavtab, "av");

    return CPL_ERROR_NONE;
}


cpl_error_code sm_scat_calcscatcorr(cpl_table *zavtab, cpl_table *rhoavtab)
{
    /*!
     * Calculates multiple scattering intensities to correct single and double
     * scattering calculations. For the former, correction factors are
     * derived, which depend on the angle between line of sight and Moon
     * \f$\rho\f$ and wavelength. For the latter, the wavelength-dependent
     * scattering intensities for high-order scattering (triple scattering and
     * higher) are computed. The input tables are modified to contain the
     * results of these calculations.
     *
     * \b INPUT:
     * \param zavtab    CPL table with scattering intensities weighted for
     *                  line-of-sight and Moon zenith distances
     * \param rhoavtab  CPL table with weighted scattering intensities
     *                  depending on wavelength
     *
     * \b OUTPUT:
     * \param zavtab    CPL table with multiple scattering corrections for
     *                  single scattering
     * \param rhoavtab  CPL table with multiple scattering intensities
     *                  depending on wavelength
     *
     * \b ERRORS:
     * - none
     */

    int nr = 0, i = 0, nz = 0, j = 0;
    double *lam, *i1s, *i2s, *ixs, *lamz, *ixsz;
    double r21 = 0., rx2 = 0.;

    /* Maximum ratio of double and single scattering intensity */
    double maxr21 = 0.9;

    /* Get number of wavelengths in rho-averaged table */
    nr = cpl_table_get_nrow(rhoavtab);

    /* Create column for multiple scattering intensity in rho-averaged
       table */
    cpl_table_new_column(rhoavtab, "I_xs", CPL_TYPE_DOUBLE);

    /* Initialise column for multiple scattering intensity in rho-averaged
       table */
    cpl_table_fill_column_window(rhoavtab, "I_xs", 0, nr, 0.);

    /* Get pointers to columns for wavelength, single, double, and multiple
       scattering intensities in rho-averaged table */
    lam = cpl_table_get_data_double(rhoavtab, "lam");
    i1s = cpl_table_get_data_double(rhoavtab, "I_1s");
    i2s = cpl_table_get_data_double(rhoavtab, "I_2s");
    ixs = cpl_table_get_data_double(rhoavtab, "I_xs");

    /* Estimate multiple scattering intensity for each wavelength */

    for (i = 0; i < nr; i++) {

        /* Calculate ratio of double and single scattering intensity */
        if (i1s[i] != 0.) {
            r21 = i2s[i] / i1s[i];
        } else {
            r21 = 0.;
        }

        /* Avoid diverging multiple scattering intensities */
        if (r21 > maxr21) {
            r21 = maxr21;
        }

        /* Estimate multiple scattering contribution by assuming a fixed
           ratio for intensities of consecutive scattering orders */
        rx2 = (1. / fabs(1. - r21)) - 1.;

        /* Get multiple scattering intensity */
        ixs[i] = i2s[i] * rx2;

    }

    /* Get number of wavelengths in z-averaged table */
    nz = cpl_table_get_nrow(zavtab);

    /* Create column for multiple scattering intensity in rho-averaged
       table */
    cpl_table_new_column(zavtab, "I_xs", CPL_TYPE_DOUBLE);

    /* Initialise column for multiple scattering intensity in z-averaged
       table */
    cpl_table_fill_column_window(zavtab, "I_xs", 0, nz, 0.);

    /* Get pointers to columns for wavelength and multiple scattering
       intensity in z-averaged table */
    lamz = cpl_table_get_data_double(zavtab, "lam");
    ixsz = cpl_table_get_data_double(zavtab, "I_xs");

    /* Copy multiple scattering intensities to z-averaged table */

    for (i = 0, j = 0; j < nz; j++) {

        /* Get correct index of wavelength in rho-averaged table */
        while (lam[i] != lamz[j]) {
            i++;
            if (i == nr) {
                i = 0;
            }
        }

        /* Copy multiple scattering intensity from rho-averaged to
           z-averaged table */
        ixsz[j] = ixs[i];

    }

    /* Calculate total scattering intensities and provide them in new column
       of z-averaged table */
    cpl_table_duplicate_column(zavtab, "I", zavtab, "I_u2s");
    cpl_table_add_columns(zavtab, "I", "I_xs");

    /* Calculate multiple scattering correction factors for single scattering
       and provide them in new column of z-averaged table */
    cpl_table_duplicate_column(zavtab, "mscorr", zavtab, "I");
    cpl_table_divide_columns(zavtab, "mscorr", "I_1s");

    return CPL_ERROR_NONE;
}


cpl_error_code sm_scat_writesinglescatcorr(const char *sscatcorfile,
                                           const cpl_table *zavtab,
                                           const cpl_table *rhotab,
                                           const smspec *lamgrid)
{
    /*!
     * Writes ASCII file of given name with multiple scattering correction
     * factors for single scattering calculations. The provided grid of
     * factors depends on the angle between line of sight and Moon \f$\rho\f$
     * and wavelength. The output file has a format that can be read as
     * ::smgrid structure.
     *
     * \b INPUT:
     * \param sscatcorname  path and name of file with multiple scattering
     *                      correction factors for single scattering
     *                      calculations
     * \param zavtab        CPL table with multiple scattering correction
     *                      factors
     * \param rhotab        grid of angles between line of sight and Moon
     * \param lamgrid       grid of wavelengths
     *
     * \b OUTPUT:
     * - none
     *
     * \b ERRORS:
     * - IOS: Invalid object structure
     * - IDG: Inconsistent data grids
     * - FOF: File opening failed
     */

    FILE *stream;
    char errtxt[SM_MAXLEN+1], lamstr[SM_MAXLEN+1], substr[SM_LENLINE+1];
    char rhostr[SM_MAXLEN+1], corrstr[SM_MAXLEN+1];
    cpl_boolean isminrho = CPL_TRUE, ismaxrho = CPL_TRUE;
    cpl_boolean isminlam = CPL_TRUE, ismaxlam = CPL_TRUE;
    int n = 0, nr = 0, nl = 0, i = 0, j = 0, k = 0, nrm = 0, nlm = 0;
    const double *rho = NULL, *lam = NULL, *corr = NULL, *rho0 = NULL;
    double minrho = 0., maxrho = 180., minlam = 0.3, maxlam = 30.;

    /* Get number of rows in data table */
    n = cpl_table_get_nrow(zavtab);

    /* Get pointers to required columns in data table */
    rho = cpl_table_get_data_double_const(zavtab, "rho");
    lam = cpl_table_get_data_double_const(zavtab, "lam");
    corr = cpl_table_get_data_double_const(zavtab, "mscorr");

    /* Get number of rho angles and wavelengths from input grid table */
    nr = cpl_table_get_nrow(rhotab);
    nl = lamgrid->n;

    /* Get pointer to grid column in rho table */
    rho0 = cpl_table_get_data_double_const(rhotab, "rho");

    /* Check agreement of sizes of rho and wavelength grids */
    if (nr * nl != n) {
        sprintf(errtxt, "%s: cpl_table *zavtab (size != product of rhotab "
                "and lamgrid sizes", SM_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IOS, "%s", errtxt);

    }

    /* Check agreement of rho and wavelength grids */
    for (i = 0, j = 0; j < nr; j++) {
        for (k = 0; k < nl; k++) {
            if (rho[i] != rho0[j] || lam[i] != lamgrid->dat[k].lam) {
                sprintf(errtxt, "%s: cpl_table *zavtab (lam and rho) != "
                        "cpl_table *rhotab (rho) and smspec *lamgrid (lam)",
                        SM_ERROR_IDG_TXT);
                return cpl_error_set_message(cpl_func, SM_ERROR_IDG, "%s",
                                             errtxt);
            }
            i++;
        }
    }

    /* Open output file */
    if ((stream = fopen(sscatcorfile, "w+")) == NULL) {
        sprintf(errtxt, "%s: %s", SM_ERROR_FOF_TXT, sscatcorfile);
        return cpl_error_set_message(cpl_func, SM_ERROR_FOF, "%s", errtxt);
    }

    /* Write file header */
    fprintf(stream, "# Multiple Rayleigh and Mie scattering correction "
            "factors:\n");

    /* Check whether required minimum and maximum values exist and increase
       number of rho or wavelengths if required */
    nrm = nr;
    nlm = nl;
    if (rho0[0] > minrho) {
        isminrho = CPL_FALSE;
        nrm++;
    }
    if (rho0[nr-1] < maxrho) {
        ismaxrho = CPL_FALSE;
        nrm++;
    }
    if (lamgrid->dat[0].lam > minlam) {
        isminlam = CPL_FALSE;
        nlm++;
    }
    if (lamgrid->dat[nl-1].lam < maxlam) {
        ismaxlam = CPL_FALSE;
        nlm++;
    }

    /* Write number of rho and wavelengths into file */
    fprintf(stream, "# Number of grid points (lambda & rho)\n");
    fprintf(stream, "%d %d\n", nlm, nrm);

    /* Get string of wavelengths */
    if (isminlam == CPL_FALSE) {
        sprintf(lamstr, "%g %g", minlam, lamgrid->dat[0].lam);
    } else {
        sprintf(lamstr, "%g", lamgrid->dat[0].lam);
    }
    for (k = 1; k < nl; k++) {
        sprintf(substr, " %g", lamgrid->dat[k].lam);
        strcat(lamstr, substr);
    }
    if (ismaxlam == CPL_FALSE) {
        sprintf(substr, " %g", maxlam);
        strcat(lamstr, substr);
    }

    /* Write grid of wavelengths */
    fprintf(stream, "# Wavelength [micron]\n");
    fprintf(stream, "%s\n", lamstr);

    /* Get string of rho */
    if (isminrho == CPL_FALSE) {
        sprintf(rhostr, "%g %g", minrho, rho0[0]);
    } else {
        sprintf(rhostr, "%g", rho0[0]);
    }
    for (j = 1; j < nr; j++) {
        sprintf(substr, " %g", rho0[j]);
        strcat(rhostr, substr);
    }
    if (ismaxrho == CPL_FALSE) {
        sprintf(substr, " %g", maxrho);
        strcat(rhostr, substr);
    }

    /* Write grid of rho */
    fprintf(stream, "# Angle between line of sight and radiation source "
            "rho [deg]\n");
    fprintf(stream, "%s\n", rhostr);

    /* Write comment on grid of multiple scattering correction factors */
    fprintf(stream, "# Correction factors from double scattering "
            "calculations\n");

    /* Write strings of rho-dependent correction factors for each
       wavelength */

    for (k = 0; k < nl; k++) {

        /* Get string for given wavelength index */
        if (isminrho == CPL_FALSE) {
            sprintf(corrstr, "%5.3f %5.3f", corr[k], corr[k]);
        } else {
            sprintf(corrstr, "%5.3f", corr[k]);
        }
        for (j = 1; j < nr; j++) {
            sprintf(substr, " %5.3f", corr[nl*j+k]);
            strcat(corrstr, substr);
        }
        if (ismaxrho == CPL_FALSE) {
            sprintf(substr, " %5.3f", corr[nl*(nr-1)+k]);
            strcat(corrstr, substr);
        }

        /* Write string (twice in the case of lacking data for the extreme
           wavelength limits */
        fprintf(stream, "%s\n", corrstr);
        if (k == 0 && isminlam == CPL_FALSE) {
             fprintf(stream, "%s\n", corrstr);
        }
        if (k == nl-1 && ismaxlam == CPL_FALSE) {
             fprintf(stream, "%s\n", corrstr);
        }

    }

    /* Close file */
    fclose(stream);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_scat_writedoublescatcorr(const char *dscatcorfile,
                                           const cpl_table *rhoavtab)
{
    /*!
     * Writes ASCII file of given name with wavelength-dependent multiple
     * scattering correction intensities for double scattering calculations.
     * The output file has a format that can be read as ::smspec structure.
     *
     * \b INPUT:
     * \param dscatcorname  path and name of file with multiple scattering
     *                      correction intensities for double scattering
     *                      calculations
     * \param rhoavtab      CPL table with multiple scattering correction
     *                      intensities
     *
     * \b OUTPUT:
     * - none
     *
     * \b ERRORS:
     * - FOF: File opening failed
     */

    FILE *stream;
    char errtxt[SM_MAXLEN+1];
    cpl_boolean isminlam = CPL_TRUE, ismaxlam = CPL_TRUE;
    int n = 0, i = 0;
    const double *lam = NULL, *ixs = NULL;
    double minlam = 0.3, maxlam = 30.;

    /* Get number of rows in data table */
    n = cpl_table_get_nrow(rhoavtab);

    /* Get pointers to required columns in data table */
    lam = cpl_table_get_data_double_const(rhoavtab, "lam");
    ixs = cpl_table_get_data_double_const(rhoavtab, "I_xs");

    /* Open output file */
    if ((stream = fopen(dscatcorfile, "w+")) == NULL) {
        sprintf(errtxt, "%s: %s", SM_ERROR_FOF_TXT, dscatcorfile);
        return cpl_error_set_message(cpl_func, SM_ERROR_FOF, "%s", errtxt);
    }

    /* Write file header */
    fprintf(stream, "# Multiple Rayleigh and Mie scattering correction "
            "intensities:\n");
    fprintf(stream, "# 1st column: wavelength [micron]\n");
    fprintf(stream, "# 2nd column: multiple scattering intensity I_xs "
            "[arcsec^-2]\n");

    /* Check whether required minimum and maximum wavelengths exist */
    if (lam[0] > minlam) {
        isminlam = CPL_FALSE;
    }
    if (lam[n-1] < maxlam) {
        ismaxlam = CPL_FALSE;
    }

    /* Write wavelengths and multiple scattering intensities */
    if (isminlam == CPL_FALSE) {
        fprintf(stream, "%6.3f %9.3e\n", minlam, ixs[0]);
    }
    for (i = 0; i < n ; i++) {
        fprintf(stream, "%6.3f %9.3e\n", lam[i], ixs[i]);
    }
    if (ismaxlam == CPL_FALSE) {
        fprintf(stream, "%6.3f %9.3e\n", maxlam, ixs[n-1]);
    }

    /* Close file */
    fclose(stream);

    return CPL_ERROR_NONE;
}

/**@}*/
