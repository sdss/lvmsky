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
 * \ingroup sm_module2_run
 */

/**@{*/

/*!
 * \file sm_run.c
 *
 * Routines for the ETC sky model stand-alone programme
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  27 May 2010
 * \date   06 Oct 2015
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <sm_run.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

cpl_error_code sm_run_readmodelpar(cpl_parameterlist *params,
                                   const char *parfile)
{
    /*!
     * Read parameters of ETC sky model parameter file and put them in a
     * CPL parameter list.
     *
     * \b INPUT:
     * \param parfile  input parameter ASCII file
     *
     * \b OUTPUT:
     * \param params   CPL parameter list containing the read parameters
     *
     * \b ERRORS:
     * - FOF: File opening failed
     * - see ::sm_param_readcheck as well
     */

    FILE *stream;
    cpl_parameter *p;
    smparam x[SM_MAXPAR];
    char errtxt[SM_MAXLEN+1];

    /* Check existence of sky emission parameter file */

    if ((stream = fopen(parfile, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", SM_ERROR_FOF_TXT, parfile);
        return cpl_error_set_message(cpl_func, SM_ERROR_FOF, "%s", errtxt);
    }

    /* Read sky emission parameter file */

    sm_param_readcheck(stream, x, "sm_h", 1);
    p = cpl_parameter_new_value("sm_h", CPL_TYPE_DOUBLE, "", "", x[1].d);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "sm_hmin", 1);
    p = cpl_parameter_new_value("sm_hmin", CPL_TYPE_DOUBLE, "", "", x[1].d);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "alt", 1);
    p = cpl_parameter_new_value("alt", CPL_TYPE_DOUBLE, "", "", x[1].d);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "alpha", 1);
    p = cpl_parameter_new_value("alpha", CPL_TYPE_DOUBLE, "", "", x[1].d);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "rho", 1);
    p = cpl_parameter_new_value("rho", CPL_TYPE_DOUBLE, "", "", x[1].d);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "altmoon", 1);
    p = cpl_parameter_new_value("altmoon", CPL_TYPE_DOUBLE, "", "", x[1].d);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "moondist", 1);
    p = cpl_parameter_new_value("moondist", CPL_TYPE_DOUBLE, "", "", x[1].d);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "pres", 1);
    p = cpl_parameter_new_value("pres", CPL_TYPE_DOUBLE, "", "", x[1].d);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "ssa", 1);
    p = cpl_parameter_new_value("ssa", CPL_TYPE_DOUBLE, "", "", x[1].d);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "calcds", 1);
    p = cpl_parameter_new_value("calcds", CPL_TYPE_STRING, "", "", x[1].c);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "o3column", 1);
    p = cpl_parameter_new_value("o3column", CPL_TYPE_DOUBLE, "", "", x[1].d);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "moonscal", 1);
    p = cpl_parameter_new_value("moonscal", CPL_TYPE_DOUBLE, "", "", x[1].d);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "lon_ecl", 1);
    p = cpl_parameter_new_value("lon_ecl", CPL_TYPE_DOUBLE, "", "", x[1].d);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "lat_ecl", 1);
    p = cpl_parameter_new_value("lat_ecl", CPL_TYPE_DOUBLE, "", "", x[1].d);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "emis_str", 1);
    p = cpl_parameter_new_value("emis_str", CPL_TYPE_STRING, "", "", x[1].c);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "temp_str", 1);
    p = cpl_parameter_new_value("temp_str", CPL_TYPE_STRING, "", "", x[1].c);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "msolflux", 1);
    p = cpl_parameter_new_value("msolflux", CPL_TYPE_DOUBLE, "", "", x[1].d);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "season", 1);
    p = cpl_parameter_new_value("season", CPL_TYPE_INT, "", "", x[1].i);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "time", 1);
    p = cpl_parameter_new_value("time", CPL_TYPE_INT, "", "", x[1].i);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "vac_air", 1);
    p = cpl_parameter_new_value("vac_air", CPL_TYPE_STRING, "", "", x[1].c);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "pwv", 1);
    p = cpl_parameter_new_value("pwv", CPL_TYPE_DOUBLE, "", "", x[1].d);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "rtcode", 1);
    p = cpl_parameter_new_value("rtcode", CPL_TYPE_STRING, "", "", x[1].c);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "resol", 1);
    p = cpl_parameter_new_value("resol", CPL_TYPE_DOUBLE, "", "", x[1].d);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "filepath", 1);
    p = cpl_parameter_new_value("filepath", CPL_TYPE_STRING, "", "", x[1].c);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "incl", 1);
    p = cpl_parameter_new_value("incl", CPL_TYPE_STRING, "", "", x[1].c);
    cpl_parameterlist_append(params, p);

    fclose(stream);

    return cpl_error_get_code();
}


cpl_error_code sm_run_readinstrupar(cpl_parameterlist *lamgrid,
                                    cpl_parameterlist *params,
                                    const char *parfile)
{
    /*!
     * Read parameters of instrument parameter file and put them in the CPL
     * parameter lists \e lamgrid and \e params. In the latter case, the read
     * parameters are appended to a list already filled by
     * ::sm_run_readmodelpar.
     *
     * \b INPUT:
     * \param params    CPL parameter list containing the sky model parameters
     * \param parfile   input parameter ASCII file
     *
     * \b OUTPUT:
     * \param lamgrid   CPL parameter list containing the wavelength grid
     *                  parameters
     * \param params    CPL parameter list containing the sky model and the
     *                  convolution parameters
     *
     * \b ERRORS:
     * - FOF: File opening failed
     * - see ::sm_param_readcheck as well
     */

    FILE *stream;
    cpl_parameter *p;
    smparam x[SM_MAXPAR];
    char errtxt[SM_MAXLEN+1];

    /* Check existence of instrument parameter file */

    if ((stream = fopen(parfile, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", SM_ERROR_FOF_TXT, parfile);
        return cpl_error_set_message(cpl_func, SM_ERROR_FOF, "%s", errtxt);
    }

    /* Read instrument parameter file */

    sm_param_readcheck(stream, x, "limlam", 2);
    p = cpl_parameter_new_value("minlam", CPL_TYPE_DOUBLE, "", "", x[1].d);
    cpl_parameterlist_append(lamgrid, p);
    p = cpl_parameter_new_value("maxlam", CPL_TYPE_DOUBLE, "", "", x[2].d);
    cpl_parameterlist_append(lamgrid, p);

    sm_param_readcheck(stream, x, "dlam", 1);
    p = cpl_parameter_new_value("dlam", CPL_TYPE_DOUBLE, "", "", x[1].d);
    cpl_parameterlist_append(lamgrid, p);

    sm_param_readcheck(stream, x, "kernrad", 1);
    p = cpl_parameter_new_value("kernrad", CPL_TYPE_INT, "", "", x[1].i);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "wbox", 1);
    p = cpl_parameter_new_value("wbox", CPL_TYPE_DOUBLE, "", "", x[1].d);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "wgauss", 1);
    p = cpl_parameter_new_value("wgauss", CPL_TYPE_DOUBLE, "", "", x[1].d);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "wlorentz", 1);
    p = cpl_parameter_new_value("wlorentz", CPL_TYPE_DOUBLE, "", "", x[1].d);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "varkern", 1);
    p = cpl_parameter_new_range("varkern", CPL_TYPE_INT, "", "", x[1].i, 0,
                                1);
    cpl_parameterlist_append(params, p);

    sm_param_readcheck(stream, x, "kernelfile", 1);
    p = cpl_parameter_new_value("kernelfile", CPL_TYPE_STRING, "", "",
                                x[1].c);
    cpl_parameterlist_append(params, p);

    fclose(stream);

    /* Add additional parameter */

    p = cpl_parameter_new_value("kernscal", CPL_TYPE_DOUBLE, "", "", 1.);
    cpl_parameterlist_append(params, p);

    return cpl_error_get_code();
}


cpl_error_code sm_run_createtable(cpl_table *skytable,
                                  const cpl_parameterlist *lamgrid)
{
    /*!
     * Creates CPL table with wavelength column from wavelength range and
     * step provided by the CPL parameter list \e lamgrid.
     *
     * \b INPUT:
     * \param lamgrid   CPL parameter list containing the wavelength grid
     *                  parameters
     *
     * \b OUTPUT:
     * \param skytable  wavelength grid as CPL table
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     */

    const cpl_parameter *p;
    char errtxt[SM_MAXLEN+1];
    int nrow = 0, i = 0;
    double minlam = 0., maxlam = 0., dlam = 0., lam = 0.;

    /* Get wavelength grid parameters from CPL parameter list */

    p = cpl_parameterlist_find_const(lamgrid, "minlam");
    minlam = cpl_parameter_get_double(p);
    p = cpl_parameterlist_find_const(lamgrid, "maxlam");
    maxlam = cpl_parameter_get_double(p);
    p = cpl_parameterlist_find_const(lamgrid, "dlam");
    dlam = cpl_parameter_get_double(p);

    /* Derive number of data points */

    if (maxlam <= minlam || dlam <= 0.) {
        nrow = 0;
    } else {
        nrow = ((maxlam - minlam) / dlam + 0.5) + 1;
    }

    /* Create table */

    cpl_table_set_size(skytable, nrow);
    cpl_table_new_column(skytable, "lam", CPL_TYPE_DOUBLE);

    /* Return table with zero data points if parameters are invalid */

    if (nrow == 0) {
        sprintf(errtxt, "%s: parameters of cpl_parameterlist *lamgrid"
                "(maxlam <= minlam || dlam <= 0.)", SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Fill table with desired wavelength grid */

    lam = minlam;
    cpl_table_set_double(skytable, "lam", 0, lam);
    for (i = 1; i < nrow; i++) {
        lam += dlam;
        cpl_table_set_double(skytable, "lam", i, lam);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_run_calcmodel(cpl_table *skytable,
                                cpl_parameterlist *params)
{
    /*!
     * \callgraph
     * Reads wavelength grid as CPL table and sky model and convolution
     * parameters as CPL parameter list and uses this information to compute a
     * corresponding radiance spectrum and transmission curve convolved to the
     * resolution given. The output CPL table includes the input wavelength
     * grid and the computed flux/transmission, errors, and individual
     * components. The main output columns are "lam", "flux", "dflux1",
     * "dflux2", "trans", "dtrans1", "dtrans2". Moreover, four transmission
     * and seven radiance columns containing the model components are added.
     *
     * \b INPUT:
     * \param skytable  input wavelength grid as CPL table
     * \param params    CPL parameter list containing the sky model and the
     *                  convolution parameters
     *
     * \b OUTPUT:
     * \param skytable  CPL table with full radiance and transmission data
     *
     * \b ERRORS:
     * - none
     */

    cpl_error_code status = CPL_ERROR_NONE;

    /* Calculate sky brightness and transmission curve */
    sm_etc_calcmodel(skytable, params);

    /* Return if error(s) occurred */
    if ((status = cpl_error_get_code()) != CPL_ERROR_NONE) {
        return status;
    }

    /* Convolve output spectra with given kernel */
    sm_run_convolvespec(skytable, params);

    /* Return error code of last error */
    return cpl_error_get_code();
}


cpl_error_code sm_run_convolvespec(cpl_table *skytable,
                                   cpl_parameterlist *params)
{
    /*!
     * \callgraph
     * Convolves a table with sky emission and transmission data with a
     * wavelength-dependent kernel consisting of boxcar, Gaussian, and
     * Lorentzian of given widths. The scaling of the kernel width is forced
     * by the assumption of constant resolution. If a constant kernel is
     * desired, the \e varkern parameter has to be set to 0. For a
     * wavelength-dependent kernel, a value of 1 is required.
     *
     * \b INPUT:
     * \param skytable  radiance and transmission data (CPL table)
     * \param params    CPL parameter list containing the sky model and the
     *                  convolution parameters
     *
     * \b OUTPUT:
     * \param skytable  convolved radiance and transmission data
     *
     * \b ERRORS:
     * - NDA: No data
     * - IOS: Invalid object structure
     */

    cpl_parameter *p;
    cpl_array *col = NULL, *kernel = NULL;
    cpl_table *temptable = NULL;
    char errtxt[SM_MAXLEN+1];
    int nrow = 0, ncol = 0, h = 0, varkern = 0, range[2] = {0, 0}, i = 0;
    double limlam[2] = {0., 0.}, reflam = 0., speedpar = 0., llam = 0.;
    double ulam = 0., clam = 0., kernscal = 0.;
    double *lam;

    /* Check number of data points in spectrum */
    nrow = cpl_table_get_nrow(skytable);
    if (nrow <= 0) {
        sprintf(errtxt, "%s: cpl_table *skytable", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    }

    /* Check number of table columns */
    ncol = cpl_table_get_ncol(skytable);
    if (ncol < 7) {
        sprintf(errtxt, "%s: cpl_table *skytable (# of columns < 7)",
                SM_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IOS, "%s", errtxt);
    }

    /* Check existence of 'lam' column */
    if (cpl_table_has_column(skytable, "lam") != 1) {
        sprintf(errtxt, "%s: cpl_table *spec (no 'lam' column)",
                SM_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IOS, "%s", errtxt);
    }

    /* Get table column names */
    col = cpl_table_get_column_names(skytable);

    /* Check position of wavelength column in table */
    if (strcmp(cpl_array_get_string(col, 0), "lam") != 0) {
        cpl_array_delete(col);
        sprintf(errtxt, "%s: cpl_table *spec ('lam' column != first column)",
                SM_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IOS, "%s", errtxt);
    }

    /* Get pointer to wavelength column */
    lam = cpl_table_get_data_double(skytable, "lam");

    /* Create temporary CPL table */
    temptable = cpl_table_duplicate(skytable);

    /* Set values to 0 in initial table (exception: wavelength) */
    for (h = 1; h < ncol; h++) {
        cpl_table_fill_column_window(skytable, cpl_array_get_string(col, h),
                                     0, nrow, 0.);
    }

    /* Get central wavelength */
    limlam[0] = cpl_table_get(skytable, "lam", 0, NULL);
    limlam[1] = cpl_table_get(skytable, "lam", nrow - 1, NULL);
    reflam = (limlam[0] + limlam[1]) / 2;

    /* Kernel: linearly increasing with wavelength or constant? */
    p = cpl_parameterlist_find(params, "varkern");
    varkern = cpl_parameter_get_int(p);
    if (varkern == 1) {
        /* variable kernel */
        speedpar = SM_LIMRELLAMVAR;
    } else {
        /* constant kernel */
        speedpar = HUGE_VAL;
    }

    /* CPL array for kernel */
    kernel = cpl_array_new(0, CPL_TYPE_DOUBLE);

    /* Convolve spectral ranges with kernel */

    do {

        /* Get lower and upper wavelengths of range */
        llam = lam[range[0]];
        ulam = llam * (1 + speedpar);

        /* Find upper range pixel */
        i = range[0];
        while (i < nrow && lam[i] < ulam) {
            i++;
        }
        if (i == nrow) {
            range[1] = nrow - 1;
        } else {
            range[1] = i;
        }

        /* Get central wavelength of range */
        ulam = lam[range[1]];
        clam = (ulam + llam) / 2;

        /* Derive FWHM scaling factor */
        if (speedpar == HUGE_VAL) {
            kernscal = 1.;
        } else {
            kernscal =  clam / reflam;
        }

        /* Write FWHM scaling factor into parameter list */
        p = cpl_parameterlist_find(params, "kernscal");
        cpl_parameter_set_double(p, kernscal);

        /* Calculate kernel */
        sm_run_calckernel(kernel, params);

        /* Write kernel for central wavelength to file or on stdout */
        if (llam < reflam && ulam > reflam) {
            sm_run_writekernel(kernel, params);
        }

        /* Convolve spectrum with combined kernel */
        sm_run_convolvewindow(skytable, temptable, range, kernel);

        /* Set next lower range pixel */
        range[0] = range[1] + 1;

    } while (range[0] < nrow);

    /* Free memory */
    cpl_array_delete(col);
    cpl_array_delete(kernel);
    cpl_table_delete(temptable);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_run_calckernel(cpl_array *kernel,
                                 const cpl_parameterlist *params)
{
    /*!
     * \callgraph
     * Calculates combined boxcar, Gaussian, and Lorentzian kernel.
     *
     * \b INPUT:
     * \param params  CPL parameter list containing the sky model and the
     *                convolution parameters
     *
     * \b OUTPUT:
     * \param kernel  CPL array containing the kernel elements
     *
     * \b ERRORS:
     * - IOV: Invalid object value(s)
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_array *voigtkern, *boxkern;
    char errtxt[SM_MAXLEN+1];
    int margpix = 0, nkpix = 0, k = 0;
    double sum = 0.;
    double *vkern, *kern;

    /* Calculate Voigt profile kernel */
    voigtkern = cpl_array_new(0, CPL_TYPE_DOUBLE);
    status = sm_run_calcvoigtkernel(voigtkern, params);
    if (status != CPL_ERROR_NONE) {
        cpl_array_delete(voigtkern);
        return status;
    }

    /* Calculate boxcar kernel */
    boxkern = cpl_array_new(0, CPL_TYPE_DOUBLE);
    status = sm_run_calcboxkernel(boxkern, params);
    if (status != CPL_ERROR_NONE) {
        cpl_array_delete(voigtkern);
        cpl_array_delete(boxkern);
        return status;
    }

    /* Convolve Voigt profile kernel with boxcar kernel */
    sm_run_convolvearray(voigtkern, boxkern);

    /* Get number of marginal pixels to be cut */
    margpix = 0.5 * (cpl_array_get_size(boxkern) - 1);

    /* Set size of output kernel pixels */
    nkpix = cpl_array_get_size(voigtkern) - 2 * margpix;
    cpl_array_set_size(kernel, nkpix);
    cpl_array_fill_window_double(kernel, 0, nkpix, 0.);

    /* Get pointers to CPL array */
    vkern = cpl_array_get_data_double(voigtkern);
    kern = cpl_array_get_data_double(kernel);

    /* Transfer kernel values */
    for (k = 0; k < nkpix; k++) {
        kern[k] = vkern[k + margpix];
    }

    /* Free memory */
    cpl_array_delete(voigtkern);
    cpl_array_delete(boxkern);

    /* Add all kernel values */
    for (sum = 0, k = 0; k < nkpix; k++) {
        if (kern[k] < 0.) {
            kern[k] = 0.;
        }
        sum += kern[k];
    }

    /* Normalise kernel values */
    if (sum > 0) {
        for (k = 0; k < nkpix; k++) {
            kern[k] /= sum;
        }
    } else {
        sprintf(errtxt, "%s: cpl_array *kernel (sum of kernel elements <= 0)",
                SM_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IOV, "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_run_calcvoigtkernel(cpl_array *kernel,
                                      const cpl_parameterlist *params)
{
    /*!
     * Calculates kernel based on a Voigt profile approximation that depends
     * on the FWHM of the Gaussian and Lorentzian components in pixels. The
     * sum of the kernel values is normalised to 1. The number of kernel
     * pixels is taken from the input kernel radius plus the FWHM of the
     * boxcar that will be convolved with the resulting Voigt profile
     * (see ::sm_run_calckernel). If a kernel corresponding to constant
     * resolution is requested, size and FWHM of the kernel are multiplied by
     * the factor \e kernscal.
     *
     * \b INPUT:
     * \param params  CPL parameter list containing the sky model and the
     *                convolution parameters
     *
     * \b OUTPUT:
     * \param kernel    CPL array containing the kernel elements
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     * - IOV: Invalid object value(s)
     */

    const cpl_parameter *p;
    char errtxt[SM_MAXLEN+1];
    int kernrad = 0, nkpix = 0, refpix = 0, nx = 0, k = 0, i = 0;
    double *kern;
    double kernscal = 0., wbox = 0., wgauss = 0., wlorentz = 0., gamma = 0.;
    double wvoigt = 0., wlwv = 0., xmax = 0., xmin = 0., dx = 0., x = 0.;
    double xv = 0., xv2 = 0., xv225 = 0., sum = 0;

    /* Get FWHM scaling factor */
    p = cpl_parameterlist_find_const(params, "kernscal");
    kernscal = cpl_parameter_get_double(p);
    if (kernscal <= 0) {
        sprintf(errtxt, "%s: kernscal <= 0", SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Get kernel radius */
    p = cpl_parameterlist_find_const(params, "kernrad");
    kernrad = floor((cpl_parameter_get_int(p) + 0.5) * kernscal);
    if (kernrad < 0) {
        sprintf(errtxt, "%s: kernrad < 0", SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Kernel with one pixel only */
    if (kernrad == 0) {
        cpl_array_set_size(kernel, 1);
        cpl_array_set(kernel, 0, 1.);
        return CPL_ERROR_NONE;
    }

    /* Get FWHM of boxcar */
    p = cpl_parameterlist_find_const(params, "wbox");
    wbox = cpl_parameter_get_double(p) * kernscal;
    if (wbox < 0) {
        sprintf(errtxt, "%s: wbox < 0", SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Get FWHM of Gaussian */
    p = cpl_parameterlist_find_const(params, "wgauss");
    wgauss = cpl_parameter_get_double(p) * kernscal;
    if (wgauss < 0) {
        sprintf(errtxt, "%s: wgauss < 0", SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Get FWHM of Lorentzian */
    p = cpl_parameterlist_find_const(params, "wlorentz");
    wlorentz = cpl_parameter_get_double(p) * kernscal;
    if (wlorentz < 0) {
        sprintf(errtxt, "%s: wlorentz < 0", SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Calculate number of kernel pixels */
    nkpix = 2 * kernrad + 1 + 2 * ceil(0.5 * (wbox - 1));
    cpl_array_set_size(kernel, nkpix);

    /* Reference pixel for mirroring of kernel values */
    refpix = (nkpix - 1) / 2;

    /* Case: FWHM of Gaussian and Lorentzian equal zero */
    if (wgauss == 0 && wlorentz == 0) {
        cpl_array_fill_window_double(kernel, 0, nkpix, 0.);
        cpl_array_set(kernel, refpix, 1.);
        return CPL_ERROR_NONE;
    }

    /* gamma of Lorentzian */
    gamma = 0.5 * wlorentz;

    /* FWHM of Voigt profile in pixels */
    wvoigt = gamma + sqrt(gamma * gamma + wgauss * wgauss);

    /* Ratio of Lorentzian and Voigt profile FWHM */
    wlwv = wlorentz / wvoigt;

    /* Get pointer to CPL array */
    cpl_array_fill_window_double(kernel, 0, nkpix, 0.);
    kern = cpl_array_get_data_double(kernel);

    /* Integration limits for upper bin */
    xmax = 0.5 * nkpix;
    xmin = xmax - 1;

    /* Number of points per pixel for integration of Voigt profile */
    nx = ceil(SM_BINS_PER_FWHM / wvoigt);

    /* Step in pixels for integration of Voigt profile */
    dx = 1. / nx;

    /* Calculate kernel up to reference pixel */

    for (k = nkpix - 1; k >= refpix; k--) {

        if (xmax <= 0.) {

            /* Skip integration */
            kern[k] = 1.;

        } else {

            /* First point in pixels for integration */
            x = xmin + dx / 2;

            /* Perform integration */

            kern[k] = 0.;

            for (i = 0; i < nx; i++) {

                /* Get variables of approximation formula */
                xv = x / wvoigt;
                xv2 = xv * xv;
                xv225 = pow(fabs(xv), 2.25);

                /* Calculate Voigt approximation for integration point */
                kern[k] += (1. - wlwv) * exp(-2.772 * xv2)
                         + wlwv / (1. + 4. * xv2)
                         + 0.016 * (1. - wlwv) * wlwv
                         * (exp(-0.4 * xv225) - 10. / (10. + xv225));

                /* Get next integration point */
                x += dx;

            }

            kern[k] /= (double) nx;

        }

        /* Shift integration limits for next bin */
        xmax = xmin;
        xmin = xmax - 1;

    }

    /* Mirror right wing of kernel */
    for (k = refpix - 1; k >= 0; k--) {
        kern[k] = kern[nkpix - k - 1];
    }

    /* Add all kernel values */
    for (sum = 0, k = 0; k < nkpix; k++) {
        if (kern[k] < 0.) {
            kern[k] = 0.;
        }
        sum += kern[k];
    }

    /* Normalise kernel values */
    if (sum > 0) {
        for (k = 0; k < nkpix; k++) {
            kern[k] /= sum;
        }
    } else {
        sprintf(errtxt, "%s: cpl_array *kernel (sum of kernel elements <= 0)",
                SM_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IOV, "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_run_calcboxkernel(cpl_array *kernel,
                                    const cpl_parameterlist *params)
{
    /*!
     * Calculates boxcar kernel. The sum of the kernel values is normalised to
     * 1. If a kernel corresponding to constant resolution is requested, the
     * size is multiplied by the factor \e kernscal.
     *
     * \b INPUT:
     * \param params  CPL parameter list containing the sky model and the
     *                convolution parameters
     *
     * \b OUTPUT:
     * \param kernel  CPL array containing the kernel elements
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     */

    const cpl_parameter *p;
    char errtxt[SM_MAXLEN+1];
    int nkpix = 0;
    double kernscal = 0., wbox = 0., boxrad = 0., mval = 0., fval = 0.;
    double k = 0.;

    /* Get FWHM scaling factor */
    p = cpl_parameterlist_find_const(params, "kernscal");
    kernscal = cpl_parameter_get_double(p);
    if (kernscal <= 0) {
        sprintf(errtxt, "%s: kernscal <= 0", SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Get FWHM of boxcar */
    p = cpl_parameterlist_find_const(params, "wbox");
    wbox = cpl_parameter_get_double(p) * kernscal;
    if (wbox < 0) {
        sprintf(errtxt, "%s: wbox < 0", SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Calculate number of kernel pixels */
    boxrad = 0.5 * (wbox - 1);
    nkpix = 2 * ceil(boxrad) + 1;
    cpl_array_set_size(kernel, nkpix);

    /* Kernel with one pixel only */
    if (nkpix == 1) {
        cpl_array_set_size(kernel, 1);
        cpl_array_set(kernel, 0, 1.);
        return CPL_ERROR_NONE;
    }

    /* Set kernel values */

    mval = fmod(boxrad, 1.) / wbox;
    fval = 1. / wbox;

    for (k = 0; k < nkpix; k++) {
        if (mval > 0 && (k == 0 || k == nkpix-1)) {
            cpl_array_set(kernel, k, mval);
        } else {
            cpl_array_set(kernel, k, fval);
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_run_convolvearray(cpl_array *array,
                                    const cpl_array *kernel)
{
    /*!
     * Convolution of array with kernel.
     *
     * \b INPUT:
     * \param array   CPL array
     * \param kernel  kernel as CPL array (required: sum of all values = 1)
     *
     * \b OUTPUT:
     * \param array   array convolved with given kernel
     *
     * \b ERRORS:
     * - IIP: Invalid object value(s)
     */

    cpl_array *tarray;
    char errtxt[SM_MAXLEN+1];
    int n = 0, nkpix = 0, k = 0, kmin = 0, kmax = 0, i = 0, j = 0;
    const double *kern = NULL;
    double *val = NULL, *tval = NULL;
    double sum = 0.;

    /* No data points -> no convolution */
    n = cpl_array_get_size(array);
    if (n <= 0) {
        return CPL_ERROR_NONE;
    }

    /* No kernel data -> no convolution */
    nkpix = cpl_array_get_size(kernel);
    if (nkpix == 0) {
        return CPL_ERROR_NONE;
    }

    /* Get pointer to kernel */
    kern = cpl_array_get_data_double_const(kernel);

    /* Check kernel */

    for (k = 0; k < nkpix; k++) {
        if (kern[k] < 0 || kern[k] > 1) {
            sprintf(errtxt, "%s: cpl_array *kernel "
                    "(kernel element(s) < 0 or > 1)", SM_ERROR_IOV_TXT);
            return cpl_error_set_message(cpl_func, SM_ERROR_IOV, "%s",
                                         errtxt);
        }
        sum += kern[k];
    }

    if (sum < 1 - SM_TOL || sum > 1 + SM_TOL) {
        sprintf(errtxt, "%s: cpl_array *kernel (sum of kernel elements != 1)",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IOV, "%s", errtxt);
    }

    /* Skip convolution if number of kernel pixels is one */
    if (nkpix == 1) {
        return CPL_ERROR_NONE;
    }

    /* Kernel with even or odd pixel number?
       Note: centre of kernel at -0.5 pixels for even pixel number */

    if (nkpix % 2 == 0) {
        kmin = - nkpix / 2;
    } else {
        kmin = - (nkpix - 1) / 2;
    }
    kmax = kmin + nkpix - 1;

    /* Create temporary CPL array */
    tarray = cpl_array_duplicate(array);

    /* Get pointers to arrays */
    val = cpl_array_get_data_double(array);
    tval = cpl_array_get_data_double(tarray);

    /* Convolve array with kernel */

    for (i = 0; i < n; i++) {

        for (val[i] = 0., k = kmin; k <= kmax; k++) {

            j = i + k;

            /* Value of first or last valid pixel for invalid pixels */
            if (j < 0) {
                j = 0;
            } else if (j >= n) {
                j = n - 1;
            }

            val[i] += tval[j] * kern[k - kmin];

        }

    }

    /* Free memory */
    cpl_array_delete(tarray);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_run_writekernel(const cpl_array *kernel,
                                  const cpl_parameterlist *params)
{
    /*!
     * Writes kernel to file or on stdout (indicated by "stdout" as filename;
     * see parameter \e kernelfile of CPL parameter list \e params). In the
     * case of "null", no kernel is written.
     *
     * \b INPUT:
     * \param kernel  kernel for convolution (CPL array)
     * \param params  CPL parameter list containing the sky model and the
     *                convolution parameters
     *
     * \b ERRORS:
     * - FOF: File opening failed
     */

    FILE *stream;
    const cpl_parameter *p;
    char kernelfile[FILENAME_MAX], errtxt[SM_MAXLEN+1];
    int k = 0;

    p = cpl_parameterlist_find_const(params, "kernelfile");
    strncpy(kernelfile, cpl_parameter_get_string(p), SM_LENLINE+2);

    if (strncmp(kernelfile, "stdout", 6) == 0) {
        /* Output on stdout */
        stream = stdout;
    } else if (strncmp(kernelfile, "null", 4) == 0) {
        /* Do not write kernel */
        return CPL_ERROR_NONE;
    } else {
        /* Check file existence */
        if ((stream = fopen(kernelfile, "w")) == NULL) {
            sprintf(errtxt, "%s: %s", SM_ERROR_FOF_TXT, kernelfile);
            return cpl_error_set_message(cpl_func, SM_ERROR_FOF, "%s",
                                         errtxt);
        }
    }

    fprintf(stream, "# Kernel:\n");
    for (k = 0; k < cpl_array_get_size(kernel); k++) {
        fprintf(stream, "%e\n", cpl_array_get(kernel, k, NULL));
    }

    if (strncmp(kernelfile, "stdout", 6) != 0) {
        fclose(stream);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_run_convolvewindow(cpl_table *outtable,
                                     const cpl_table *intable,
                                     const int range[2],
                                     const cpl_array *kernel)
{
    /*!
     * Convolution of non-wavelength columns (all but the first column) in a
     * table with given kernel for a fixed range of pixels. The output table
     * has to exist and have the same structure as the input table. Note that
     * the convolved flux/transmission is added to the intial flux in the
     * output table.
     *
     * \b INPUT:
     * \param outtable  CPL table with same structure as input table
     * \param intable   CPL table with data to be convolved
     * \param range     range of pixels to be considered (minimum and maximum)
     * \param kernel    kernel as CPL array
     *                  (required: sum of all values = 1)
     *
     * \b OUTPUT:
     * \param outtable  output table with added convolved flux/transmission
     *                  from given range
     *
     * \b ERRORS:
     * - IOS: Invalid object structure
     * - IIP: Invalid input parameter(s)
     * - IOV: Invalid object value(s)
     */

    cpl_array *col = NULL;
    char errtxt[SM_MAXLEN+1];
    int n = 0, nkpix = 0, k = 0, kmin = 0, kmax = 0, jmin = 0, jmax = 0;
    int h = 0, j = 0, i = 0;
    const double *kern = NULL, *influx = NULL;
    double *outflux = NULL;
    double sum = 0., in0 = 0., innm1 = 0., in = 0.;

    /* No data points -> no convolution */
    n = cpl_table_get_nrow(intable);
    if (n <= 0) {
        return CPL_ERROR_NONE;
    }

    /* Check correspondence of table structures */
    if (cpl_table_compare_structure(outtable, intable) != 0) {
        sprintf(errtxt, "%s: cpl_table *outtable != cpl_table *intable "
                "(structure)", SM_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IOS, "%s", errtxt);
    }

    /* Check range values */

    if (range[0] > range[1]) {
        sm_run_copytablewindow(outtable, intable, range);
        sprintf(errtxt, "%s: range[2] (min. > max.)", SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    if (range[0] < 0 || range[1] >= n) {
        sm_run_copytablewindow(outtable, intable, range);
        sprintf(errtxt, "%s: range[2] (invalid pixel)", SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* No kernel data -> no convolution */
    nkpix = cpl_array_get_size(kernel);
    if (nkpix == 0) {
        sm_run_copytablewindow(outtable, intable, range);
        return CPL_ERROR_NONE;
    }

    /* Get pointer to kernel array */
    kern = cpl_array_get_data_double_const(kernel);

    /* Check kernel */

    for (sum = 0., k = 0; k < nkpix; k++) {
        if (kern[k] < 0 || kern[k] > 1) {
            sm_run_copytablewindow(outtable, intable, range);
            sprintf(errtxt, "%s: cpl_array *kernel "
                    "(kernel element(s) < 0 or > 1)", SM_ERROR_IOV_TXT);
            return cpl_error_set_message(cpl_func, SM_ERROR_IOV, "%s",
                                         errtxt);
        }
        sum += kern[k];
    }

    if (sum < 1 - SM_TOL || sum > 1 + SM_TOL) {
        sm_run_copytablewindow(outtable, intable, range);
        sprintf(errtxt, "%s: cpl_array *kernel (sum of kernel elements != 1)",
                SM_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IOV, "%s", errtxt);
    }

    /* Skip convolution if number of kernel pixels is one */
    if (nkpix == 1) {
        sm_run_copytablewindow(outtable, intable, range);
        return CPL_ERROR_NONE;
    }

    /* Get table column names */
    col = cpl_table_get_column_names(intable);

    /* Kernel with even or odd pixel number?
       Note: centre of kernel at -0.5 pixels for even pixel number */
    if (nkpix % 2 == 0) {
        kmin = - nkpix / 2;
    } else {
        kmin = - (nkpix - 1) / 2;
    }
    kmax = kmin + nkpix - 1;

    /* Set pixel range (add virtual pixels for marginal ranges) */

    if (range[0] == 0) {
        jmin = -kmax;
    } else {
        jmin = range[0];
    }

    if (range[1] == n-1) {
        jmax = n-1 - kmin;
    } else {
        jmax = range[1];
    }

    /* Convolve each flux/transmission column with kernel
       (skip 1st column -> wavelength) */

    for (h = 1; h < cpl_array_get_size(col); h++) {

        /* Get pointers to selected column in input and output table */
        influx = cpl_table_get_data_double_const(intable,
                                                cpl_array_get_string(col, h));
        outflux = cpl_table_get_data_double(outtable,
                                            cpl_array_get_string(col, h));

        /* Set flux of virtual input pixels */

        in0 = influx[0];
        innm1 = influx[n-1];

        /* Perform convolution for selected column */

        for (j = jmin; j <= jmax; j++) {

            /* Flux of real and virtual input pixels */

            if (j < 0) {
                in = in0;
            } else if (j >= n) {
                in = innm1;
            } else {
                in = influx[j];
            }

            /* Calculate output flux for each kernel element and add it to the
               corresponding output pixel */

            for (k = SM_MAX(kmin, -j); k <= SM_MIN(kmax, n-1 - j); k++) {

                i = j + k;

                outflux[i] += in * kern[k - kmin];

            }

        }

    }

    /* Free memory */
    cpl_array_delete(col);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_run_copytablewindow(cpl_table *outtable,
                                      const cpl_table *intable,
                                      const int range[2])
{
    /*!
     * Copies the content of a table except for the first column which is
     * expected to be the wavelength column. The output table must exist and
     * have the same structure as the input table. Only columns with
     * CPL_TYPE_DOUBLE can be handled.
     *
     * \b INPUT:
     * \param outtab  output CPL table with same structure as input table
     * \param intab   input CPL table
     *
     * \b OUTPUT:
     * \param outtab  output table with copied content of input table
     *
     * \b ERRORS:
     * - IOS: Invalid object structure
     * - IIP: Invalid input parameter(s)
     * - IOV: Invalid object value(s)
     */

    cpl_array *col = NULL;
    char errtxt[SM_MAXLEN+1];
    int n = 0, h = 0, i = 0;
    const double *influx = NULL;
    double *outflux = NULL;

    /* No data points -> no convolution */
    n = cpl_table_get_nrow(intable);
    if (n <= 0) {
        return CPL_ERROR_NONE;
    }

    /* Check correspondence of table structures */
    if (cpl_table_compare_structure(outtable, intable) != 0) {
        sprintf(errtxt, "%s: cpl_table *outtable != cpl_table *intable "
                "(structure)", SM_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IOS, "%s", errtxt);
    }

    /* Check range values */

    if (range[0] > range[1]) {
        sprintf(errtxt, "%s: range[2] (min. > max.)", SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    if (range[0] < 0 || range[1] >= n) {
        sprintf(errtxt, "%s: range[2] (invalid pixel)", SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Get table column names */
    col = cpl_table_get_column_names(intable);

    /* Copy each flux/transmission column (skip 1st column -> wavelength) */

    for (h = 1; h < cpl_array_get_size(col); h++) {

        /* Get pointers to selected column in input and output table */
        influx = cpl_table_get_data_double_const(intable,
                                                cpl_array_get_string(col, h));
        outflux = cpl_table_get_data_double(outtable,
                                            cpl_array_get_string(col, h));

        /* Copy selected column */
        for (i = range[0]; i <= range[1]; i++) {
            outflux[i] = influx[i];
        }

    }

    /* Free memory */
    cpl_array_delete(col);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_run_writefits(const cpl_table *skytable,
                                const char *radfile, const char *transfile)
{
    /*!
     * Writes sky model CPL table to a FITS file for radiance and to a FITS
     * file for transmission.
     *
     * \b INPUT:
     * \param skytable   CPL table with radiance and transmission data
     * \param radfile    name of sky emission spectrum file
     * \param transfile  name of transmission curve file
     *
     * \b ERRORS:
     * - IOS: Invalid object structure
     */

    cpl_array *colnames;
    cpl_table *spec, *trans;
    char errtxt[SM_MAXLEN+1], colname[SM_MAXLEN+1];
    int ncol = 0, j = 0;

    /* Check the existence of the mandatory columns */
    if (cpl_table_get_ncol(skytable) < 7 ||
        cpl_table_has_column(skytable, "lam") != 1 ||
        cpl_table_has_column(skytable, "flux") != 1 ||
        cpl_table_has_column(skytable, "dflux1") != 1 ||
        cpl_table_has_column(skytable, "dflux2") != 1 ||
        cpl_table_has_column(skytable, "trans") != 1 ||
        cpl_table_has_column(skytable, "dtrans1") != 1 ||
        cpl_table_has_column(skytable, "dtrans2") != 1) {
        sprintf(errtxt, "%s: cpl_table *skytable (number of columns < 7 "
                "and/or invalid column names)", SM_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IOS, "%s", errtxt);
    }

    /* Duplicate the input table for sky emission and transmission files */
    spec = cpl_table_duplicate(skytable);
    trans = cpl_table_duplicate(skytable);

    /* Get number and names of spec table columns */
    ncol = cpl_table_get_ncol(spec);
    colnames = cpl_table_get_column_names(spec);

    /* Erase transmission columns in spec table */
    for (j = 0; j < ncol; j++) {
        sprintf(colname, "%s", cpl_array_get_string(colnames, j));
        if (strstr(colname, "trans") != NULL) {
            cpl_table_erase_column(spec, colname);
        }
    }

    /* Free memory */
    cpl_array_delete(colnames);

    /* Get number and names of trans table columns */
    ncol = cpl_table_get_ncol(trans);
    colnames = cpl_table_get_column_names(trans);

    /* Erase flux columns in trans table */
    for (j = 0; j < ncol; j++) {
        sprintf(colname, "%s", cpl_array_get_string(colnames, j));
        if (strstr(colname, "flux") != NULL) {
            cpl_table_erase_column(trans, colname);
        }
    }

    /* Free memory */
    cpl_array_delete(colnames);

    /* Write output FITS files */
    cpl_table_save(spec, NULL, NULL, radfile, CPL_IO_CREATE);
    cpl_table_save(trans, NULL, NULL, transfile, CPL_IO_CREATE);

    /* Free memory */
    cpl_table_delete(spec);
    cpl_table_delete(trans);

    return CPL_ERROR_NONE;
}

/**@}*/
