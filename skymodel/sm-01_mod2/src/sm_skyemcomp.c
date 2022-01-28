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
 * \file sm_skyemcomp.c
 *
 * Routines for the ETC sky model
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  02 Oct 2009
 * \date   06 Oct 2015
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <sm_skyemcomp.h>
#include <sm_scatmoonlight.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

cpl_error_code sm_etc_calcmodel(cpl_table *skytable,
                                const cpl_parameterlist *params)
{
    /*!
     * \callgraph
     * Reads wavelength grid as CPL table and sky model parameters as CPL
     * parameter list and uses this information to compute a corresponding
     * radiance spectrum and transmission curve which are committed to the
     * ETC as a CPL table including the input wavelength grid and the
     * resulting flux/transmission, errors, and individual components.
     *
     * \b INPUT:
     * \param skytable  input wavelength grid as CPL table with column "lam"
     * \param params    CPL sky model parameter list with the parameters
     *                  "sm_h" (D), "sm_hmin" (D),
     *                  "alt" (D), "alpha" (D), "rho" (D), "altmoon" (D),
     *                  "moondist" (D), "pres" (D), "ssa" (D), "calcds" (S),
     *                  "o3column" (D), "moonscal" (D), "lon_ecl" (D),
     *                  "lat_ecl" (D), "emis_str" (S), temp_str (S),
     *                  "msolflux" (D), "season" (I), "time" (I),
     *                  "vac_air" (S), "pwv" (D), "rtcode" (S), "resol" (D),
     *                  "filepath" (S), and "incl" (S)
     *
     * \b OUTPUT:
     * \param skytable  CPL table with full radiance and transmission data in
     *                  the columns "lam", "flux", "dflux1", "dflux2",
     *                  "trans", "dtrans1", and "dtrans2", "trans_ma",
     *                  "trans_o3", "trans_rs", "trans_ms", "flux_sml",
     *                  "flux_ssl", "flux_zl", "flux_tie", "flux_tme",
     *                  "flux_ael", and "flux_arc"
     *
     * \b ERRORS:
     * - see subroutines
     */

    cpl_table *radcomp, *transcomp;
    smparmodel modelpar;
    smspec radiance, transmission;
    int error = 0; /* non-zero upon error */

    /* Fill smparmodel structure with data of the input CPL parameter list
       and a special file containing crucial paths and file names */
    if (sm_etc_getparams(&modelpar, params))
        return cpl_error_set_where(cpl_func);
    if (sm_etc_readfilenames(&modelpar, params))
        return cpl_error_set_where(cpl_func);

    /* Initialise CPL tables for sky model components (radiance and
       transmission) */
    radcomp = cpl_table_new(0);
    transcomp = cpl_table_new(0);

    /* Write input wavelength grids to smspec structures for radiance and
     * transmission
     * Compute sky model
     * Write radiance spectrum and transmission curve to output CPL table
     */
    error = sm_spec_readcpl(&radiance, skytable) ||
            sm_spec_readcpl(&transmission, skytable) ||
            sm_comp_skyemcomp(&radiance, &transmission, radcomp, transcomp,
                              modelpar) ||
            sm_etc_writetable(skytable, &radiance, &transmission, radcomp,
                              transcomp);

    /* Free allocated memory */
    cpl_table_delete(radcomp);
    cpl_table_delete(transcomp);
    sm_spec_free(&radiance);
    sm_spec_free(&transmission);

    /* Return error code of last error if any */
    return error ? cpl_error_set_where(cpl_func) : CPL_ERROR_NONE;
}


cpl_error_code sm_etc_getparams(smparmodel *modelpar,
                                const cpl_parameterlist *params)
{
    /*!
     * Reads sky model parameters from a CPL parameter list and put them in a
     * ::smparmodel structure.
     * The strings of comma-separated emissivity and temperature values are
     * converted into double arrays.
     * As only non-input parameter the airmass is added to the ::smparmodel
     * structure. It is computed from the altitude angle by using the formula
     * of Rozenberg (1966).
     *
     * \b INPUT:
     * \param params    CPL sky model parameter list
     *
     * \b OUTPUT:
     * \param modelpar  ::smparmodel structure containing the read parameters
     *
     * \b ERRORS:
     * - see subroutines
     */

    const cpl_parameter *p;
    cpl_array *emis = NULL, *temp = NULL;
    char emis_str[SM_MAXLEN+1], temp_str[SM_MAXLEN+1];
    int nemis = 0, ntemp = 0, i = 0;

    /* Read CPL parameter list and fill smparmodel structure */

    p = sm_etc_parameterlist_find_const(params, "sm_h", CPL_TYPE_DOUBLE);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    modelpar->sm_h = cpl_parameter_get_double(p);

    p = sm_etc_parameterlist_find_const(params, "sm_hmin", CPL_TYPE_DOUBLE);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    modelpar->sm_hmin = cpl_parameter_get_double(p);

    p = sm_etc_parameterlist_find_const(params, "alt", CPL_TYPE_DOUBLE);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    modelpar->alt = cpl_parameter_get_double(p);

    p = sm_etc_parameterlist_find_const(params, "alpha", CPL_TYPE_DOUBLE);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    modelpar->alpha = cpl_parameter_get_double(p);

    p = sm_etc_parameterlist_find_const(params, "rho", CPL_TYPE_DOUBLE);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    modelpar->rho = cpl_parameter_get_double(p);

    p = sm_etc_parameterlist_find_const(params, "altmoon", CPL_TYPE_DOUBLE);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    modelpar->altmoon = cpl_parameter_get_double(p);

    p = sm_etc_parameterlist_find_const(params, "moondist", CPL_TYPE_DOUBLE);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    modelpar->moondist = cpl_parameter_get_double(p);

    p = sm_etc_parameterlist_find_const(params, "pres", CPL_TYPE_DOUBLE);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    modelpar->pres = cpl_parameter_get_double(p);

    p = sm_etc_parameterlist_find_const(params, "ssa", CPL_TYPE_DOUBLE);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    modelpar->ssa = cpl_parameter_get_double(p);

    p = sm_etc_parameterlist_find_const(params, "calcds", CPL_TYPE_STRING);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    strncpy(modelpar->calcds, cpl_parameter_get_string(p), SM_LENLINE+1);

    p = sm_etc_parameterlist_find_const(params, "o3column", CPL_TYPE_DOUBLE);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    modelpar->o3column = cpl_parameter_get_double(p);

    p = sm_etc_parameterlist_find_const(params, "moonscal", CPL_TYPE_DOUBLE);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    modelpar->moonscal = cpl_parameter_get_double(p);

    p = sm_etc_parameterlist_find_const(params, "lon_ecl", CPL_TYPE_DOUBLE);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    modelpar->lon_ecl = cpl_parameter_get_double(p);

    p = sm_etc_parameterlist_find_const(params, "lat_ecl", CPL_TYPE_DOUBLE);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    modelpar->lat_ecl = cpl_parameter_get_double(p);

    p = sm_etc_parameterlist_find_const(params, "emis_str", CPL_TYPE_STRING);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    strcpy(emis_str, cpl_parameter_get_string(p));

    p = sm_etc_parameterlist_find_const(params, "temp_str", CPL_TYPE_STRING);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    strcpy(temp_str, cpl_parameter_get_string(p));

    p = sm_etc_parameterlist_find_const(params, "msolflux", CPL_TYPE_DOUBLE);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    modelpar->msolflux = cpl_parameter_get_double(p);

    p = sm_etc_parameterlist_find_const(params, "season", CPL_TYPE_INT);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    modelpar->season = cpl_parameter_get_int(p);

    p = sm_etc_parameterlist_find_const(params, "time", CPL_TYPE_INT);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    modelpar->time = cpl_parameter_get_int(p);

    p = sm_etc_parameterlist_find_const(params, "vac_air", CPL_TYPE_STRING);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    strncpy(modelpar->vac_air, cpl_parameter_get_string(p), SM_LENLINE+1);

    p = sm_etc_parameterlist_find_const(params, "pwv", CPL_TYPE_DOUBLE);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    modelpar->pwv = cpl_parameter_get_double(p);

    p = sm_etc_parameterlist_find_const(params, "rtcode", CPL_TYPE_STRING);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    strncpy(modelpar->rtcode, cpl_parameter_get_string(p), SM_LENLINE+1);

    p = sm_etc_parameterlist_find_const(params, "resol", CPL_TYPE_DOUBLE);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    modelpar->resol = cpl_parameter_get_double(p);

    p = sm_etc_parameterlist_find_const(params, "incl", CPL_TYPE_STRING);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    strncpy(modelpar->incl, cpl_parameter_get_string(p), SM_LENLINE+1);

    /* Handle emissivity and temperature strings */

    emis = cpl_array_new(0, CPL_TYPE_DOUBLE);
    sm_etc_splitstring(emis, emis_str);
    nemis = cpl_array_get_size(emis);
    temp = cpl_array_new(0, CPL_TYPE_DOUBLE);
    sm_etc_splitstring(temp, temp_str);
    ntemp = cpl_array_get_size(temp);

    if (nemis != 0 && ntemp == nemis) {
        modelpar->ncomp = nemis;
        for (i = 0; i < modelpar->ncomp; i++) {
            modelpar->eps[i] = cpl_array_get(emis, i, NULL);
            modelpar->temp[i] = cpl_array_get(temp, i, NULL);
        }
    } else {
        modelpar->ncomp = 0;
    }

    cpl_array_delete(emis);
    cpl_array_delete(temp);

    /* Compute airmass (Rozenberg 1966) and write result to smparmodel
       structure */
    modelpar->airmass = sm_comp_alttoairmass(modelpar->alt);

    return cpl_error_get_code();
}


const cpl_parameter
*sm_etc_parameterlist_find_const(const cpl_parameterlist *self,
                                 const char *name, cpl_type type)
{
    /*!
     * Error-handling wrapper around cpl_parameterlist_find_const(),
     * catching errors of mismatching parameter name and/or type.
     *
     * \b INPUT:
     * \param self  CPL parameter list
     * \param name  parameter name
     * \param type  CPL type of parameter
     *
     * \b RETURN:
     * - CPL parameter
     *
     * \b ERRORS:
     * - CPL_ERROR_DATA_NOT_FOUND
     * - CPL_ERROR_TYPE_MISMATCH
     */

    const cpl_parameter *p = cpl_parameterlist_find_const(self, name);

    if (p == NULL) {
        (void) cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                     "parameter not found: %s", name);
    } else if (cpl_parameter_get_type(p) != type) {
        (void) cpl_error_set_message(cpl_func, CPL_ERROR_TYPE_MISMATCH,
                                 "parameter %s has type %s, expected type %s",
                                 name,
                                 cpl_type_get_name(cpl_parameter_get_type(p)),
                                 cpl_type_get_name(type));
        p = NULL;
    }

    return p;
}


cpl_error_code sm_etc_splitstring(cpl_array *outval, char *instring)
{
    /*!
     * Splits a string of comma-separated values and converts the substrings
     * into double values.
     *
     * \b INPUT:
     * \param instring  comma-separated list of values
     *
     * \b OUTPUT:
     * \param outval    CPL array of double values
     *
     * \b ERRORS:
     * - ISM: Insufficient memory
     * - NDA: No data
     */

    char **str = NULL, errtxt[SM_MAXLEN+1];
    int i = 0, nval = 0;

    /* Allocate memory to string array */

    str = (char **) calloc(SM_MAXPAR, sizeof(char *));
    if (str == NULL) {
        sprintf(errtxt, "%s: char **str", SM_ERROR_ISM_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_ISM, "%s", errtxt);
    }

    /* Read parameter values in strings */

    str[0] = strtok(instring, ", \n");
    if (str[0] == NULL) {
        free(str);
        str = NULL;
        sprintf(errtxt, "%s: char *instring", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    }

    for (i = 1; i < SM_MAXPAR; i++) {
        str[i] = strtok(NULL, ", \n");
        if (str[i] == NULL) {
            nval = i;
            break;
        }
    }

    /* Create CPL double array */

    cpl_array_set_size(outval, nval);
    if (outval == NULL) {
        free(str);
        str = NULL;
        sprintf(errtxt, "%s: cpl_array *outval", SM_ERROR_ISM_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_ISM, "%s", errtxt);
    }

    /* Convert strings in double numbers */

    for (i = 0; i < nval; i++) {
        cpl_array_set_double(outval, i, strtod(str[i], NULL));
    }

    /* Free memory space of string array */

    free(str);
    str = NULL;

    return CPL_ERROR_NONE;
}


cpl_error_code sm_etc_readfilenames(smparmodel *modelpar,
                                    const cpl_parameterlist *params)
{
    /*!
     * Reads data paths and file names from a file (path from CPL parameter
     * list and name from header file) and put them in a ::smparmodel
     * structure.
     *
     * \b INPUT:
     * \param params    CPL sky model parameter list
     *
     * \b OUTPUT:
     * \param modelpar  ::smparmodel structure containing the read parameters
     *
     * \b ERRORS:
     * - FOF: File opening failed
     * - see ::sm_param_readcheck as well
     */

    cpl_error_code status = CPL_ERROR_NONE;
    FILE *stream;
    smparam x[SM_MAXPAR];
    const cpl_parameter *p;
    char filepath[FILENAME_MAX], parfile[FILENAME_MAX], errtxt[SM_MAXLEN+1];
    char libpath[FILENAME_MAX], file[FILENAME_MAX], datapath[FILENAME_MAX];

    /* Check existence of parameter file */

    p = sm_etc_parameterlist_find_const(params, "filepath", CPL_TYPE_STRING);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    strncpy(filepath, cpl_parameter_get_string(p), FILENAME_MAX);
    sprintf(parfile, "%s/%s", filepath, SM_FILENAMELIST);
    if ((stream = fopen(parfile, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", SM_ERROR_FOF_TXT, parfile);
        return cpl_error_set_message(cpl_func, SM_ERROR_FOF, "%s", errtxt);
    }

    /* Read parameter file */

    sm_param_readcheck(stream, x, "libpath", 1);
    if (x[1].c[0] != '/') {
        sprintf(libpath, "%s/%s", filepath, x[1].c);
    } else {
        sprintf(libpath, "%s", x[1].c);
    }
    if ((status = sm_basic_access(libpath, F_OK)) != CPL_ERROR_NONE) {
        fclose(stream);
        return status;
    }
    strcpy(modelpar->libpath, libpath);

    sm_param_readcheck(stream, x, "libstruct", 2);
    strcpy(modelpar->libstruct1, x[1].c);
    strcpy(modelpar->libstruct2, x[2].c);
    p = sm_etc_parameterlist_find_const(params, "pwv", CPL_TYPE_DOUBLE);
    if (p == NULL)  {
        fclose(stream);
        return cpl_error_set_where(cpl_func);
    }
    if (cpl_parameter_get_double(p) >= 0.) {
        /* For positive PWV take PWV-dependent library */
        strcpy(modelpar->libstruct, x[2].c);
        sprintf(file, "%s/%s", libpath, x[2].c);
    } else {
        /* For negative PWV take time-dependent library */
        strcpy(modelpar->libstruct, x[1].c);
        sprintf(file, "%s/%s", libpath, x[1].c);
    }
    sm_basic_access(file, F_OK);

    sm_param_readcheck(stream, x, "datapath", 1);
    if (x[1].c[0] != '/') {
        sprintf(datapath, "%s/%s", filepath, x[1].c);
    } else {
        sprintf(datapath, "%s", x[1].c);
    }
    if ((status = sm_basic_access(datapath, F_OK)) != CPL_ERROR_NONE) {
       fclose(stream);
       return status;
    }
    strcpy(modelpar->datapath, datapath);

    sm_param_readcheck(stream, x, "solspecname", 1);
    strcpy(modelpar->solspecname, x[1].c);
    sprintf(file, "%s/%s", datapath, x[1].c);
    sm_basic_access(file, F_OK);

    sm_param_readcheck(stream, x, "mieextname", 1);
    strcpy(modelpar->mieextname, x[1].c);
    if (strncmp(x[1].c, "NONE", 4) != 0) {
        sprintf(file, "%s/%s", datapath, x[1].c);
        sm_basic_access(file, F_OK);
    }

    sm_param_readcheck(stream, x, "lunirrname", 1);
    strcpy(modelpar->lunirrname, x[1].c);
    sprintf(file, "%s/%s", datapath, x[1].c);
    sm_basic_access(file, F_OK);

    sm_param_readcheck(stream, x, "miephasename", 1);
    strcpy(modelpar->miephasename, x[1].c);
    sprintf(file, "%s/%s", datapath, x[1].c);
    sm_basic_access(file, F_OK);

    sm_param_readcheck(stream, x, "sscatcorname", 1);
    strcpy(modelpar->sscatcorname, x[1].c);
    sprintf(file, "%s/%s", datapath, x[1].c);
    sm_basic_access(file, F_OK);

    sm_param_readcheck(stream, x, "dscatcorname", 1);
    strcpy(modelpar->dscatcorname, x[1].c);
    sprintf(file, "%s/%s", datapath, x[1].c);
    sm_basic_access(file, F_OK);

    sm_param_readcheck(stream, x, "o3transname", 1);
    strcpy(modelpar->o3transname, x[1].c);
    sprintf(file, "%s/%s", datapath, x[1].c);
    sm_basic_access(file, F_OK);

    sm_param_readcheck(stream, x, "starspecname", 1);
    strcpy(modelpar->starspecname, x[1].c);
    sprintf(file, "%s/%s", datapath, x[1].c);
    sm_basic_access(file, F_OK);

    sm_param_readcheck(stream, x, "zodtabname", 1);
    strcpy(modelpar->zodtabname, x[1].c);
    sprintf(file, "%s/%s", datapath, x[1].c);
    sm_basic_access(file, F_OK);

    sm_param_readcheck(stream, x, "linetabname", 1);
    strcpy(modelpar->linetabname, x[1].c);
    sprintf(file, "%s/%s", datapath, x[1].c);
    sm_basic_access(file, F_OK);

    sm_param_readcheck(stream, x, "vardatname", 1);
    strcpy(modelpar->vardatname, x[1].c);
    sprintf(file, "%s/%s", datapath, x[1].c);
    sm_basic_access(file, F_OK);

    sm_param_readcheck(stream, x, "acontname", 1);
    strcpy(modelpar->acontname, x[1].c);
    sprintf(file, "%s/%s", datapath, x[1].c);
    sm_basic_access(file, F_OK);

    fclose(stream);

    return cpl_error_get_code();
}


cpl_error_code sm_etc_writetable(cpl_table *skytable, const smspec *radiance,
                                 const smspec *transmission,
                                 const cpl_table *radcomp,
                                 const cpl_table *transcomp)
{
    /*!
     * Writes radiance spectrum and transmission curve which are provided as
     * ::smspec structures to the CPL table of the ETC interface. In addition,
     * the individual sky model components are also written into the output
     * table.
     *
     * \b INPUT:
     * \param skytable      input wavelength grid as CPL table
     * \param radiance      radiance spectrum as ::smspec structure
     * \param transmission  transmission curve as ::smspec structure
     * \param radcomp       sky radiance components of sky model as CPL table
     * \param transcomp     transmission components of sky model as CPL table
     *
     * \b OUTPUT:
     * \param skytable      CPL table with full radiance and transmission data
     *
     * \b ERRORS:
     * - NDA: No wavelength data
     * - IDG: Inconsistent data grids
     */

    cpl_array *colnames;
    char errtxt[SM_MAXLEN+1], colname[SM_MAXLEN+1];
    int n = 0, nrow = 0, i = 0, ncol = 0, j = 0;
    double lam = 0.;

    /* Check wavelength grids of input spectra */

    if (radiance->n == 0 || transmission->n == 0) {
        sprintf(errtxt, "%s: smspec *radiance and/or smspec *transmission",
                SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    } else if ((int) sm_spec_compgrids(radiance, transmission) ==
            (int) SM_ERROR_IDG) {
        sprintf(errtxt, "%s: smspec *radiance != smspec *transmission",
                SM_ERROR_IDG_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IDG, "%s", errtxt);
    }

    /* Check wavelength grid of output CPL table */

    n = cpl_table_get_nrow(skytable);
    if (n == 0) {
        sprintf(errtxt, "%s: cpl_table *skytable", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    }

    for (i = 0; i < n; i++) {
        lam = cpl_table_get_double(skytable, "lam", i, NULL);
        if (lam != radiance->dat[i].lam) {
            sprintf(errtxt, "%s: cpl_table *skytable != smspec *radiance",
                    SM_ERROR_IDG_TXT);
            return cpl_error_set_message(cpl_func, SM_ERROR_IDG, "%s",
                                         errtxt);
        }
    }

    /* Check number of wavelengths in sky radiance component table */

    nrow = cpl_table_get_nrow(radcomp);
    if (nrow != n) {
        sprintf(errtxt, "%s: cpl_table *radcomp != cpl_table *skytable "
                "(# of rows)", SM_ERROR_IDG_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IDG, "%s", errtxt);
    }

    /* Check number of wavelengths in transmission component table */

    nrow = cpl_table_get_nrow(transcomp);
    if (nrow != n) {
        sprintf(errtxt, "%s: cpl_table *transcomp != cpl_table *skytable "
                "(# of rows)", SM_ERROR_IDG_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IDG, "%s", errtxt);
    }

    /* Create missing columns of CPL table */

    cpl_table_new_column(skytable, "flux", CPL_TYPE_DOUBLE);
    cpl_table_new_column(skytable, "dflux1", CPL_TYPE_DOUBLE);
    cpl_table_new_column(skytable, "dflux2", CPL_TYPE_DOUBLE);
    cpl_table_new_column(skytable, "trans", CPL_TYPE_DOUBLE);
    cpl_table_new_column(skytable, "dtrans1", CPL_TYPE_DOUBLE);
    cpl_table_new_column(skytable, "dtrans2", CPL_TYPE_DOUBLE);

    /* Fill new columns of CPL table with data from the input spectra */

    for (i = 0; i < n; i++) {
        cpl_table_set_double(skytable, "flux", i, radiance->dat[i].flux);
        cpl_table_set_double(skytable, "dflux1", i, radiance->dat[i].dflux1);
        cpl_table_set_double(skytable, "dflux2", i, radiance->dat[i].dflux2);
        cpl_table_set_double(skytable, "trans", i, transmission->dat[i].flux);
        cpl_table_set_double(skytable, "dtrans1", i,
                             transmission->dat[i].dflux1);
        cpl_table_set_double(skytable, "dtrans2", i,
                             transmission->dat[i].dflux2);
    }

    /* Copy columns of sky radiance component table to output table */

    ncol = cpl_table_get_ncol(radcomp);
    colnames = cpl_table_get_column_names(radcomp);
    for (j = 0; j < ncol; j++) {
        sprintf(colname, "%s", cpl_array_get_string(colnames, j));
        cpl_table_duplicate_column(skytable, colname, radcomp, colname);
    }
    cpl_array_delete(colnames);

    /* Copy columns of transmission component table to output table */

    ncol = cpl_table_get_ncol(transcomp);
    colnames = cpl_table_get_column_names(transcomp);
    for (j = 0; j < ncol; j++) {
        sprintf(colname, "%s", cpl_array_get_string(colnames, j));
        cpl_table_duplicate_column(skytable, colname, transcomp, colname);
    }
    cpl_array_delete(colnames);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_comp_skyemcomp(smspec *radiance, smspec *transmission,
                                 cpl_table *radcomp, cpl_table *transcomp,
                                 const smparmodel modelpar)
{
    /*!
     * \callgraph
     * Creates sky emission spectrum and transmission curve (with
     * uncertainties) at object position for given wavelength grid and input
     * parameters.
     * Imports solar spectrum, LBLRTM or RFM radiance and transmission spectra
     * (ASCII or FITS files), and airglow emission line and continuum table,
     * estimates lunar and zodical sky brightness, considers scattered
     * starlight, and calculates thermal telescope/instrument emission.
     * Using the "incl" parameter of the ::smparmodel structure, the selection
     * of the seven different components for the computation of the output
     * radiance spectrum can be modified.
     * In addition to the combined specta, the routine also provides the
     * individual emission and transmission components (without uncertainties)
     * by two CPL tables. The number of columns of these tables corresponds to
     * the number of components.
     *
     * \b INPUT:
     * \param radiance      wavelength grid of sky emission spectrum
     * \param transmission  wavelength grid of transmission curve
     * \param radcomp       empty CPL table
     * \param transcomp     empty CPL table
     * \param modelpar      input parameters related to sky emission
     *                      (see typedef of ::smparmodel)
     *
     * \b OUTPUT:
     * \param radiance      sky emission spectrum
     *                      (::smspec structure, errors considered),
     *                      units: \f$\mu{\rm m}\f$,
     *    \f${\rm phot\,s}^{-1}{\rm m}^{-2}\mu{\rm m}^{-1}{\rm arcsec}^{-2}\f$
     * \param transmission  transmission curve
     *                      (::smspec structure, errors considered),
     *                      units: \f$\mu{\rm m}\f$, [0,1]
     * \param radcomp       sky radiance components of sky model as CPL table
     * \param transcomp     transmission components of sky model as CPL table
     *
     * \b ERRORS:
     * - NDA: No data
     * - IDG: Inconsistent data grids
     * - see subroutines as well
     */

    cpl_error_code status = CPL_ERROR_NONE;
    smspec lamgrid, solspec, intsolspec, libtrans, molabstrans, abstrans;
    smspec o3trans, rscattrans, mscattrans, lunskybright, scatstarlight;
    smspec zodskybright, telem, librad, skylines, airglowcont;
    char errtxt[SM_MAXLEN+1], solspecfile[FILENAME_MAX];
    double lim[2] = {0, HUGE_VAL};

    /* Check input wavelength grids */
    if (radiance->n == 0 || transmission->n == 0) {
        sprintf(errtxt, "%s: smspec *radiance and/or smspec *transmission",
                SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    } else if ((int) sm_spec_compgrids(radiance, transmission) ==
            (int) SM_ERROR_IDG) {
        sprintf(errtxt, "%s: smspec *radiance != smspec *transmission",
                SM_ERROR_IDG_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IDG, "%s", errtxt);
    }

    /* Initialise output wavelength grid */
    sm_spec_copy(&lamgrid, radiance);
    lamgrid.type = 2;
    radiance->type = 4;
    transmission->type = 4;

    /* Read and interpolate solar spectrum */
    sprintf(solspecfile, "%s/%s", modelpar.datapath, modelpar.solspecname);
    sm_spec_read(&solspec, solspecfile);
    sm_spec_copy(&intsolspec, &lamgrid);
    sm_spec_interpol(&intsolspec, &solspec);
    sm_spec_free(&solspec);

    /* Read, rebin, and scale (to airmass = 1) LBLRTM/RFM transmission
       spectrum */
    sm_spec_copy(&libtrans, &lamgrid);
    sm_comp_extrapoltrans(&libtrans, modelpar, lim);

    /* Exit if errors occurred */
    if ((status = cpl_error_get_code()) != CPL_ERROR_NONE) {
        sm_spec_free(&lamgrid);
        sm_spec_free(&intsolspec);
        sm_spec_free(&libtrans);
        return status;
    }

    /* Separate UV/optical ozone absorption from molecular absorption and
       scale it to given relative ozone column density */
    sm_comp_scaleo3(&molabstrans, &abstrans, &o3trans, &libtrans, modelpar);

    /* Calculate Rayleigh and Mie (aerosol) extinction curves */
    sm_spec_copy(&rscattrans, &lamgrid);
    sm_comp_calcrayleighscat(&rscattrans, modelpar);
    sm_spec_copy(&mscattrans, &lamgrid);
    sm_comp_calcmiescat(&mscattrans, modelpar);

    /* Merge curves for molecular absorption, Rayleigh scattering, and aerosol
       extinction */
    sm_spec_calc(transmission, '=', &molabstrans);
    sm_spec_calc(transmission, '*', &rscattrans);
    sm_spec_calc(transmission, '*', &mscattrans);

    /* Scale transmission curve depending on object airmass */
    sm_comp_scaletranscurv(transmission, modelpar);

    /* Write transmission curve components for object airmass into CPL
       table */
    sm_comp_writetranscomp(transcomp, "trans_ma", &abstrans, modelpar);
    sm_comp_writetranscomp(transcomp, "trans_o3", &o3trans, modelpar);
    sm_comp_writetranscomp(transcomp, "trans_rs", &rscattrans, modelpar);
    sm_comp_writetranscomp(transcomp, "trans_ms", &mscattrans, modelpar);

    /* Get error code */
    status = cpl_error_get_code();

    /* Add sky emission components */

    if (status == CPL_ERROR_NONE &&
        (modelpar.incl[0] == 'Y' || modelpar.incl[0] == 'y')) {
        /* Compute lunar sky brightness */
        sm_spec_copy(&lunskybright, &lamgrid);
        status = sm_comp_lunskybright(&lunskybright, modelpar, &intsolspec,
                                      &abstrans, &o3trans, &rscattrans,
                                      &mscattrans);
        sm_spec_calc(radiance, '+', &lunskybright);
        sm_spec_writecplcolumn(radcomp, "flux_sml", &lunskybright, 2);
        sm_spec_free(&lunskybright);
    } else {
        sm_spec_writecplcolumn(radcomp, "flux_sml", &lamgrid, 2);
    }

    if (status == CPL_ERROR_NONE &&
        (modelpar.incl[1] == 'Y' || modelpar.incl[1] == 'y')) {
        /* Get spectrum of scattered starlight */
        sm_spec_copy(&scatstarlight, &lamgrid);
        status = sm_comp_scatstarlight(&scatstarlight, modelpar,
                                       &molabstrans);
        sm_spec_calc(radiance, '+', &scatstarlight);
        sm_spec_writecplcolumn(radcomp, "flux_ssl", &scatstarlight, 2);
        sm_spec_free(&scatstarlight);
    } else {
        sm_spec_writecplcolumn(radcomp, "flux_ssl", &lamgrid, 2);
    }

    if (status == CPL_ERROR_NONE &&
        (modelpar.incl[2] == 'Y' || modelpar.incl[2] == 'y')) {
        /* Compute zodiacal sky brightness */
        sm_spec_copy(&zodskybright, &lamgrid);
        status = sm_comp_zodskybright(&zodskybright, modelpar, &intsolspec,
                                      &molabstrans, &rscattrans, &mscattrans);
        sm_spec_calc(radiance, '+', &zodskybright);
        sm_spec_writecplcolumn(radcomp, "flux_zl", &zodskybright, 2);
        sm_spec_free(&zodskybright);
    } else {
        sm_spec_writecplcolumn(radcomp, "flux_zl", &lamgrid, 2);
    }

    if (status == CPL_ERROR_NONE &&
        (modelpar.incl[3] == 'Y' || modelpar.incl[3] == 'y')) {
        /* Compute telescope/instrument emission */
        sm_spec_copy(&telem, &lamgrid);
        status = sm_comp_telem(&telem, modelpar);
        sm_spec_calc(radiance, '+', &telem);
        sm_spec_writecplcolumn(radcomp, "flux_tie", &telem, 2);
        sm_spec_free(&telem);
    } else {
        sm_spec_writecplcolumn(radcomp, "flux_tie", &lamgrid, 2);
    }

    if (status == CPL_ERROR_NONE &&
        (modelpar.incl[4] == 'Y' || modelpar.incl[4] == 'y')) {
        /* Read and rebin LBLRTM/RFM radiance spectrum.
           Neglect short wavelengths with very low fluxes. */
        sm_spec_copy(&librad, &lamgrid);
        librad.type = 4; /* consider asymmetric errors */
        lim[0] = SM_RADMINLAM;
        status = sm_comp_extrapolrad(&librad, &molabstrans, modelpar, lim);
        sm_spec_calc(radiance, '+', &librad);
        sm_spec_writecplcolumn(radcomp, "flux_tme", &librad, 2);
        sm_spec_free(&librad);
    } else {
        sm_spec_writecplcolumn(radcomp, "flux_tme", &lamgrid, 2);
    }

    if (status == CPL_ERROR_NONE &&
        (modelpar.incl[5] == 'Y' || modelpar.incl[5] == 'y')) {
        /* Read airglow emission lines and compute spectrum */
        sm_spec_copy(&skylines, &lamgrid);
        skylines.type = 3; /* consider uncertainties */
        status = sm_comp_getlinespec(&skylines, modelpar, &rscattrans,
                                     &mscattrans);
        sm_spec_calc(radiance, '+', &skylines);
        sm_spec_writecplcolumn(radcomp, "flux_ael", &skylines, 2);
        sm_spec_free(&skylines);
    } else {
        sm_spec_writecplcolumn(radcomp, "flux_ael", &lamgrid, 2);
    }

    if (status == CPL_ERROR_NONE &&
        (modelpar.incl[6] == 'Y' || modelpar.incl[6] == 'y')) {
        /* Compute airglow continuum */
        sm_spec_copy(&airglowcont, &lamgrid);
        airglowcont.type = 3; /* consider uncertainties */
        status = sm_comp_airglowcont(&airglowcont, modelpar, &molabstrans,
                                     &rscattrans, &mscattrans);
        sm_spec_calc(radiance, '+', &airglowcont);
        sm_spec_writecplcolumn(radcomp, "flux_arc", &airglowcont, 2);
        sm_spec_free(&airglowcont);
    } else {
        sm_spec_writecplcolumn(radcomp, "flux_arc", &lamgrid, 2);
    }

    /* Free remaining routine pointers */
    sm_spec_free(&lamgrid);
    sm_spec_free(&intsolspec);
    sm_spec_free(&libtrans);
    sm_spec_free(&molabstrans);
    sm_spec_free(&abstrans);
    sm_spec_free(&o3trans);
    sm_spec_free(&rscattrans);
    sm_spec_free(&mscattrans);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_comp_extrapoltrans(smspec *spec, const smparmodel modelpar,
                                     const double lim[2])
{
    /*!
     * \callgraph
     * Extracts and rebins library LBLRTM/RFM molecular transmission spectrum
     * and extrapolates it to airmass = 1. For the airmass calculations, the
     * formula of Rozenberg (1966) is used.
     *
     * \b INPUT:
     * \param spec       desired wavelength grid
     * \param modelpar   sky emission parameters (see typedef of ::smparmodel)
     * \param lim        minimum and maximum wavelength extraction limits
     *
     * \b OUTPUT:
     * \param spec       extracted, rebinned, and scaled LBLRTM/RFM
     *                   transmission curve
     *
     * \b ERRORS:
     * - NDA: No data
     * - IOD: Invalid order of data points
     * - IIP: Invalid input parameter(s)
     */

    cpl_error_code err = CPL_ERROR_NONE;
    cpl_parameter *p;
    cpl_parameterlist *libfilepar;
    char errtxt[SM_MAXLEN+1];
    double xz = 1.;

    /* spec = NULL */
    if (spec->n <= 0) {
        sprintf(errtxt, "%s: smspec *spec", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    }

    /* Check input wavelength limits */
    if (lim[1] <= lim[0]) {
        sprintf(errtxt, "%s: lim[1] <= lim[0] (input wavelength limits)",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Valid target altitude? */
    if (modelpar.alt < 0 || modelpar.alt > 90) {
        sprintf(errtxt, "%s: alt of smparmodel modelpar < 0 or > 90 degrees",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Search in the library of molecular spectra for a file with the desired
       parameters */
    libfilepar = cpl_parameterlist_new();
    if ((err = sm_comp_readlibstruct(libfilepar, modelpar, "T")) !=
        CPL_ERROR_NONE) {
        cpl_parameterlist_delete(libfilepar);
        return err;
    }

    /* Get transmission curve from library */
    if ((err = sm_comp_getmolspec(spec, libfilepar, lim)) !=
        CPL_ERROR_NONE) {
        cpl_parameterlist_delete(libfilepar);
        return err;
    }

    /* Get airmass of selected library file */
    if ((p = cpl_parameterlist_find(libfilepar, "airmass")) != NULL) {
        xz = cpl_parameter_get_double(p);
    } else if ((p = cpl_parameterlist_find(libfilepar, "alt")) != NULL) {
        /* Airmass according to Rozenberg (1966) */
        xz = sm_comp_alttoairmass(cpl_parameter_get_double(p));
    } else {
        /* Assume airmass = 1 */
        xz = 1;
    }

    /* Scale transmission curve to airmass = 1 */
    sm_spec_scale(spec, '^', 1. / xz);

    /* Free memory */
    cpl_parameterlist_delete(libfilepar);

    return err;
}


double sm_comp_alttoairmass(const double alt)
{
    /*!
     * Calculates airmass from altitude angle in deg using the formula of
     * Rozenberg (1966). The maximum airmass is 40.
     *
     * \b INPUT:
     * \param alt  altitude angle in deg

     * \b RETURN:
     * - airmass
     *
     * \b ERRORS:
     * - none
     */

    double altc = 0., z = 0., xz = 0.;

    /* Limit alitude angles */
    if (alt < 0.) {
        altc = 0.;
    } else if (alt > 90.) {
        altc = 90.;
    } else {
        altc = alt;
    }

    /* Convert altitude in deg into zenith distance in rad */
    z = (90. - altc) * CPL_MATH_RAD_DEG;

    /* Use airmass approximation of Rozenberg (1966) */
    xz = 1 / (cos(z) + 0.025 * exp(-11 * cos(z)));

    /* Return airmass */
    return xz;
}


cpl_error_code sm_comp_getmolspec(smspec *spec,
                                  const cpl_parameterlist *libfilepar,
                                  const double lim[2])
{
    /*!
     * \callgraph
     * Extracts and rebins library LBLRTM/RFM molecular radiance or
     * transmission spectra.
     * The file name and the decision between vacuum and air wavelengths have
     * to be provided by a CPL parameter list.
     * Each library file has to consist of the four columns wavelength,
     * radiance/transmission for average atmospheric profiles, and
     * radiance/transmission for minus and plus sigma atmospheric profiles.
     * This structure is converted to the ::smspec structure, i.e. "lam",
     * "flux", "dflux1" (negative error), and "dflux2" (positive error).
     * The extraction-related wavelength range can be delimitated regardless
     * of the desired grid in order to reduce computing time by avoiding the
     * extraction of wavelengths that do not contribute significantly to the
     * composite sky spectrum.
     *
     * \b INPUT:
     * \param spec        desired wavelength grid
     * \param libfilepar  CPL parameter list containing file name (+ path) and
     *                    vacuum/air parameter
     * \param lim         minimum and maximum wavelength extraction limits
     *
     * \b OUTPUT:
     * \param spec        extracted and rebinned LBLRTM/RFM spectrum (radiance
     *                    or transmission)
     *
     * \b ERRORS:
     * - NDA: No data
     * - IOD: Invalid order of data points
     * - IIP: Invalid input parameter(s)
     * - IOS: Invalid object structure
     * - FOF: File opening failed
     * - IFE: Invalid file name extension
     * - CPL: Access out of range [warning]
     */

    FILE *stream;
    cpl_error_code err = CPL_ERROR_NONE;
    const cpl_parameter *p;
    smspec libspec;
    smbool isoutrange = F;
    char errtxt[SM_MAXLEN+1], vac_air[4], libfile[FILENAME_MAX];
    int i = 0;
    double speclim[2] = {0, HUGE_VAL}, mlim[2] = {0, HUGE_VAL};
    double limlam[2] = {0, HUGE_VAL}, mflux = 0., pflux = 0.;

    /* spec = NULL */
    if (spec->n <= 0) {
        sprintf(errtxt, "%s: smspec *spec", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    }

    /* Check wavelength limits */
    speclim[0] = spec->dat[0].lam;
    speclim[1] = spec->dat[(spec->n)-1].lam;
    if (speclim[1] <= speclim[0]) {
        sprintf(errtxt, "%s: smspec *spec (wavelength grid)",
                SM_ERROR_IOD_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IOD, "%s", errtxt);
    }

    /* Check input wavelength limits */
    if (lim[1] <= lim[0]) {
        sprintf(errtxt, "%s: lim[1] <= lim[0] (input wavelength limits)",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Vacuum or air wavelengths? */
    p = sm_etc_parameterlist_find_const(libfilepar, "vac_air",
                                        CPL_TYPE_STRING);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    strncpy(vac_air, cpl_parameter_get_string(p), 4);

    /* Derive limits for extraction */
    if (speclim[0] < lim[0]) {
        mlim[0] = lim[0];
    } else {
        mlim[0] = speclim[0];
    }
    if (speclim[1] > lim[1]) {
        mlim[1] = lim[1];
    } else {
        mlim[1] = speclim[1];
    }

    /* Modify limits if air wavelengths are desired */
    if (strncmp(vac_air, "air", 3) == 0) {
        limlam[0] = sm_comp_airtovac_single(mlim[0]);
        limlam[1] = sm_comp_airtovac_single(mlim[1]);
    } else {
        limlam[0] = mlim[0];
        limlam[1] = mlim[1];
    }

    if (limlam[1] <= limlam[0]) {
        /* Desired wavelength range beyond extraction limits
           -> return 0 for all wavelengths */
        return CPL_ERROR_NONE;
    }

    /* Get file name from parameter list */
    p = sm_etc_parameterlist_find_const(libfilepar, "filename",
                                        CPL_TYPE_STRING);
    if (p == NULL) return cpl_error_set_where(cpl_func);
    strncpy(libfile, cpl_parameter_get_string(p), FILENAME_MAX);

    /* Check file existence */
    if ((stream = fopen(libfile, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", SM_ERROR_FOF_TXT, libfile);
        return cpl_error_set_message(cpl_func, SM_ERROR_FOF, "%s", errtxt);
    }
    fclose(stream);

    /* Check file name extension and read LBLRTM/RFM spectrum
       (each file has to contain a wavelength column and three flux
       columns) */

    if (strstr(libfile, ".dat") != NULL ||
        strstr(libfile, ".ascii") != NULL) {
        sm_spec_readrange(&libspec, libfile, limlam, SM_RRSTEP);
    } else if (strstr(libfile, ".fits") != NULL ||
               strstr(libfile, ".mt") != NULL) {
        sm_spec_readfitsrange(&libspec, libfile, limlam, SM_RRSTEP);
    } else {
        libspec.n = 0;
        sprintf(errtxt, "%s: %s ('dat', 'ascii', 'fits', or 'mt' only)",
                SM_ERROR_IFE_TXT, libfile);
        err = cpl_error_set_message(cpl_func, SM_ERROR_IFE, "%s", errtxt);
    }

    /* Convert vacuum to air wavelengths if required */
    if (strncmp(vac_air, "air", 3) == 0) {
        sm_comp_vactoair_spec(&libspec);
    }

    /* Rebin LBLRTM/RFM spectrum */
    sm_spec_rebin(spec, &libspec);
    sm_spec_free(&libspec);

    /* Change library spectrum structure to smspec structure */

    for (i = 0; i < spec->n; i++) {
        if (spec->dat[i].dflux1 > spec->dat[i].dflux2) {
            mflux = spec->dat[i].dflux2;
            pflux = spec->dat[i].dflux1;
        } else {
            mflux = spec->dat[i].dflux1;
            pflux = spec->dat[i].dflux2;
        }
        spec->dat[i].dflux1 = spec->dat[i].flux - mflux;
        spec->dat[i].dflux2 = pflux - spec->dat[i].flux;
        if (spec->dat[i].dflux1 < 0) {
            spec->dat[i].dflux1 = 0;
            /*isoutrange = T;*/
        }
        if (spec->dat[i].dflux2 < 0) {
            spec->dat[i].dflux2 = 0;
            /*isoutrange = T;*/
        }
    }

    /* Set warning in the case of negative errors */
    if (isoutrange == T) {
        sprintf(errtxt, "smspec *spec (negative error(s)) [warning]");
        err = cpl_error_set_message(cpl_func, CPL_ERROR_ACCESS_OUT_OF_RANGE,
                                    "%s", errtxt);
    }

    return err;
}


cpl_error_code sm_comp_readlibstruct(cpl_parameterlist *libfilepar,
                                     const smparmodel modelpar,
                                     const char *spectype)
{
    /*!
     * Builds file name of library spectrum or line list depending on input
     * parameters such as library path, radiative transfer code, target
     * altitude or airmass, period of the year, period of the night,
     * precipitable water vapour, resolution, or file type. The required
     * parameters are provided by a ::smparmodel structure except for the file
     * type which has to be given directly and can be either "R" for radiance
     * spectrum, "T" for transmission spectrum, or "L" for list of line
     * transmissions. The file name format is given by a library structure
     * file. The routine returns a CPL parameter list with the file name and
     * the values of all parameters building the file name plus a parameter
     * for vacuum/air wavelengths. The latter is added to indicate whether the
     * library spectrum (always in vacuum wavelengths) has to be converted to
     * air wavelengths.
     *
     * \b INPUT:
     * \param modelpar    sky emission parameters (see typedef of
     *                    ::smparmodel)
     * \param spectype    radiance ("R"), transmission ("T"), or line
     *                    transmission ("L")?
     *
     * \b OUTPUT:
     * \param libfilepar  CPL parameter list containing file name (+ path) of
     *                    library spectrum and its constituents
     *
     * \b ERRORS:
     * - FOF: File opening failed
     * - UFS: Unexpected file structure
     * - IIP: Invalid input parameter(s)
     */

    FILE *stream;
    cpl_parameter *p;
    smparam nameform[SM_MAXPAR], pardef[SM_MAXPAR], val[SM_MAXPAR];
    char errtxt[SM_MAXLEN+1], libstructfile[FILENAME_MAX];
    char filename[FILENAME_MAX], libfile[FILENAME_MAX];
    int i = 0, ascii0 = 48, j = 0, fpos[7] = {-1,-1,-1,-1,-1,-1,-1};
    int flen[7] = {0}, npar = 0, k0 = -1, k = 0;
    double alt = 0., prev = 0., mean = 0., airmass = 0., pwv = 0., resol = 0.;

    /* Check existence of library structure file */
    sprintf(libstructfile, "%s/%s", modelpar.libpath, modelpar.libstruct);
    if ((stream = fopen(libstructfile, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", SM_ERROR_FOF_TXT, libstructfile);
        return cpl_error_set_message(cpl_func, SM_ERROR_FOF, "%s", errtxt);
    }

    /* Read structure of file names */
    sm_param_read(stream, nameform);
    strcpy(filename, nameform[0].c);

    /* Find positions and lengths of the parameter fields in the file name */

    while (filename[i] != '\0' && i < SM_LENLINE) {
        if (filename[i] - ascii0 == j + 1 && fpos[j] < 0) {
            fpos[j] = i;
        }
        if (filename[i] - ascii0 != j + 1 && fpos[j] >= 0) {
            flen[j] = i - fpos[j];
            j++;
            if (j >= 7) {
                break;
            }
            if (filename[i] - ascii0 == j + 1) {
                fpos[j] = i;
            }
        }
        i++;
    }

    /* Number of parameters found */
    npar = j;

    /* 1 to 7 library parameters are allowed */
    if (npar < 1 || npar > 7) {
        fclose(stream);
        sprintf(errtxt, "%s: %s (N_fields < 1 or > 7)", SM_ERROR_UFS_TXT,
                libstructfile);
        return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s", errtxt);
    }

    /* Read parameter format and values and write fitting values in file name
       template */

    for (j = 0; j < npar; j++) {

        /* Parameter format */
        sm_param_read(stream, pardef);

        /* Check parameter description and parameter ID */
        if (pardef[0].n != 3 || pardef[0].i != j + 1) {
            fclose(stream);
            if (pardef[0].n != 3) {
                sprintf(errtxt, "%s: %s (# of parameter description values "
                        "!= 3)", SM_ERROR_UFS_TXT, libstructfile);
            } else {
                sprintf(errtxt, "%s: %s (ID != sequence number)",
                        SM_ERROR_UFS_TXT, libstructfile);
            }
            return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s",
                                         errtxt);
        }

        /* Parameter values */
        sm_param_read(stream, val);

        k0 = -1;

        if (strncmp(pardef[1].c, "rtcode", 6) == 0) {

            /* Valid radiative transfer code? */
            if (modelpar.rtcode[0] != 'L' && modelpar.rtcode[0] != 'R') {
                fclose(stream);
                sprintf(errtxt,
                        "%s: rtcode of smparmodel modelpar != L and != R",
                        SM_ERROR_IIP_TXT);
                return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s",
                                             errtxt);
            } else if (modelpar.rtcode[0] == 'R') {
                cpl_msg_warning(cpl_func, "RFM is no more supported");
            }

            /* Find label of selected radiative transfer code */
            for (k = 0; k < val[0].n; k++) {
                if (val[k].c[0] == modelpar.rtcode[0]) {
                    k0 = k;
                }
            }

            /* Add parameter to list */
            p = cpl_parameter_new_value("rtcode", CPL_TYPE_STRING, "", "",
                                        val[k0].c);
            cpl_parameterlist_append(libfilepar, p);

        } else if (strncmp(pardef[1].c, "season", 6) == 0) {

            /* Valid seasonal bin? */
            if (modelpar.season < 0 || modelpar.season > 6) {
                fclose(stream);
                sprintf(errtxt,
                        "%s: season of smparmodel modelpar < 0 or > 6",
                        SM_ERROR_IIP_TXT);
                return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s",
                                             errtxt);
            }

            /* Find season */
            for (k = 0; k < val[0].n; k++) {
                if (val[k].i == modelpar.season) {
                    k0 = k;
                }
            }

            /* Add parameter to list */
            p = cpl_parameter_new_value("season", CPL_TYPE_INT, "", "",
                                        val[k0].i);
            cpl_parameterlist_append(libfilepar, p);

        } else if (strncmp(pardef[1].c, "time", 4) == 0) {

            /* Valid time bin? */
            if (modelpar.time < 0 || modelpar.time > 3) {
                fclose(stream);
                sprintf(errtxt,
                        "%s: time of smparmodel modelpar < 0 or > 3",
                        SM_ERROR_IIP_TXT);
                return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s",
                                             errtxt);
            }

            /* Find period of the night */
            for (k = 0; k < val[0].n; k++) {
                if (val[k].i == modelpar.time) {
                    k0 = k;
                }
            }

            /* Add parameter to list */
            p = cpl_parameter_new_value("time", CPL_TYPE_INT, "", "",
                                        val[k0].i);
            cpl_parameterlist_append(libfilepar, p);

        } else if (strncmp(pardef[1].c, "alt", 3) == 0) {

            /* Valid target altitude? */
            if (modelpar.alt < 0 || modelpar.alt > 90) {
                fclose(stream);
                sprintf(errtxt,
                        "%s: alt of smparmodel modelpar < 0 or > 90 degrees",
                        SM_ERROR_IIP_TXT);
                return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s",
                                             errtxt);
            }

            /* Find listed altitude with best correspondence to target
               altitude */

            prev = 0.;

            for (k = 0; k < val[0].n; k++) {
                alt = val[k].d * pow(10., pardef[2].d);
                if (alt >= modelpar.alt) {
                    break;
                } else {
                    prev = alt;
                }
            }

            if (k == val[0].n) {
                k0 = k - 1;
            } else {
                mean = (prev + alt) / 2;
                if (k > 0 && mean > modelpar.alt) {
                    k0 = k - 1;
                } else {
                    k0 = k;
                }
            }

            /* Add parameter to list */
            p = cpl_parameter_new_value("alt", CPL_TYPE_DOUBLE, "", "",
                                        val[k0].d * pow(10., pardef[2].d));
            cpl_parameterlist_append(libfilepar, p);

        } else if (strncmp(pardef[1].c, "airmass", 7) == 0) {

            /* Valid target altitude? */
            if (modelpar.alt < 0 || modelpar.alt > 90) {
                fclose(stream);
                sprintf(errtxt,
                        "%s: alt of smparmodel modelpar < 0 or > 90 degrees",
                        SM_ERROR_IIP_TXT);
                return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s",
                                             errtxt);
            }

            /* Find airmass with best correspondence to target altitude */

            prev = 0.;

            for (k = 0; k < val[0].n; k++) {
                airmass = val[k].d * pow(10., pardef[2].d);
                if (airmass >= modelpar.airmass) {
                    break;
                } else {
                    prev = airmass;
                }
            }

            if (k == val[0].n) {
                k0 = k - 1;
            } else {
                mean = (prev + airmass) / 2;
                if (k > 0 && mean >= modelpar.airmass) {
                    k0 = k - 1;
                } else {
                    k0 = k;
                }
            }

            /* Add parameter to list */
            p = cpl_parameter_new_value("airmass", CPL_TYPE_DOUBLE, "", "",
                                        val[k0].d * pow(10., pardef[2].d));
            cpl_parameterlist_append(libfilepar, p);

        } else if (strncmp(pardef[1].c, "pwv", 3) == 0) {

            /* Find listed PWV value with best correspondence to given
               value */

            prev = 0.;

            for (k = 0; k < val[0].n; k++) {
                pwv = val[k].d * pow(10., pardef[2].d);
                if (pwv >= modelpar.pwv) {
                    break;
                } else {
                    prev = pwv;
                }
            }

            if (k == val[0].n) {
                k0 = k - 1;
            } else {
                mean = (prev + pwv) / 2;
                if (k > 0 && mean > modelpar.pwv) {
                    k0 = k - 1;
                } else {
                    k0 = k;
                }
            }

            /* Add parameter to list */
            p = cpl_parameter_new_value("pwv", CPL_TYPE_DOUBLE, "", "",
                                        val[k0].d * pow(10., pardef[2].d));
            cpl_parameterlist_append(libfilepar, p);

        } else if (strncmp(pardef[1].c, "resol", 5) == 0) {

            /* Valid resolution? */
            if (modelpar.resol <= 0) {
                fclose(stream);
                sprintf(errtxt,
                        "%s: resol of smparmodel modelpar <= 0",
                        SM_ERROR_IIP_TXT);
                return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s",
                                             errtxt);
            }

            /* Find closest listed resolution which is at least as high as the
               given value (if not available, take highest value in list) */

            for (k = 0; k < val[0].n; k++) {
                resol = val[k].d * pow(10., pardef[2].d);
                if ((1 + SM_TOL) * resol >= modelpar.resol) {
                    break;
                }
            }

            if (k == val[0].n) {
                k0 = k - 1;
            } else {
                k0 = k;
            }

            /* Add parameter to list */
            p = cpl_parameter_new_value("resol", CPL_TYPE_DOUBLE, "", "",
                                        val[k0].d * pow(10., pardef[2].d));
            cpl_parameterlist_append(libfilepar, p);

        } else if (strncmp(pardef[1].c, "spectype", 8) == 0) {

            /* Radiance or transmission? */
            for (k = 0; k < val[0].n; k++) {
                if (val[k].c[0] == spectype[0]) {
                    k0 = k;
                }
            }

            /* Add parameter to list */
            p = cpl_parameter_new_value("spectype", CPL_TYPE_STRING, "", "",
                                        val[k0].c);
            cpl_parameterlist_append(libfilepar, p);

        } else {

            /* Unknown parameter label */
            fclose(stream);
            sprintf(errtxt, "%s: %s (parameter labels 'rtcode', 'season', "
                    "'time', 'alt', 'airmass', 'pwv', 'resol', or 'spectype' "
                    "only)", SM_ERROR_UFS_TXT, libstructfile);
            return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s",
                                         errtxt);

        }

        /* Error if the selected value is not present in the library
           structure file */
        if (k0 == -1) {
            fclose(stream);
            sprintf(errtxt, "%s: %s (library structure inconsistent with "
                    "selected parameter values)", SM_ERROR_UFS_TXT,
                    libstructfile);
            return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s",
                                         errtxt);
        }

        /* Copy best-fit parameter value to file name template */
        for (i = 0; i < flen[j]; i++) {
            filename[i+fpos[j]] = val[k0].c[i];
        }

    }

    /* Close file stream */
    fclose(stream);

    /* Compose name of LBLRTM/RFM library spectrum */
    sprintf(libfile, "%s/%s", modelpar.libpath, filename);

    /* Add file name to parameter list */
    p = cpl_parameter_new_value("filename", CPL_TYPE_STRING, "", "", libfile);
    cpl_parameterlist_append(libfilepar, p);

    /* Add vacuum/air wavelength parameter to parameter list */
    if (strncmp(modelpar.vac_air, "air", 3) == 0) {
        p = cpl_parameter_new_value("vac_air", CPL_TYPE_STRING, "", "",
                                    "air");
    } else {
        p = cpl_parameter_new_value("vac_air", CPL_TYPE_STRING, "", "",
                                    "vac");
    }
    cpl_parameterlist_append(libfilepar, p);

    return CPL_ERROR_NONE;
}


double sm_comp_vactoair_single(const double lam)
{
    /*!
     * Converts vacuum wavelength [\f$\mu{\rm m}\f$] to air wavelength by
     * using the formula of Edlen (1966)
     *
     * \b INPUT:
     * \param lam  vacuum wavelength in \f$\mu{\rm m}\f$
     *
     * \b RETURN:
     * - air wavelength in \f$\mu{\rm m}\f$
     *
     * \b ERRORS:
     * - none
     */

    double sig2 = 0., n = 0., airlam = 0.;

    sig2 = pow(lam, -2);
    n = 8342.13 + 2406030. / (130. - sig2) + 15997. / (38.9 - sig2);
    n = 1. + 1e-8 * n;
    airlam = lam / n;

    return airlam;
}


cpl_error_code sm_comp_vactoair_spec(smspec *spec)
{
    /*!
     * Converts spectrum of vacuum wavelengths to spectrum of air wavelengths
     * by using the formula of Edlen (1966)
     *
     * \b INPUT:
     * \param spec  spectrum with vacuum wavelengths
     *
     * \b OUTPUT:
     * \param spec  spectrum with air wavelengths
     *
     * \b ERRORS:
     * - none
     */

    int i = 0;
    double lam = 0.;

    for (i = 0; i < spec->n; i++) {
        lam = sm_comp_vactoair_single(spec->dat[i].lam);
        spec->dat[i].lam = lam;
    }

    return CPL_ERROR_NONE;
}


double sm_comp_airtovac_single(const double lam)
{
    /*!
     * Converts air wavelength [\f$\mu{\rm m}\f$] to vacuum wavelength by
     * using the formula of Edlen (1966)
     *
     * \b INPUT:
     * \param lam  air wavelength in \f$\mu{\rm m}\f$
     *
     * \b RETURN:
     * - vacuum wavelength in \f$\mu{\rm m}\f$
     *
     * \b ERRORS:
     * - none
     */

    double sig2 = 0., n = 0., airlam = 0.;

    sig2 = pow(lam, -2);
    n = 8342.13 + 2406030. / (130. - sig2) + 15997. / (38.9 - sig2);
    n = 1. + 1e-8 * n;
    airlam = lam * n;

    return airlam;
}


cpl_error_code sm_comp_airtovac_spec(smspec *spec)
{
    /*!
     * Converts spectrum of air wavelengths to spectrum of vacuum wavelengths
     * by using the formula of Edlen (1966)
     *
     * \b INPUT:
     * \param spec  spectrum with air wavelengths
     *
     * \b OUTPUT:
     * \param spec  spectrum with vacuum wavelengths
     *
     * \b ERRORS:
     * - none
     */

    int i = 0;
    double lam = 0.;

    for (i = 0; i < spec->n; i++) {
        lam = sm_comp_airtovac_single(spec->dat[i].lam);
        spec->dat[i].lam = lam;
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_comp_scaleo3(smspec *molabstrans, smspec *abstrans,
                               smspec *o3trans, const smspec *libtrans,
                               const smparmodel modelpar)
{
    /*!
     * Separates UV/optical ozone absorption from molecular absorption (as
     * calculated by the radiative transfer code) and scales the ozone
     * transmission curve using the given relative ozone column density. The
     * result is then combined again with the separated transmission curve
     * for all other absorption.
     *
     * \b INPUT:
     * \param libtrans     transmission curve from radiative transfer code
     *                     (scaled to airmass = 1)
     * \param modelpar     sky emission parameters
     *                     (see typedef of ::smparmodel)
     *
     * \b OUTPUT:
     * \param molabstrans  corrected transmission curve from radiative
     *                     transfer code
     * \param abstrans     transmission curve from radiative transfer code
     *                     without ozone UV/optical absorption
     * \param o3trans      transmission curve for ozone UV/optical absorption
     *
     * \b ERRORS:
     * - none
     */

    smspec ino3trans;
    char o3transfile[FILENAME_MAX];

    /* Read transmission curve for UV/optical ozone absorption
       (airmass = 1) */
    sprintf(o3transfile, "%s/%s", modelpar.datapath, modelpar.o3transname);
    sm_spec_read(&ino3trans, o3transfile);

    /* Convert vacuum to air wavelengths if required */
    if (strncmp(modelpar.vac_air, "air", 3) == 0) {
        sm_comp_vactoair_spec(&ino3trans);
    }

    /* Wavelength grid for ozone transmission curve */
    sm_spec_copy(o3trans, libtrans);
    o3trans->type = 2;

    /* Adapt ozone transmission curve to output wavelength grid */
    sm_spec_interpol(o3trans, &ino3trans);

    /* Divide transmission curve for molecular absorption (scaled to
       airmass = 1) by UV/optical ozone absorption */
    sm_spec_copy(abstrans, libtrans);
    abstrans->type = 2;
    sm_spec_calc(abstrans, '/', o3trans);

    /* Eliminate numerical noise in the bluest end of transmission
       and force trans=1 when trans>1 everywhere */
    for (int index = 0; index < abstrans->n; ++index) {
        if (abstrans->dat[index].lam < 0.360 || abstrans->dat[index].flux > 1.) {
            abstrans->dat[index].flux = 1.;
        }
    }

    /* Change ozone column density */
    sm_spec_scale(o3trans, '^', modelpar.o3column);

    /* Recombine molecular and ozone absorption and curve */
    sm_spec_copy(molabstrans, abstrans);
    sm_spec_calc(molabstrans, '*', o3trans);

    /* Free memory */
    sm_spec_free(&ino3trans);

    return CPL_ERROR_NONE;
}


double sm_comp_calcrayleighscat1(const double lam, const double pres, const double sm_h)
{
    /*!
     * Calculates optical depth by Rayleigh scattering for given wavelength.
     * The parametrisation of Liou (2002, P. 352). The pressure for Cerro
     * Paranal is an input parameter.
     *
     * \b INPUT:
     * \param lam   wavelength in \f$\mu{\rm m}\f$
     * \param pres  atmospheric pressure at Cerro Paranal in hPa
     *
     * \b RETURN:
     * - optical depth at zenith for given wavelength
     *
     * \b ERRORS:
     * - none
     */

    /* Fixed scattering curve parameters */
    double a = 0.00864, b = 6.5e-6, c = 3.916, d = 0.074, e = 0.050;
    /* Paranal-related parameters for pressure profile */
    double ps = 1013.25; // standard pressure at sea level in mbar

    /* Calculate optical depth at zenith for given wavelength */
    return (a + b * sm_h) * pow(lam, -c - d * lam - e / lam) * pres / ps;
}


cpl_error_code sm_comp_calcrayleighscat(smspec *rscattrans,
                                        const smparmodel modelpar)
{
    /*!
     * Calculates Rayleigh scattering transmission curve by using the
     * parametrisation given by Liou (2002, P. 352). The pressure for Cerro
     * Paranal is taken from the ::smparmodel structure.
     *
     * \b INPUT:
     * \param rscattrans  desired wavelength grid in \f$\mu{\rm m}\f$
     * \param modelpar    sky emission parameters (see typedef of
     *                    ::smparmodel)
     *
     * \b OUTPUT:
     * \param rscattrans  Rayleigh scattering curve (transmission at zenith)
     *
     * \b ERRORS:
     * - none
     */

    int i = 0;
    double taur = 0.;

    /* Calculate transmission at zenith for each wavelength */
    for (i = 0; i < rscattrans->n; i++) {
        taur = sm_comp_calcrayleighscat1(rscattrans->dat[i].lam,
                                         modelpar.pres,
                                         modelpar.sm_h);
        rscattrans->dat[i].flux = exp(-taur);
    }

    return CPL_ERROR_NONE;
}


double sm_comp_calcmiescat1(const double lam)
{
    /*!
     * Calculates optical depth my Mie (aerosol) extinction for given
     * wavelength. The parametrisation of Patat (2011) is used and a cut in
     * the optical depth for wavelengths below \f$0.4\,\mu{\rm m}\f$ is
     * applied.
     *
     * \b INPUT:
     * \param lam  wavelength in \f$\mu{\rm m}\f$
     *
     * \b RETURN:
     * - optical depth at zenith for given wavelength
     *
     * \b ERRORS:
     * - none
     */

    double k0 = 0., lamx = 0.;
    double lam0 = 0.4, k0mag = 0.014, alpha = -1.38;

    /* Convert extinction coefficient */
    k0 = 0.4 * log(10) * k0mag;

    /* Consider cut in optical depth */
    if (lam < lam0) {
        lamx = lam0;
    } else {
        lamx = lam;
    }

    /* Calculate optical depth at zenith for given wavelength */
    return k0 * pow(lamx, alpha);
}


cpl_error_code sm_comp_calcmiescat(smspec *mscattrans,
                                   const smparmodel modelpar)
{
    /*!
     * Calculates transmission curve for Mie (aerosol) extinction depending
     * on the "mieextname" parameter of the ::smparmodel structure. The
     * string "NONE" selects the parametrisation given by Patat et al. (2011)
     * and a cut in the optical depth for wavelengths below
     * \f$0.4\,\mu{\rm m}\f$. For all other names, the corresponding aerosol
     * extinction table in the "data" folder is loaded. The file has to
     * provide wavelengths in \f$\mu{\rm m}\f$ and optical depths in natural
     * units.
     *
     * \b INPUT:
     * \param mscattrans  desired wavelength grid in \f$\mu{\rm m}\f$
     * \param modelpar    sky emission parameters (see typedef of
     *                    ::smparmodel)
     *
     * \b OUTPUT:
     * \param mscattrans  Mie extinction curve (transmission at zenith)
     *
     * \b ERRORS:
     * - none
     */

    smspec mieext;
    char mieextfile[FILENAME_MAX];
    int i = 0;
    double taum = 0.;

    /* Calculate transmission at zenith for each wavelength */

    if (strncmp(modelpar.mieextname, "NONE", 4) == 0) {

        /* Use standard parameterisation by Patat et al. (2011) */
        for (i = 0; i < mscattrans->n; i++) {
            taum = sm_comp_calcmiescat1(mscattrans->dat[i].lam);
            mscattrans->dat[i].flux = exp(-taum);
        }

    } else {

        /* Load extinction table, convert optical depth into transmission,
           and map resulting spectrum to input wavelength grid */
        sprintf(mieextfile, "%s/%s", modelpar.datapath, modelpar.mieextname);
        sm_spec_read(&mieext, mieextfile);
        sm_spec_scale(&mieext, '*', -1.);
        sm_spec_funct(&mieext, "exp", 'e');
        sm_spec_interpol(mscattrans, &mieext);
        sm_spec_free(&mieext);

    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_comp_scaletranscurv(smspec *trans,
                                      const smparmodel modelpar)
{
    /*!
     * Scales transmission curve depending on airmass as calculated by the
     * formula of Rozenberg (1966)
     *
     * \b INPUT:
     * \param trans     transmission curve at airmass = 1
     * \param modelpar  sky emission parameters (see typedef of ::smparmodel)
     *
     * \b OUTPUT:
     * \param trans     transmission curve for given altitude
     *
     * \b ERRORS:
     * - NDA: No data
     * - IIP: Invalid input parameter(s)
     */

    char errtxt[SM_MAXLEN+1];

    /* No data points */
    if (trans->n <= 0) {
        sprintf(errtxt, "%s: smspec *trans", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    }

    /* Check validity of object altitude and, hence, airmass */
    if (modelpar.alt < 0. || modelpar.alt > 90.) {
        sprintf(errtxt, "%s: alt of smparmodel modelpar < 0 or > 90 degrees",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Airmass correction of transmission curve */
    sm_spec_scale(trans, '^', modelpar.airmass);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_comp_writetranscomp(cpl_table *transcomp,
                                      const char *colname,
                                      const smspec *trans,
                                      const smparmodel modelpar)
{
    /*!
     * Writes input transmission curve (airmass = 1 is required) for given
     * object airmass into a CPL table for transmission components. The table
     * column is an input parameter. It will be created if it does not exist.
     *
     * \b INPUT:
     * \param transcomp  CPL table for transmission components of sky model
     * \param colname    name of CPL table column
     * \param trans      input transmission curve as ::smspec structure
     * \param modelpar   sky emission parameters (see typedef of ::smparmodel)
     *
     * \b OUTPUT:
     * \param transcomp  component table with added transmission data
     *
     * \b ERRORS:
     * - none
     */

    smspec tmptrans;

    /* Copy transmission spectrum */
    sm_spec_copy(&tmptrans, trans);

    /* Scale transmission curve (airmass = 1) to object airmass */
    sm_comp_scaletranscurv(&tmptrans, modelpar);

    /* Write transmission values into desired column of CPL table */
    sm_spec_writecplcolumn(transcomp, colname, &tmptrans, 2);

    /* Free memory */
    sm_spec_free(&tmptrans);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_comp_lunskybright_o(smspec *spec, const smparmodel modelpar,
                                      const smspec *solspec,
                                      const smspec *molabstrans,
                                      const smspec *rscattrans,
                                      const smspec *mscattrans)
{
    /*!
     * \callgraph
     * Evaluates predicted lunar part of sky brightness following
     * Krisciunas & Schaefer (1991, PASP 103, 1033).
     *
     * Original implementation by J. Thorstensen.
     * Modified by F. Patat to generalise to UBVRI filters.
     * Spectroscopic version by S. Noll.
     *
     * This routine uses the solar spectrum as a proxy for the unattenuated
     * lunar spectrum and two different transmission curves: one with and one
     * wihout molecular absorption. The second one is important for the
     * estimate of the amount of scattered moonlight at the object position
     * and includes Rayleigh scattering and aerosol extinction. The model
     * depends on object and moon altitude, their angular separation, moon
     * phase (expressed by angular distance of moon and sun), and moon
     * distance.
     *
     * \note
     * - The model becomes unreliable for moon-object distances around and
     *   below 30°. There, the main contribution comes from aerosol
     *   scattering, which is highly variable and, therefore, difficult to
     *   predict.
     * - The routine uses the library transmission curve selected
     *   depending on the object altitude and not the moon altitude (see
     *   ::sm_comp_getmolspec). This reduces the computing time significantly,
     *   but especially the flux in centres of telluric absorption features
     *   can be overestimated if the moon is much lower than the object.
     *   However, for the composite sky emission spectrum the flux differences
     *   should be lower than 10\% even in the A-band.
     *
     * \b INPUT:
     * \param spec         desired wavelength grid
     * \param modelpar     sky emission parameters (see typedef of
     *                     ::smparmodel)
     * \param solspec      solar spectrum
     * \param molabstrans  transmission curve from radiative transfer code
     * \param rscattrans   transmission curve for Rayleigh scattering
     * \param mscattrans   transmission curve for Mie scattering
     *
     * \b OUTPUT:
     * \param spec         spectrum of lunar sky brightness
     *
     * \b ERRORS:
     * - NDA: No data
     * - IIP: Invalid input parameter(s)
     * - IDG: Inconsistent data grids
     * - see subroutines as well
     */

    smspec trans, iscatr, iscatm;
    char errtxt[SM_MAXLEN+1];
    int i = 0;
    double phi = 0., rho = 0., rho_rad = 0., zmoon = 0., z = 0.;
    double fr = 0., fm = 0., xzm = 0., xo = 0.;

    /* Scale factor for solar spectrum (Colina et al. 1996)
       -> flux in V-band */
    double scale_solspec = 1846.;
    /* Flux conversion: W m^-2 mum^-1 arcsec^-2 * lam [mum]
                     -> phot s^-1 m^-2 mum^-1 arcsec^-2 */
    double conv = SM_LAM_UNIT / (CPL_PHYS_H * CPL_PHYS_C);

    /* No error treatment */
    spec->type = 2;

    /* spec = NULL */
    if (spec->n <= 0) {
        sprintf(errtxt, "%s: smspec *spec", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    }

    /* Check validity of input parameters */
    if (modelpar.alpha < 0. || modelpar.alpha > 360. ||
        modelpar.rho < 0. || modelpar.rho > 180. ||
        modelpar.altmoon < -90. || modelpar.altmoon > 90. ||
        modelpar.alt < 0. || modelpar.alt > 90. ||
        modelpar.moondist < 0.91 || modelpar.moondist > 1.08) {
        sprintf(errtxt, "%s: alpha, rho, altmoon, alt, and/or moondist of "
                "smparmodel modelpar out of range", SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Conversions */
    phi = fabs(180. - modelpar.alpha);
    rho = modelpar.rho;
    rho_rad = modelpar.rho * CPL_MATH_RAD_DEG;
    zmoon = (90. - modelpar.altmoon) * CPL_MATH_RAD_DEG;
    z = (90. - modelpar.alt) * CPL_MATH_RAD_DEG;

    /* Test agreement of wavelength grids in input spectra */
    if ((int) sm_spec_compgrids(solspec, spec) == (int) SM_ERROR_IDG ||
        (int) sm_spec_compgrids(molabstrans, spec) == (int) SM_ERROR_IDG ||
        (int) sm_spec_compgrids(rscattrans, spec) == (int) SM_ERROR_IDG ||
        (int) sm_spec_compgrids(mscattrans, spec) == (int) SM_ERROR_IDG) {
        sprintf(errtxt, "%s: smspec *solspec, *molabstrans, *rscattrans, or "
                "*mscattrans != smspec *spec", SM_ERROR_IDG_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IDG, "%s", errtxt);
    }

    /* Scale solar spectrum to full moon illuminance
       (using Schaefer's S&T paper constants and V-band) */
    sm_spec_calc(spec, '=', solspec);
    sm_spec_scale(spec, '*', 2.31e-20 / scale_solspec);

    /* Compute apparent magnitude of moon */
    sm_spec_scale(spec, '*', pow(10., -0.4 * (0.026 * phi +
                                              4.0e-9 * pow(phi, 4.))) /
                                      pow(modelpar.moondist, 2));

    /* Correct for opposition effect:
       35% brighter at full, effect tapering linearly to
       zero at 7 degrees away from full. Mentioned peripherically in
       Krisciunas and Schaefer, p. 1035. */
    if (phi < 7.) {
        sm_spec_scale(spec, '*', 1.35 - 0.05 * phi);
    }

    /* Compute Rayleigh scattering function */
    fr = 229087. * (1.06 + pow(cos(rho_rad), 2));

    /* Compute Mie scattering function */
    if (rho >= 10.) {
        fm = pow(10., 6.15 - rho / 40.); /* eqn 21 */
    } else if (rho < 10. && rho > 0.25) {
        fm = 6.2e7 / pow(rho, 2); /* eqn 19 */
    } else {
        fm = 9.9e8; /* for 1/4 degree -- radius of moon! */
    }

    /* Compute moon airmass */
    if (modelpar.altmoon <= 0.) {
        xzm = 1e4;
    } else {
        xzm = 1 / sqrt(1.0 - 0.96 * pow(sin(zmoon), 2));
    }

    /* Compute target airmass */
    xo = 1 / sqrt(1.0 - 0.96 * pow(sin(z), 2));

    /* Compute relative scattering intensities */
    sm_spec_copy(&iscatr, rscattrans);
    sm_spec_copy(&iscatm, mscattrans);
    for (i = 0; i < spec->n; i++) {
        /* Rayleigh scattering */
        iscatr.dat[i].flux = (1. - pow(rscattrans->dat[i].flux, xo));
        /* Mie scattering */
        iscatm.dat[i].flux = (1. - pow(mscattrans->dat[i].flux, xo));
    }

    /* Get absolute scattering intensities */
    sm_spec_scale(&iscatr, '*', 2.2 * fr);
    sm_spec_scale(&iscatm, '*', 10. * fm);

    /* Get total extinction curve */
    sm_spec_copy(&trans, rscattrans);
    sm_spec_calc(&trans, '*', mscattrans);
    sm_spec_calc(&trans, '*', molabstrans);

    /* Compute Bmoon (in W m-2 mum-1 arcsec-2) */
    for (i = 0; i < spec->n; i++) {
        spec->dat[i].flux *= (iscatr.dat[i].flux + iscatm.dat[i].flux)
                             * pow(trans.dat[i].flux, xzm);
    }

    /* Flux in phot s^-1 m^-2 mum^-1 arcsec^-2 */
    sm_spec_convunits(spec, conv, 1);

    /* Free memory */
    sm_spec_free(&trans);
    sm_spec_free(&iscatr);
    sm_spec_free(&iscatm);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_comp_lunskybright(smspec *spec, const smparmodel modelpar,
                                    const smspec *solspec,
                                    const smspec *abstrans,
                                    const smspec *o3trans,
                                    const smspec *rscattrans,
                                    const smspec *mscattrans)
{
    /*!
     * \callgraph
     * Calculates scattered moonlight. The Moon illuminance is approximated
     * from Kieffer & Stone (2005) based on ROLO data. The scattering in the
     * atmosphere is obtained from single (plus optional double) scattering
     * calculations with multiple scattering correction (see ::sm_scat_moon
     * and subroutines) and the Cerro Paranal extinction curve of Patat et al.
     * (2011). A simple estimate which depends on the optical depth dependent
     * effective airmass is taken for the molecular absorption of the
     * moonlight. The absorption by the stratospheric ozone layer is
     * considered separately. The model depends on object and Moon altitude,
     * their angular separation, Moon phase (expressed by angular distance of
     * Moon and Sun), and Moon distance.
     *
     * \b INPUT:
     * \param spec        desired wavelength grid
     * \param modelpar    sky emission parameters
     *                    (see typedef of ::smparmodel)
     * \param solspec     solar spectrum
     * \param abstrans    transmission curve from radiative transfer code
     *                    without ozone UV/optical absorption
     * \param o3trans     transmission curve for ozone UV/optical absorption
     * \param rscattrans  transmission curve for Rayleigh scattering
     * \param mscattrans  transmission curve for Mie scattering
     *
     * \b OUTPUT:
     * \param spec        spectrum of lunar sky brightness
     *
     * \b ERRORS:
     * - NDA: No data
     * - IIP: Invalid input parameter(s)
     * - IDG: Inconsistent data grids
     * - see subroutines as well
     */

    smgrid miephase, sscatcor;
    smspec albedo, iscat, o3transv, dscatcor;
    char errtxt[SM_MAXLEN+1], miephasefile[FILENAME_MAX];
    char sscatcorfile[FILENAME_MAX], dscatcorfile[FILENAME_MAX];
    double rho = 0., zmoon = 0., z = 0., xo3 = 0.;

    /* Solid angle of the Moon in sr */
    double omegamoon = 6.4177e-5;
    /* Flux conversion: W m^-2 mum^-1 arcsec^-2 * lam [mum]
                     -> phot s^-1 m^-2 mum^-1 arcsec^-2 */
    double conv = SM_LAM_UNIT / (CPL_PHYS_H * CPL_PHYS_C);

    /* No error treatment */
    spec->type = 2;

    /* spec = NULL */
    if (spec->n <= 0) {
        sprintf(errtxt, "%s: smspec *spec", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    }

    /* Check validity of input parameters */
    if (modelpar.alpha < 0. || modelpar.alpha > 360. ||
        modelpar.rho < 0. || modelpar.rho > 180. ||
        modelpar.altmoon < -90. || modelpar.altmoon > 90. ||
        modelpar.alt < 0. || modelpar.alt > 90. ||
        modelpar.moondist < 0.91 || modelpar.moondist > 1.08 ||
        modelpar.ssa < 0. || modelpar.ssa > 1.) {
        sprintf(errtxt, "%s: alpha, rho, altmoon, alt, moondist, and/or ssa "
                "of smparmodel modelpar out of range", SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Check validity of calcds string */
    if (modelpar.calcds[0] != 'Y' && modelpar.calcds[0] != 'y' &&
        modelpar.calcds[0] != 'N' && modelpar.calcds[0] != 'n') {
        sprintf(errtxt, "%s: calcds of smparmodel != 'Y' and 'N'",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Check validity of sm_h */
    if (modelpar.sm_h < 2.4 || modelpar.sm_h > 3.060) { // lasilla - armazones
        sprintf(errtxt, "%s: sm_h of smparmodel is out of range",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Check validity of sm_hmin */
    if (modelpar.sm_hmin < 2.0 || modelpar.sm_h > 3.060) { // lasilla - armazones
        sprintf(errtxt, "%s: sm_hmin of smparmodel is out of range",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Skip calculations in the case of dark time */
    if (modelpar.altmoon < -18.) {
        return CPL_ERROR_NONE;
    }

    /* Conversions */
    rho = modelpar.rho;
    zmoon = 90. - modelpar.altmoon;
    z = 90. - modelpar.alt;

    /* Test agreement of wavelength grids in input spectra */
    if ((int) sm_spec_compgrids(solspec, spec) == (int) SM_ERROR_IDG ||
        (int) sm_spec_compgrids(abstrans, spec) == (int) SM_ERROR_IDG ||
        (int) sm_spec_compgrids(o3trans, spec) == (int) SM_ERROR_IDG ||
        (int) sm_spec_compgrids(rscattrans, spec) == (int) SM_ERROR_IDG ||
        (int) sm_spec_compgrids(mscattrans, spec) == (int) SM_ERROR_IDG) {
        sprintf(errtxt, "%s: smspec *solspec, *abstrans, *o3trans, "
                "*rscattrans, or *mscattrans != smspec *spec",
                SM_ERROR_IDG_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IDG, "%s", errtxt);
    }

    /* Wavelength grid for lunar albedo spectrum */
    sm_spec_copy(&albedo, spec);

    /* Wavelength grid for scattering intensity curve */
    sm_spec_copy(&iscat, spec);

    /* Take solar spectrum from Colina et al. (1996) in W m^-2 mum^-1 */
    sm_spec_calc(spec, '=', solspec);

    /* Compute lunar irradiance according to Kieffer & Stone (2005) model */
    sm_comp_getmoonalbedo(&albedo, modelpar);
    sm_spec_calc(spec, '*', &albedo);
    sm_spec_scale(spec, '*', omegamoon / CPL_MATH_PI);
    sm_spec_scale(spec, '/', modelpar.moondist * modelpar.moondist);

    /* Get "airmass" for ozone layer (van Rhijn formula for height of
       25 km) */
    if (modelpar.altmoon < 0.) {
        xo3 = 1. / sqrt(0.008);
    } else {
        xo3 = 1. / sqrt(1.0 - 0.992 * pow(sin(zmoon * CPL_MATH_RAD_DEG), 2));
    }

    /* Copy UV/optical ozone transmission curve */
    sm_spec_copy(&o3transv, o3trans);

    /* Estimate ozone transmission for given airmass */
    sm_spec_scale(&o3transv, '^', xo3);

    /* Correct intensity for ozone absorption */
    sm_spec_calc(spec, '*', &o3transv);

    /* Read table for Mie scattering phase functions */
    sprintf(miephasefile, "%s/%s", modelpar.datapath, modelpar.miephasename);
    sm_grid_read(&miephase, miephasefile);

    /* Read table for multiple scattering correction for single scattering */
    sprintf(sscatcorfile, "%s/%s", modelpar.datapath, modelpar.sscatcorname);
    sm_grid_read(&sscatcor, sscatcorfile);

    /* Read spectrum for multiple scattering correction for double
       scattering */
    sprintf(dscatcorfile, "%s/%s", modelpar.datapath, modelpar.dscatcorname);
    sm_spec_read(&dscatcor, dscatcorfile);

    /* Get scattering intensities from 3D single scattering calculations with
       multiple scattering correction */
    sm_scat_moon(&iscat, modelpar.sm_h, modelpar.sm_hmin, z, zmoon, rho, modelpar.ssa, abstrans, rscattrans,
                 mscattrans, &miephase, &sscatcor, &dscatcor,
                 modelpar.calcds);

    /* Multiply lunar spectrum by scattering curve and convert integrated flux
       into flux per arcsec^-2 */
    sm_spec_calc(spec, '*', &iscat);

    /* Intensity in phot s^-1 m^-2 mum^-1 arcsec^-2 */
    sm_spec_convunits(spec, conv, 1);

    /* Adaption to observed data (FORS sample) */
    sm_spec_scale(spec, '*', modelpar.moonscal);

    /* Free memory */
    sm_spec_free(&albedo);
    sm_spec_free(&iscat);
    sm_spec_free(&o3transv);
    sm_grid_free(&miephase);
    sm_grid_free(&sscatcor);
    sm_spec_free(&dscatcor);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_comp_getmoonalbedo(smspec *albedo,
                                     const smparmodel modelpar)
{
    /*!
     * Calculates the wavelength-dependent disc-equivalent lunar albedo from
     * the Moon illuminance model of Kieffer & Stone (2005) based on ROLO
     * data. Execpt for the libration-dependent terms, the full model is
     * implemented. The required model coefficients are read from a file.
     * As user parameter, the angle between Sun and Moon as seen from Earth is
     * provided by a ::smparmodel structure. For the calculation of the Moon
     * albedo, it is converted into absolute Moon phase angle and
     * selenographic longitude of the Sun. The latter is required for
     * considering differences in the brightness of the waxing and the waning
     * Moon.
     *
     * \b INPUT:
     * \param albedo    desired wavelength grid
     * \param modelpar  sky emission parameters (see typedef of ::smparmodel)
     *
     * \b OUTPUT:
     * \param albedo    spectrum of lunar albedo
     *
     * \b ERRORS:
     * - IIP: Invalid input parameters
     * - FOF: File opening failed
     */

    FILE *stream;
    smspec lunalb;
    smparam p[SM_MAXPAR], n[SM_MAXPAR], c[SM_MAXPAR];
    char errtxt[SM_MAXLEN+1], lunirrfile[FILENAME_MAX];
    int i = 0;
    double gdeg = 0., g = 0., phi0 = 0., phi = 0., aterm = 0., bterm = 0.;
    double dterm = 0.;

    /* Check validity of input angle alpha */
    if (modelpar.alpha < 0. || modelpar.alpha > 360.) {
        sprintf(errtxt, "%s: alpha of smparmodel modelpar out of range",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Moon phase angle (absolute value) */
    gdeg = fabs(180. - modelpar.alpha);
    g = gdeg * CPL_MATH_RAD_DEG;

    /* Estimate of the selenographic longitude of the Sun (decreases during a
       synodic month) */
    phi0 = (180. - modelpar.alpha) * CPL_MATH_RAD_DEG;
    if (gdeg < 97.) {
        phi = phi0;
    } else {
        /* Avoid extrapolation of b-term beyond ROLO data */
        phi = 97. * phi0 / gdeg;
    }

    /* Check existence of Moon albedo data file */
    sprintf(lunirrfile, "%s/%s", modelpar.datapath, modelpar.lunirrname);
    if ((stream = fopen(lunirrfile, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", SM_ERROR_FOF_TXT, lunirrfile);
        return cpl_error_set_message(cpl_func, SM_ERROR_FOF, "%s", errtxt);
    }

    /* Read data file and calculate albedo spectrum */

    /* Read constant coefficients of lunar irradiance model */
    sm_param_read(stream, p);

    /* Get number of wavelengths  */
    sm_param_read(stream, n);

    /* Allocate memory for lunar albedo spectrum */
    sm_spec_malloc(&lunalb, n[0].i);
    if (lunalb.n == 0) {
        fclose(stream);
        sprintf(errtxt, "%s: smspec lunalb", SM_ERROR_ISM_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_ISM, "%s", errtxt);
    }

    for (i = 0; i < n[0].i; i++) {

        /* Read wavelength-dependent coefficients of lunar irradiance model */
        sm_param_read(stream, c);

        /* Convert read wavelength from nm to mum and write it into
           spectrum */
        lunalb.dat[i].lam = 0.001 * c[0].d;

        /* Calculate lunar albedo model and write it into spectrum */
        aterm = c[1].d + c[2].d * g + c[3].d * g * g + c[4].d * pow(g, 3);
        bterm = c[5].d * phi + c[6].d * pow(phi, 3) + c[7].d * pow(phi, 5);
        dterm = c[8].d * exp(-gdeg / p[0].d) + c[9].d * exp(-gdeg / p[1].d) +
                c[10].d * cos((gdeg - p[2].d) / p[3].d);
        lunalb.dat[i].flux = exp(aterm + bterm + dterm);

        /* Old approximation */
        //lunalb.dat[i].flux = lunalb.dat[i].lam * 2.125e-1 *
        //                     exp(-2.832e-2 * gdeg) +
        //                     1.344e-2 * exp(-4.285e-2 * gdeg);

    }

    /* Close data file */
    fclose(stream);

    /* Correction of possible underluminosity of ROLO data
       (Velikodsky et al. 2011) */
    sm_spec_scale(&lunalb, '/', 0.87);

    /* Interpolate tabulated model to input wavelength grid */
    sm_spec_interpol(albedo, &lunalb);

    /* Free memory */
    sm_spec_free(&lunalb);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_comp_scatstarlight(smspec *spec, const smparmodel modelpar,
                                     const smspec *molabstrans)
{
    /*!
     * \callgraph
     * Reads representative spectrum of scattered starlight, rebins it to
     * given wavelength grid by means of linear interpolation, and considers
     * telluric absorption for an effective airmass of 1.25. The spectrum was
     * calculated by means of integrated starlight photometric data based on
     * Pioneer 10/11 observations (see Toller et al. 1987;
     * Leinert et al. 1998; Melchior et al. 2007), spectroscopic data by
     * Mattila (1980), and 3D atmospheric scattering calculations (see
     * Wolstencroft & van Breda 1967; Staude 1975; Bernstein et al. 2002).
     *
     * \b INPUT:
     * \param spec         desired wavelength grid
     * \param modelpar     sky emission parameters (see typedef of
     *                     ::smparmodel)
     * \param molabstrans  transmission curve from radiative transfer code
     *                     (scaled to airmass = 1)
     *
     * \b OUTPUT:
     * \param spec         spectrum of scattered starlight
     *
     * \b ERRORS:
     * - IDG: Inconsistent data grids
     * - see ::sm_spec_rebin
     */

    cpl_error_code err = CPL_ERROR_NONE;
    smspec starspec, atrans;
    char errtxt[SM_MAXLEN+1], starspecfile[FILENAME_MAX];
    double limlam[2] = {0., 0.};
    double airmass = 1.25; // effective airmass for input starlight spectrum

    /* Read mean spectrum of scattered starlight */
    sprintf(starspecfile, "%s/%s", modelpar.datapath, modelpar.starspecname);
    if ((err = sm_spec_read(&starspec, starspecfile)) != CPL_ERROR_NONE) {
        return err;
    }

    /* Test agreement of wavelength grids in input spectra */
    if ((int) sm_spec_compgrids(molabstrans, spec) == (int) SM_ERROR_IDG) {
        sprintf(errtxt, "%s: smspec *molabstrans != smspec *spec",
                SM_ERROR_IDG_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IDG, "%s", errtxt);
    }

    /* Rebin read spectrum to input wavelength grid */
    sm_spec_interpol(spec, &starspec);

    /* Set flux = 0 beyond 2.5 mum */
    limlam[0] = 2.5 + SM_TOL;
    limlam[1] = spec->dat[(spec->n)-1].lam;
    if (limlam[1] >= limlam[0]) {
        sm_spec_scalerange(spec, limlam, '=', 0.);
    }

    /* Correction for molecular absorption of starlight in the lower
       atmosphere */
    sm_spec_copy(&atrans, molabstrans);
    atrans.type = 2;
    sm_spec_scale(&atrans, '^', airmass);
    sm_spec_calc(spec, '*', &atrans);

    /* Free memory */
    sm_spec_free(&starspec);
    sm_spec_free(&atrans);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_comp_zodskybright(smspec *spec, const smparmodel modelpar,
                                    const smspec *solspec,
                                    const smspec *molabstrans,
                                    const smspec *rscattrans,
                                    const smspec *mscattrans)
{
    /*!
     * \callgraph
     * Derives sky brightness caused by the zodiacal light.
     * The procedure is based on recipes provided in Leinert et al. (1998)
     * and uses a table containing the brightness of the zodiacal light at
     * \f$0.5\,\mu{\rm m}\f$
     * [\f${\rm phot\,s}^{-1}{\rm m}^{-2}\mu{\rm m}^{-1}{\rm arcsec}^{-2}\f$]
     * dependent on the ecliptic coordinates of the object. For the zodiacal
     * light spectrum a slightly reddened solar spectrum is taken.
     *
     * The extinction of zodiacal light by molecular absorption, Rayleigh
     * scattering, and aerosol extinction is considered by input transmission
     * curves provided for airmass = 1. For the scattering curves the
     * effective optical depth is reduced by recipes depending on the
     * top-of-atmosphere line-of-sight zodiacal light intensity. The
     * parametrisations used were obtained by 3D atmospheric scattering
     * calculations for extended emission (see, e.g., Bernstein et al. 2002).
     *
     * \b INPUT:
     * \param spec         desired wavelength grid
     * \param modelpar     sky emission parameters (see typedef of
     *                     ::smparmodel)
     * \param solspec      solar spectrum
     *                     (note: normalisation for Colina et al. 1996)
     * \param molabstrans  transmission curve from radiative transfer code
     *                     (scaled to airmass = 1)
     * \param rscattrans   transmission curve for Rayleigh scattering
     * \param mscattrans   transmission curve for Mie scattering
     *
     * \b OUTPUT:
     * \param spec         spectrum of zodiacal light
     *
     * \b ERRORS:
     * - NDA: No data
     * - IIP: Invalid input parameter(s)
     * - IDG: Inconsistent data grids
     * - see subroutines as well
     */

    smgrid zodtab;
    smspec atrans, rtrans, mtrans;
    char errtxt[SM_MAXLEN+1], zodtabfile[FILENAME_MAX];
    int i = 0;
    double i0 = 0., elong = 0., cs = 0., cl = 0., lgi0 = 0., fr = 0.;
    double cxr = 0., fm = 0., cxm = 0.;

    /* Reference wavelength in mum */
    double lam0 = 0.5;
    /* Scale factor for solar spectrum (Colina et al. 1996)
       -> flux at reference wavelength */
    double scale_solspec = 1921.;
    /* Coefficients for elongation-dependent reddening correction */
    double cs30 = 1.2, cl30 = 0.8, cs90 = 0.9, cl90 = 0.6;
    /* Flux conversion: 10^-8 * W m^-2 mum^-1 sr^-1 * lam [mum]
                     -> phot s^-1 m^-2 mum^-1 arcsec^-2 */
    double conv = 1e-8 * SM_LAM_UNIT /
                  (CPL_PHYS_H * CPL_PHYS_C * SM_SR_IN_ARCSEC2); // = 1.18324
    /* Coefficients for effective tau correction */
    double lgilimr = 2.2444, mr[2] = {1.4072, 0.5267};
    double cr[2] = {-2.6916, -0.7152};
    double lgilimm = 2.2546, mm[2] = {1.3088, 0.4683};
    double cm[2] = {-2.5976, -0.7025};

    /* No error treatment */
    spec->type = 2;

    /* spec = NULL */
    if (spec->n <= 0) {
        sprintf(errtxt, "%s: smspec *spec", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    }

    /* Check validity of input parameters */
    if (modelpar.lon_ecl < -180. || modelpar.lon_ecl > 180. ||
        modelpar.lat_ecl < -90. || modelpar.lat_ecl > 90. ||
        modelpar.alt < 0. || modelpar.alt > 90.) {
        sprintf(errtxt, "%s: lon_ecl, lat_ecl, and/or alt of smparmodel "
                "modelpar out of range", SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Test agreement of wavelength grids in input spectra */
    if ((int) sm_spec_compgrids(solspec, spec) == (int) SM_ERROR_IDG ||
        (int) sm_spec_compgrids(molabstrans, spec) == (int) SM_ERROR_IDG ||
        (int) sm_spec_compgrids(rscattrans, spec) == (int) SM_ERROR_IDG ||
        (int) sm_spec_compgrids(mscattrans, spec) == (int) SM_ERROR_IDG) {
        sprintf(errtxt, "%s: smspec *solspec, *molabstrans, *rscattrans, or "
                "*mscattrans != smspec *spec", SM_ERROR_IDG_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IDG, "%s", errtxt);
    }

    /* Read table of zodiacal light brightness at 0.5 mum */
    sprintf(zodtabfile, "%s/%s", modelpar.datapath, modelpar.zodtabname);
    sm_grid_read(&zodtab, zodtabfile);

    /* Extract intensity for given ecpliptic coordinates */
    sm_grid_extract(&i0, &zodtab, fabs(modelpar.lon_ecl),
                    fabs(modelpar.lat_ecl));

    /* Free memory occupied by zodiacal data table */
    sm_grid_free(&zodtab);

    /* Scale solar spectrum to brightness of zodiacal light */
    sm_spec_calc(spec, '=', solspec);
    sm_spec_scale(spec, '*', i0 / scale_solspec);

    /* Consideration of elongation-dependent reddening of zodiacal light
       (compared to solar spectrum) */

    elong = acos(cos(modelpar.lon_ecl * CPL_MATH_RAD_DEG)
                 * cos(modelpar.lat_ecl * CPL_MATH_RAD_DEG))
            / CPL_MATH_RAD_DEG;

    if (elong <= 30) {
        cs = cs30;
        cl = cl30;
    } else if (elong >= 90) {
        cs = cs90;
        cl = cl90;
    } else {
        cs = ((cs90 - cs30) / 60) * (elong - 30) + cs30;
        cl = ((cl90 - cl30) / 60) * (elong - 30) + cl30;
    }

    for (i = 0; i < spec->n; i++) {
        if (spec->dat[i].lam < lam0) {
            spec->dat[i].flux *= 1.0 + cs * log10(spec->dat[i].lam / lam0);
        } else {
            spec->dat[i].flux *= 1.0 + cl * log10(spec->dat[i].lam / lam0);
        }
    }

    /* Flux in phot s^-1 m^-2 mum^-1 arcsec^-2 */
    sm_spec_convunits(spec, conv, 1);

    /* Correction for molecular absorption of zodiacal light in the lower
       atmosphere */
    sm_spec_copy(&atrans, molabstrans);
    atrans.type = 2;
    sm_spec_scale(&atrans, '^', modelpar.airmass);
    sm_spec_calc(spec, '*', &atrans);
    sm_spec_free(&atrans);

    /* Extinction correction considering reduced loss of scattered light
       for extended emission compared to point source */

    /* Rayleigh scattering (coefficients from 3D scattering calculations) */
    lgi0 = log10(i0);
    if (lgi0 > lgilimr) {
        fr = mr[1] * lgi0 + cr[1];
    } else {
        fr = mr[0] * lgi0 + cr[0];
    }
    cxr = fr * modelpar.airmass;
    sm_spec_copy(&rtrans, rscattrans);
    sm_spec_scale(&rtrans, '^', cxr);
    sm_spec_calc(spec, '*', &rtrans);
    sm_spec_free(&rtrans);

    /* Mie scattering/aerosol extinction (coeff. from 3D scattering calc.) */
    if (lgi0 > lgilimm) {
        fm = mm[1] * lgi0 + cm[1];
    } else {
        fm = mm[0] * lgi0 + cm[0];
    }
    cxm = fm * modelpar.airmass;
    sm_spec_copy(&mtrans, mscattrans);
    sm_spec_scale(&mtrans, '^', cxm);
    sm_spec_calc(spec, '*', &mtrans);
    sm_spec_free(&mtrans);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_comp_telem(smspec *spec, const smparmodel modelpar)
{
    /*!
     * \callgraph
     * Computes thermal emission by telesope/instrument.
     * Several components are possible. The routine assumes grey bodies
     * depending on emissivity and temperature for each component. The order
     * of the components matters. The first emissivity and temperature have to
     * be for the main mirror. Any emission or absorption which is not related
     * to the listed optical components is not considered.
     *
     * The routine assumes that the absorption inside the instrument is
     * handled by a response curve, which is applied somewhere else. For this
     * reason, the emission of the different components is raised depending
     * on the position along the light path. For example, the emission of the
     * main mirror is increased by the reciprocal of the mirror transmission.
     * This factor makes sure that the telescope emission can be added to the
     * other sky model components, for which the flux is reduced by absorption
     * of the main mirror and all other optical components. The routine
     * simulates the apparent telescope/instrument emission that is not
     * removed in the course of typical astronomical flux calibration.
     *
     * \b INPUT:
     * \param spec       desired wavelength grid
     * \param modelpar   sky emission parameters (see typedef of ::smparmodel)
     *
     * \b OUTPUT:
     * \param spec       apparent telescope/instrument emission
     *
     * \b ERRORS:
     * - NDA: No data
     * - IIP: Invalid input parameter(s)
     * - IDR: Invalid data range
     */

    char errtxt[SM_MAXLEN+1];
    int j = 0, i = 0;
    double c1 = 0., c2 = 0., c3 = 0., c4 = 0., c5 = 0., corr = 0.;

    /* Constants and unit conversions */

    c1 = 2 * CPL_PHYS_H * pow(CPL_PHYS_C, 2);
    c2 = CPL_PHYS_H * CPL_PHYS_C / CPL_PHYS_K;
    c3 = SM_LAM_UNIT / (CPL_PHYS_H * CPL_PHYS_C * SM_SR_IN_ARCSEC2);
    c4 = c1 * c3 / pow(SM_LAM_UNIT, 4);
    c5 = c2 / SM_LAM_UNIT;

    /* No error treatment */
    spec->type = 2;

    /* spec = NULL */
    if (spec->n <= 0) {
        sprintf(errtxt, "%s: smspec *spec", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    }

    /* Flux = 0 */
    sm_spec_scale(spec, '=', 0);

    /* Check number of emissivity/temperature pairs */

    if (modelpar.ncomp <= 0) {
        sprintf(errtxt, "%s: ncomp of smparmodel modelpar <= 0",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Check limits of input parameters */

    for (j = 0; j < modelpar.ncomp; j++) {
        if (modelpar.eps[j] < 0 || modelpar.eps[j] > 1 ||
            modelpar.temp[j] < 0) {
            sprintf(errtxt, "%s: components of eps and/or temp of "
                    "smparmodel modelpar out of range", SM_ERROR_IIP_TXT);
            return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s",
                                         errtxt);
        }
    }

    /* Compute and add grey body emission of each component */

    for (corr = 1., j = 0; j < modelpar.ncomp; j++) {

        /* Compute correction factor for reduced absorption losses */
        corr *= 1. / (1. - modelpar.eps[j]);

        for (i = 0; i < spec->n; i++) {
            if (spec->dat[i].lam <= 0) {
                sprintf(errtxt, "%s: smspec *spec (wavelengths <= 0)",
                        SM_ERROR_IDR_TXT);
                return cpl_error_set_message(cpl_func, SM_ERROR_IDR, "%s",
                                             errtxt);
            }
            if (modelpar.temp[j] > 0) {
                spec->dat[i].flux += corr * modelpar.eps[j] * c4 /
                                     (pow(spec->dat[i].lam, 4) * (exp(c5 /
                                     (spec->dat[i].lam * modelpar.temp[j]))
                                      - 1));
            }
        }

    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_comp_extrapolrad(smspec *spec, const smspec *trans,
                                   const smparmodel modelpar,
                                   const double lim[2])
{
    /*!
     * \callgraph
     * Extracts and rebins library LBLRTM/RFM molecular radiance spectrum
     * and extrapolates it to the requested airmass, which is computed by
     * means of the formula of Rozenberg (1966). For the optical depth
     * dependent extrapolation, the corresponding transmission curve as
     * provided by ::sm_comp_extrapoltrans is used.
     *
     * \b INPUT:
     * \param spec       desired wavelength grid
     * \param trans      transmission curve belonging to radiance spectrum
     *                   (scaled to airmass = 1)
     * \param modelpar   sky emission parameters (see typedef of ::smparmodel)
     * \param lim        minimum and maximum wavelength extraction limits
     *
     * \b OUTPUT:
     * \param spec       extracted, rebinned, and scaled LBLRTM/RFM radiance
     *                   spectrum
     *
     * \b ERRORS:
     * - NDA: No data
     * - IOD: Invalid order of data points
     * - IIP: Invalid input parameter(s)
     * - IDG: Inconsistent data grids
     */

    cpl_error_code err = CPL_ERROR_NONE;
    cpl_parameter *p;
    cpl_parameterlist *libfilepar;
    smspec tmptrans;
    char errtxt[SM_MAXLEN+1];
    double xz = 1.;

    /* spec = NULL */
    if (spec->n <= 0) {
        sprintf(errtxt, "%s: smspec *spec", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    }

    /* trans = NULL */
    if (spec->n <= 0) {
        sprintf(errtxt, "%s: smspec *trans", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    }

    /* Check input wavelength limits */
    if (lim[1] <= lim[0]) {
        sprintf(errtxt, "%s: lim[1] <= lim[0] (input wavelength limits)",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Valid target altitude? */
    if (modelpar.alt < 0 || modelpar.alt > 90) {
        sprintf(errtxt, "%s: alt of smparmodel modelpar < 0 or > 90 degrees",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Search in the library of molecular spectra for a file with the desired
       parameters */
    libfilepar = cpl_parameterlist_new();
    if ((err = sm_comp_readlibstruct(libfilepar, modelpar, "R")) !=
        CPL_ERROR_NONE) {
        cpl_parameterlist_delete(libfilepar);
        return err;
    }

    /* Get radiance spectrum from library */
    if ((err = sm_comp_getmolspec(spec, libfilepar, lim)) !=
        CPL_ERROR_NONE) {
        cpl_parameterlist_delete(libfilepar);
        return err;
    }

    /* Compare wavelength grids of radiance and transmission spectrum */
    if ((int) sm_spec_compgrids(trans, spec) == (int) SM_ERROR_IDG) {
        sprintf(errtxt, "%s: smspec *trans != smspec *spec",
                SM_ERROR_IDG_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IDG, "%s", errtxt);
    }

    /* Get airmass of selected library file */
    if ((p = cpl_parameterlist_find(libfilepar, "airmass")) != NULL) {
        xz = cpl_parameter_get_double(p);
    } else if ((p = cpl_parameterlist_find(libfilepar, "alt")) != NULL) {
        /* Airmass according to Rozenberg (1966) */
        xz = sm_comp_alttoairmass(cpl_parameter_get_double(p));
    } else {
        /* Assume airmass = 1 */
        xz = 1;
    }

    /* Estimate radiance spectrum for requested airmass by means of
       transmission curve (scaled to airmass = 1) */
    sm_spec_copy(&tmptrans, trans);
    sm_spec_scale(&tmptrans, '^', xz);
    sm_spec_scale(&tmptrans, '*', fabs(modelpar.airmass - xz) / xz);
    if (modelpar.airmass < xz) {
        sm_spec_scale(&tmptrans, '_', 1.);
    } else {
        sm_spec_scale(&tmptrans, '+', 1.);
    }
    sm_spec_calc(spec, '*', &tmptrans);

    /* Free memory */
    cpl_parameterlist_delete(libfilepar);
    sm_spec_free(&tmptrans);

    return err;
}


cpl_error_code sm_comp_getlinespec(smspec *spec, const smparmodel modelpar,
                                   const smspec *rscattrans,
                                   const smspec *mscattrans)
{
    /*!
     * \callgraph
     * Converts table of airglow emission lines (see Hanuschik 2003 and
     * Rousselot et al. 2000) to spectrum (::smspec structure) by applying
     * wavelength-dependent correction factors which depend on object
     * altitude, emission layer height, solar radio flux, period of the year,
     * and period of the night (see ::sm_comp_readvarpar). The wavelengths are
     * converted from vacuum to air if required. Moreover, the change of line
     * intensity by molecular absorption, Rayleigh scattering, and aerosol
     * extinction is considered. The former correction is realised by a
     * precalculated library of lists of line-specific transmission correction
     * factors for each radiative transfer code transmission spectrum. The
     * latter corrections are performed by provided extinction curves for
     * which the effective optical depth is reduced by airmass-dependent
     * recipes which were obtained by 3D atmospheric scattering calculations
     * for extended airglow emission.
     *
     * \b INPUT:
     * \param spec        desired wavelength grid
     * \param modelpar    airglow parameters (see typedef of ::smparmodel)
     * \param rscattrans  transmission curve for Rayleigh scattering
     * \param mscattrans  transmission curve for Mie scattering
     *
     * \b OUTPUT:
     * \param spec        airglow emission line spectrum
     *
     * \b ERRORS:
     * - IIP: Invalid input parameter(s)
     * - IDG: Inconsistent data grids
     * - see subroutines as well
     */

    FILE *stream;
    cpl_error_code err = CPL_ERROR_NONE;
    cpl_parameterlist *libfilepar;
    cpl_parameter *p;
    cpl_table *varpar, *linedat;
    smspec linetrans, linetab, rtrans, mtrans;
    char libfile[FILENAME_MAX], linetabfile[FILENAME_MAX];
    char errtxt[SM_MAXLEN+1];
    int i = 0;
    double xz = 0., z = 0., xag = 0., lgx = 0., fr = 0., cxr = 0., fm = 0.;
    double cxm = 0.;

    /* Coefficients for effective tau correction */
    double mr = 1.6693, cr = -0.1460, mm = 1.7320, cm = -0.3175;

    /* Search in the library of molecular spectra for a file with the desired
       parameters */
    libfilepar = cpl_parameterlist_new();
    if ((err = sm_comp_readlibstruct(libfilepar, modelpar, "L")) !=
        CPL_ERROR_NONE) {
        cpl_parameterlist_delete(libfilepar);
        return err;
    }

    /* Get file name from parameter list */
    p = cpl_parameterlist_find(libfilepar, "filename");
    strcpy(libfile, cpl_parameter_get_string(p));

    /* Check file existence */
    if ((stream = fopen(libfile, "r")) == NULL) {
        cpl_parameterlist_delete(libfilepar);
        sprintf(errtxt, "%s: %s", SM_ERROR_FOF_TXT, libfile);
        return cpl_error_set_message(cpl_func, SM_ERROR_FOF, "%s", errtxt);
    }
    fclose(stream);

    /* Check file name extension and read line list */
    if (strstr(libfile, ".dat") != NULL ||
        strstr(libfile, ".ascii") != NULL) {
        sm_spec_read(&linetrans, libfile);
    } else if (strstr(libfile, ".fits") != NULL ||
               strstr(libfile, ".mt") != NULL) {
        sm_spec_readfits(&linetrans, libfile);
    } else {
        cpl_parameterlist_delete(libfilepar);
        sprintf(errtxt, "%s: %s ('dat', 'ascii', 'fits', or 'mt' only)",
                SM_ERROR_IFE_TXT, libfile);
        return cpl_error_set_message(cpl_func, SM_ERROR_IFE, "%s", errtxt);
    }

    /* Get airmass of selected library file */
    if ((p = cpl_parameterlist_find(libfilepar, "airmass")) != NULL) {
        xz = cpl_parameter_get_double(p);
    } else if ((p = cpl_parameterlist_find(libfilepar, "alt")) != NULL) {
        /* Airmass according to Rozenberg (1966) */
        xz = sm_comp_alttoairmass(cpl_parameter_get_double(p));
    } else {
        /* Assume airmass = 1 */
        xz = 1;
    }

    /* Free memory */
    cpl_parameterlist_delete(libfilepar);

    /* Read table of airglow emission lines */
    sprintf(linetabfile, "%s/%s", modelpar.datapath, modelpar.linetabname);
    sm_spec_read(&linetab, linetabfile);

    /* Test agreement of line wavelengths */
    if ((int) sm_spec_compgrids(&linetab, &linetrans) == (int) SM_ERROR_IDG) {
        sm_spec_free(&linetrans);
        sm_spec_free(&linetab);
        sprintf(errtxt, "%s: smspec *linetab != smspec *linetrans",
                SM_ERROR_IDG_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IDG, "%s", errtxt);
    }

    /* Check validity of object altitude and, hence, airmass */
    if (modelpar.alt < 0. || modelpar.alt > 90.) {
        sm_spec_free(&linetrans);
        sm_spec_free(&linetab);
        sprintf(errtxt, "%s: alt of smparmodel modelpar < 0 or > 90 degrees",
                SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Scale line transmission values to given airmass */
    sm_spec_scale(&linetrans, '^', modelpar.airmass / xz);

    /* Correct airglow line strengths for molecular absorption in the lower
       atmosphere */
    for (i = 0; i < linetab.n; i++) {
        linetab.dat[i].flux *= linetrans.dat[i].flux;
    }

    /* Convert vacuum to air wavelengths if required */
    if (strncmp(modelpar.vac_air, "air", 3) == 0) {
        sm_comp_vactoair_spec(&linetab);
    }

    /* Read airglow scaling parameters from file */
    varpar = cpl_table_new(0);
    sm_comp_readvarpar(varpar, modelpar);

    /* Scale line fluxes depending on correction factors provided by varpar */
    linedat = cpl_table_new(0);
    sm_comp_scalelinetab(linedat, &linetab, varpar);

    /* Convert line table to spectrum */
    sm_comp_convertlinetab(spec, linedat);

    /* Free memory */
    sm_spec_free(&linetrans);
    sm_spec_free(&linetab);
    cpl_table_delete(varpar);
    cpl_table_delete(linedat);

    /* Test agreement of wavelength grids in input spectra */
    if ((int) sm_spec_compgrids(rscattrans, spec) == (int) SM_ERROR_IDG) {
        sprintf(errtxt, "%s: smspec *rscattrans != smspec *spec",
                SM_ERROR_IDG_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IDG, "%s", errtxt);
    }
    if ((int) sm_spec_compgrids(mscattrans, spec) == (int) SM_ERROR_IDG) {
        sprintf(errtxt, "%s: smspec *mscattrans != smspec *spec",
                SM_ERROR_IDG_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IDG, "%s", errtxt);
    }

    /* Extinction correction considering reduced loss of scattered light
       for extended emission compared to point source */

    /* Calculate airglow airmass (based on van Rhijn formula for 90 km) */
    z = (90. - modelpar.alt) * CPL_MATH_RAD_DEG;
    xag = 1 / sqrt(1.0 - 0.972 * pow(sin(z), 2));

    /* Rayleigh scattering (coefficients from 3D scattering calculations) */
    lgx = log10(xag);
    fr = mr * lgx + cr;
    cxr = fr * xag;
    sm_spec_copy(&rtrans, rscattrans);
    sm_spec_scale(&rtrans, '^', cxr);
    sm_spec_calc(spec, '*', &rtrans);
    sm_spec_free(&rtrans);

    /* Mie scattering/aerosol extinction (coeff. from 3D scattering calc.) */
    fm = mm * lgx + cm;
    cxm = fm * xag;
    sm_spec_copy(&mtrans, mscattrans);
    sm_spec_scale(&mtrans, '^', cxm);
    sm_spec_calc(spec, '*', &mtrans);
    sm_spec_free(&mtrans);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_comp_readvarpar(cpl_table *varpar,
                                  const smparmodel modelpar)
{
    /*!
     * Reads data related to airglow scaling from file and fills a CPL table
     * with feature-related correction factors and uncertainties.
     * The resulting factors are related to the standard strengths of
     * airglow emission features of the upper atmosphere given by Hanuschik
     * (2003) and Rousselot et al. (2000).
     * The scaling factors depend on the emission layer width (related
     * to airmass), the monthly-averaged solar flux in sfu, the season of
     * the year, and the time of the day.
     * In addition, the relative Doppler widths \f$\sigma_{\rm D}\f$ of the
     * different features are provided by the output table. They are derived
     * from the read mol masses and temperatures of the emission layers.
     *
     * \b INPUT:
     * \param modelpar  sky emission parameters (see typedef of ::smparmodel)
     *
     * \b OUTPUT:
     * \param varpar    airglow variability data as CPL table with columns
     *   - \e N:        feature number
     *   - \e fac:      scaling factor for emission lines depending on feature
     *   - \e dfac:     uncertainties of factors
     *   - \e relwidth: relative Doppler width
     *
     * \b ERRORS:
     * - FOF: File opening failed
     * - UFS: Unexpected file structure
     * - IIP: Invalid input parameter(s)
     */

    FILE *stream;
    smparam x[SM_MAXPAR], m[SM_MAXPAR];
    char vardatfile[FILENAME_MAX], errtxt[SM_MAXLEN+1];
    int nfeat = 0, nseason = 0, ntime = 0;
    int i = 0, j = 0;
    double molmass = 0., temp = 0., height = 0., scale = 0., cons = 0.;
    double slope = 0., relwid = 0., mean = 0., sig = 0., z = 0., cvr = 0.;
    double csol = 0., fac = 0.;
    double cdopp = 9.618e-9; // constant for calculation of Doppler width

    /* Check existence of airglow scaling parameter file */
    sprintf(vardatfile, "%s/%s", modelpar.datapath, modelpar.vardatname);
    if ((stream = fopen(vardatfile, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", SM_ERROR_FOF_TXT, vardatfile);
        return cpl_error_set_message(cpl_func, SM_ERROR_FOF, "%s", errtxt);
    }

    /* Read airglow scaling parameter file */

    /* Read numbers for data set size (N_bin = (nseason + 1) * (ntime + 1)) */
    sm_param_read(stream, x);
    nfeat = x[0].i;
    sm_param_read(stream, x);
    nseason = x[0].i;
    ntime = x[1].i;

    /* Return NULL vectors in the case of one or more zero for the data set
       size */
    if (nfeat == 0 || nseason == 0 || ntime == 0) {
        fclose(stream);
        sprintf(errtxt, "%s: %s (nfeat == 0 || nseason == 0 || ntime == 0)",
                SM_ERROR_UFS_TXT, vardatfile);
        return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s", errtxt);
    }

    /* Check modelpar input parameters */
    if (modelpar.alt < 0. || modelpar.alt > 90. || modelpar.msolflux < 0. ||
        modelpar.season < 0 || modelpar.season > nseason ||
        modelpar.time < 0 || modelpar.time > ntime) {
        fclose(stream);
        sprintf(errtxt, "%s: alt, msolflux, season, and/or time of "
                "smparmodel modelpar out of range", SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Prepare CPL table for feature scaling data */
    cpl_table_set_size(varpar, nfeat);
    cpl_table_new_column(varpar, "N", CPL_TYPE_INT);
    cpl_table_new_column(varpar, "fac", CPL_TYPE_DOUBLE);
    cpl_table_new_column(varpar, "dfac", CPL_TYPE_DOUBLE);
    cpl_table_new_column(varpar, "relwidth", CPL_TYPE_DOUBLE);

    /* Read data for each feature, extract bin data for selected time, and
       fill smspec structure with these variability data */

    for (i = 0; i < nfeat; i++) {

        /* Time-independent parameters */
        sm_param_read(stream, x);
        molmass = x[0].d;
        sm_param_read(stream, x);
        temp = x[0].d;
        sm_param_read(stream, x);
        height = x[0].d;
        sm_param_read(stream, x);
        scale = x[0].d;
        sm_param_read(stream, x);
        cons = x[0].d;
        slope = x[1].d;

        /* Calculate relative Doppler width (sigma_D) */
        relwid = cdopp * sqrt(temp / molmass);

        /* Set feature number and Doppler width in output table */
        cpl_table_set(varpar, "N", i, i + 1);
        cpl_table_set(varpar, "relwidth", i, relwid);

        /* Mean value for selected time */
        for (j = 0; j < ntime + 1; j++) {
            sm_param_read(stream, m);
            if (j == modelpar.time) {
                mean = m[modelpar.season].d;
            }
        }

        /* Standard deviation for selected time */
        for (j = 0; j < ntime + 1; j++) {
            sm_param_read(stream, m);
            if (j == modelpar.time) {
                sig = m[modelpar.season].d;
            }
        }

        /* Emission line brightness depending on airmass and layer height
           (van Rhijn 1921) */
        z = (90. - modelpar.alt) * CPL_MATH_RAD_DEG;
        cvr = 1 / sqrt(1 - pow((SM_ERAD / (SM_ERAD + height)) * sin(z), 2));

        /* Influence of solar radio flux [sfu] */
        csol = slope * modelpar.msolflux + cons;

        /* Set feature-specific correction factors and their uncertainties */
        fac = scale * cvr * csol * mean;
        cpl_table_set(varpar, "fac", i, fac);
        cpl_table_set(varpar, "dfac", i, fac * sig / mean);

    }

    /* Close file */
    fclose(stream);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_comp_scalelinetab(cpl_table *linedat, const smspec *linetab,
                                    const cpl_table *varpar)
{
    /*!
     * Writes line data in CPL table and corrects line fluxes by factors
     * provided by \e varpar for the different features.
     *
     * \b INPUT:
     * \param linetab  table of emission lines (provided as ::smspec
     *                 structure; dflux1 contains feature number)
     * \param varpar   line parameters (scaling factors and Doppler widths)
     *                 provided as CPL table
     *
     * \b OUTPUT:
     * \param linedat  CPL table with columns "lam", "flux", "dflux", "feat",
     *                 and "width"
     *
     * \b ERRORS:
     * - NDA: No data
     * - ISD: Insufficient data points
     * - IOV: Invalid object value(s)
     */

    smbool fl_feat = T;
    char errtxt[SM_MAXLEN+1];
    int nfeat = 0, i = 0, feat = 0, j = 0;
    double lam = 0., flux = 0., dflux = 0., width = 0.;

    /* Initialise output CPL table */
    cpl_table_set_size(linedat, linetab->n);
    cpl_table_new_column(linedat, "lam", CPL_TYPE_DOUBLE);
    cpl_table_new_column(linedat, "flux", CPL_TYPE_DOUBLE);
    cpl_table_new_column(linedat, "dflux", CPL_TYPE_DOUBLE);
    cpl_table_new_column(linedat, "width", CPL_TYPE_DOUBLE);
    cpl_table_new_column(linedat, "feat", CPL_TYPE_INT);

    /* Check number of features in parameter list */
    nfeat = cpl_table_get_nrow(varpar);
    if (nfeat == 0) {
        sprintf(errtxt, "%s: cpl_table *varpar", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    }

    for (i = 0; i < linetab->n; i++) {

        /* Get feature number */
        feat = (int) linetab->dat[i].dflux1;
        if (feat < 1 || feat > nfeat) {
            /* Take first feature in the case of invalid input */
            fl_feat = F;
            feat = 1;
        }
        j = feat - 1;

        /* Get correction factors and line widths for given feature */
        lam = linetab->dat[i].lam;
        flux = cpl_table_get(varpar, "fac", j, NULL) * linetab->dat[i].flux;
        dflux = cpl_table_get(varpar, "dfac", j, NULL) * linetab->dat[i].flux;
        width = cpl_table_get(varpar, "relwidth", j, NULL) * lam;

        /* Write line data into CPL table */
        cpl_table_set(linedat, "lam", i, lam);
        cpl_table_set(linedat, "flux", i, flux);
        cpl_table_set(linedat, "dflux", i, dflux);
        cpl_table_set(linedat, "width", i, width);
        cpl_table_set(linedat, "feat", i, feat);

    }

    if (fl_feat == F) {
        sprintf(errtxt, "%s: invalid feature number(s) in smspec *linetab",
                SM_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IOV, "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_comp_convertlinetab(smspec *spec, const cpl_table *linedat)
{
    /*!
     * Converts table of emission lines with central wavelength, line flux,
     * flux uncertainty, and line width into line spectrum with given
     * wavelength grid. Gaussian line shape is assumed.
     *
     * \b INPUT:
     * \param spec     desired wavelength grid as ::smspec structure
     * \param linedat  CPL table with line data
     *
     * \b OUTPUT:
     * \param spec     sky line spectrum
     *
     * \b ERRORS:
     * - NDA: No data
     * - ISD: Insufficient data points
     */

    cpl_table *gauss;
    cpl_array *pixlim, *dpix;
    char errtxt[SM_MAXLEN+1];
    int i = 0, nsigbin = 0, ngauss = 0, h = 0, nlin = 0, imin = 0, j = 0;
    double *llim = NULL, *dlam = NULL, *sig = NULL, *psig = NULL;
    const double *lam = NULL, *flux = NULL, *dflux = NULL, *width = NULL;
    double minlam = 0., maxlam = 0., dl = 0., dlmin = 0., maxwidth = 0.;
    double sigfrac = 0., x = 0., pxsum = 0., lamfac = 1e-5;
    double lmax = 0., lmin = 0., psum = 0., l0 = 0.;

    /* Allow errors in output spectrum */
    spec->type = 3;
    sm_spec_scale(spec, '=', 0.);

    /* Check number of data points in spectrum */
    if (spec->n <= 0) {
        sprintf(errtxt, "%s: smspec *spec", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    } else if (spec->n == 1) {
        sprintf(errtxt, "%s: smspec *spec (spec->n == 1)", SM_ERROR_ISD_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_ISD, "%s", errtxt);
    }

    /* Create CPL array for wavelength limits of pixels in output spectrum */
    pixlim = cpl_array_new((spec->n)+1, CPL_TYPE_DOUBLE);
    cpl_array_set(pixlim, 0, 1.5 * spec->dat[0].lam - 0.5 * spec->dat[1].lam);
    for (i = 1; i < spec->n; i++) {
        cpl_array_set(pixlim, i,
                      0.5 * (spec->dat[i-1].lam + spec->dat[i].lam));
    }
    cpl_array_set(pixlim, spec->n, 1.5 * spec->dat[(spec->n)-1].lam
                  - 0.5 * spec->dat[(spec->n)-2].lam);
    llim = cpl_array_get_data_double(pixlim);

    /* Get wavelength range for line selection */
    minlam = llim[0] * (1 - lamfac);
    maxlam = llim[spec->n] * (1 + lamfac);

    /* Create CPL array for wavelength ranges of pixels in output spectrum */
    dpix = cpl_array_new(spec->n, CPL_TYPE_DOUBLE);
    for (dlmin = HUGE_VAL, i = 0; i < spec->n; i++) {
        dl = llim[i+1] - llim[i];
        cpl_array_set(dpix, i, dl);
        if (dl < dlmin) {
            dlmin = dl;
        }
    }
    dlam = cpl_array_get_data_double(dpix);

    /* Get suitable resolution for Gaussian */
    maxwidth = cpl_table_get_column_max(linedat, "width");
    nsigbin = (int) ceil(SM_NSIGBIN * maxwidth / dlmin);
    if (nsigbin < SM_NSIGBIN) {
        nsigbin = SM_NSIGBIN;
    }
    sigfrac = 1. / ((double) nsigbin);

    /* Create CPL table for Gaussian data */
    ngauss = 2 * (int) (SM_SIGMAX / sigfrac);
    gauss = cpl_table_new(ngauss);
    cpl_table_new_column(gauss, "sig", CPL_TYPE_DOUBLE);
    cpl_table_new_column(gauss, "psig", CPL_TYPE_DOUBLE);

    /* Get pointers to CPL table columns */
    sig = cpl_table_get_data_double(gauss, "sig");
    psig = cpl_table_get_data_double(gauss, "psig");

    /* Calculate Gaussian */
    x = - SM_SIGMAX + 0.5 * sigfrac;
    for (pxsum = 0, h = 0; h < ngauss; h++) {
        sig[h] = x;
        psig[h] = exp(-0.5 * pow(x, 2));
        pxsum += psig[h];
        x += sigfrac;
    }
    for (h = 0; h < ngauss; h++) {
        psig[h] /= pxsum;
    }

    /* Get number of lines */
    nlin = cpl_table_get_nrow(linedat);

    /* Get pointers to CPL table columns */
    lam = cpl_table_get_data_double_const(linedat, "lam");
    flux = cpl_table_get_data_double_const(linedat, "flux");
    dflux = cpl_table_get_data_double_const(linedat, "dflux");
    width = cpl_table_get_data_double_const(linedat, "width");

    /* Add lines to output spectrum */

    for (imin = 0, j = 0; j < nlin; j++) {

        /* Skip lines if outside required range */
        if (lam[j] < minlam || lam[j] > maxlam) {
            continue;
        }

        /* Skip lines shortwards of lower wavelength limit */
        lmax = lam[j] + SM_SIGMAX * width[j];
        if (lmax < llim[0]) {
            continue;
        }

        /* Find lowest spectrum pixel relevant for given line */
        lmin = lam[j] - SM_SIGMAX * width[j];
        while (imin < spec->n && lmin >= llim[imin+1]) {
            imin++;
        }

        /* Skip lines longwards of upper wavelength limit */
        if (imin == spec->n) {
            break;
        }

        /* Treat lines that cover only one pixel in output spectrum */
        if (lmax < llim[imin+1]) {
            spec->dat[imin].flux += flux[j] / dlam[imin];
            spec->dat[imin].dflux1 += dflux[j] / dlam[imin];
            continue;
        }

        /* Treat lines that cover several pixels in output spectrum */
        for (i = imin, psum = 0., h = 0; h < ngauss; h++) {
            l0 = lam[j] + sig[h] * width[j];
            if (l0 >= llim[0] && l0 < llim[i+1]) {
                psum += psig[h];
            }
            if (l0 >= llim[i+1] || h == ngauss - 1) {
                /* psum = ratio of pixel-related flux to total line flux */
                spec->dat[i].flux += psum * flux[j] / dlam[i];
                spec->dat[i].dflux1 += psum * dflux[j] / dlam[i];
                i++;
                if (i == spec->n) {
                    break;
                }
                psum = psig[h];
            }
        }

    }

    /* Copy lower error to upper error */
    sm_spec_changetype(spec, 4);
    sm_spec_changetype(spec, 3);

    /* Free memory */
    cpl_array_delete(pixlim);
    cpl_array_delete(dpix);
    cpl_table_delete(gauss);

    return CPL_ERROR_NONE;
}


cpl_error_code sm_comp_convertlinetabo(smspec *spec, const cpl_table *linedat)
{
    /*!
     * Converts table of emission lines with central wavelength, line flux,
     * flux uncertainty, and line width into line spectrum with given
     * wavelength grid.
     *
     * \b INPUT:
     * \param spec     desired wavelength grid as ::smspec structure
     * \param linedat  CPL table with line data
     *
     * \b OUTPUT:
     * \param spec     sky line spectrum
     *
     * \b ERRORS:
     * - NDA: No data
     * - ISD: Insufficient data points
     */

    char errtxt[SM_MAXLEN+1];
    int nlin = 0, i = 0, j = 0;
    const double *lam = NULL, *flux = NULL, *dflux = NULL;
    double lmin = 0, lmax = 0, dl = 0;

    /* Allow errors in output spectrum */
    spec->type = 3;

    /* Check number of data points in spectrum */
    if (spec->n <= 0) {
        sprintf(errtxt, "%s: smspec *spec", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    } else if (spec->n == 1) {
        sprintf(errtxt, "%s: smspec *spec (spec->n == 1)", SM_ERROR_ISD_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_ISD, "%s", errtxt);
    }

    /* Get number of lines */
    nlin = cpl_table_get_nrow(linedat);

    /* Get pointers to CPL table columns */
    lam = cpl_table_get_data_double_const(linedat, "lam");
    flux = cpl_table_get_data_double_const(linedat, "flux");
    dflux = cpl_table_get_data_double_const(linedat, "dflux");

    for (i = 0; i < spec->n; i++) {

        /* Limits of wavelength bin in spectrum */

        if (i == 0) {
            lmin = 1.5 * spec->dat[i].lam - 0.5 * spec->dat[i+1].lam;
        } else {
            lmin = lmax;
        }

        if (i == spec->n - 1) {
            lmax = 1.5 * spec->dat[i].lam - 0.5 * spec->dat[i-1].lam;
        } else {
            lmax = 0.5 * (spec->dat[i].lam + spec->dat[i+1].lam);
        }

        dl = lmax - lmin;

        /* Add line fluxes in wavelength bin */

        spec->dat[i].flux = 0.;
        spec->dat[i].dflux1 = 0.;
        spec->dat[i].dflux2 = 0.;

        while (i == 0 && j < nlin && lam[j] <= lmin) {
            j++;
        }

        while (j < nlin && lam[j] <= lmax) {
            spec->dat[i].flux += flux[j];
            spec->dat[i].dflux1 += dflux[j];
            j++;
        }

        /* Divide line flux sum by bin width */
        spec->dat[i].flux /= dl;
        spec->dat[i].dflux1 /= dl;

        /* Update upper error for possible future use */
        spec->dat[i].dflux2 = spec->dat[i].dflux1;

    }

    return CPL_ERROR_NONE;
}


cpl_error_code sm_comp_airglowcont(smspec *spec, const smparmodel modelpar,
                                   const smspec *molabstrans,
                                   const smspec *rscattrans,
                                   const smspec *mscattrans)
{
    /*!
     * \callgraph
     * Reads airglow continuum spectrum (with uncertainties) and scales it
     * depending on the emission layer width (related to airmass), the
     * monthly-averaged solar flux in sfu, the season of the year, and the
     * time of the day. Moreover, the change of continuum emission by
     * molecular absorption, Rayleigh scattering, and aerosol extinction is
     * considered by provided transmission curves. For the scattering curves
     * the effective optical depth is reduced by airmass-dependent recipes
     * which were obtained by 3D atmospheric scattering calculations for
     * extended airglow emission.
     *
     * \b INPUT:
     * \param spec         desired wavelength grid
     * \param modelpar     sky emission parameters (see typedef of
     *                     ::smparmodel)
     * \param molabstrans  transmission curve from radiative transfer code
     *                     (scaled to airmass = 1)
     * \param rscattrans   transmission curve for Rayleigh scattering
     * \param mscattrans   transmission curve for Mie scattering
     *
     * \b OUTPUT:
     * \param spec         airglow continuum
     *
     * \b ERRORS:
     * - NDA: No data
     * - IDG: Inconsistent data grids
     * - FOF: File opening failed
     * - UFS: Unexpected file structure
     * - IIP: Invalid input parameter(s)
     * - ISM: Insufficient memory
     */

    FILE *stream;
    smspec acont, atrans, rtrans, mtrans;
    smparam x[SM_MAXPAR], m[SM_MAXPAR];
    char errtxt[SM_MAXLEN+1], acontfile[FILENAME_MAX];
    int nseason = 0, ntime = 0, ndat = 0, j = 0, i = 0;
    double height = 0., scale = 0., cons = 0., slope = 0., mean = 0.;
    double sig = 0., z = 0., cvr = 0., csol = 0., corr = 0., rsig = 0.;
    double rx = 0., xag = 0., lgx = 0., fr = 0., cxr = 0., fm = 0., cxm = 0.;

    /* Coefficients for effective tau correction */
    double mr = 1.6693, cr = -0.1460, mm = 1.7320, cm = -0.3175;

    /* spec = NULL */
    if (spec->n <= 0) {
        sprintf(errtxt, "%s: smspec *spec", SM_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_NDA, "%s", errtxt);
    }

    /* Test agreement of wavelength grids in input spectra */
    if ((int) sm_spec_compgrids(molabstrans, spec) == (int) SM_ERROR_IDG) {
        sprintf(errtxt, "%s: smspec *molabstrans != smspec *spec",
                SM_ERROR_IDG_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IDG, "%s", errtxt);
    }
    if ((int) sm_spec_compgrids(rscattrans, spec) == (int) SM_ERROR_IDG) {
        sprintf(errtxt, "%s: smspec *rscattrans != smspec *spec",
                SM_ERROR_IDG_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IDG, "%s", errtxt);
    }
    if ((int) sm_spec_compgrids(mscattrans, spec) == (int) SM_ERROR_IDG) {
        sprintf(errtxt, "%s: smspec *mscattrans != smspec *spec",
                SM_ERROR_IDG_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IDG, "%s", errtxt);
    }

    /* Check existence of airglow continuum data file */
    sprintf(acontfile, "%s/%s", modelpar.datapath, modelpar.acontname);
    if ((stream = fopen(acontfile, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", SM_ERROR_FOF_TXT, acontfile);
        return cpl_error_set_message(cpl_func, SM_ERROR_FOF, "%s", errtxt);
    }

    /* Read airglow continuum data file */

    /* Read numbers for data set size (N_bin = (nseason + 1) * (ntime + 1)) */
    sm_param_read(stream, x);
    nseason = x[0].i;
    ntime = x[1].i;
    sm_param_read(stream, x);
    ndat = x[0].i;

    /* Check data set size */
    if (nseason == 0 || ntime == 0 || ndat == 0) {
        fclose(stream);
        sprintf(errtxt, "%s: %s (nseason == 0 || ntime == 0 || ndat == 0)",
                SM_ERROR_UFS_TXT, acontfile);
        return cpl_error_set_message(cpl_func, SM_ERROR_UFS, "%s", errtxt);
    }

    /* Check modelpar input parameters */
    if (modelpar.alt < 0. || modelpar.alt > 90. || modelpar.msolflux < 0. ||
        modelpar.season < 0 || modelpar.season > nseason ||
        modelpar.time < 0 || modelpar.time > ntime) {
        fclose(stream);
        sprintf(errtxt, "%s: alt, msolflux, season, and/or time of "
                "smparmodel modelpar out of range", SM_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_IIP, "%s", errtxt);
    }

    /* Allocate memory for airglow continuum data */
    sm_spec_malloc(&acont, ndat);
    if (acont.n == 0) {
        fclose(stream);
        sprintf(errtxt, "%s: smspec acont", SM_ERROR_ISM_TXT);
        return cpl_error_set_message(cpl_func, SM_ERROR_ISM, "%s", errtxt);
    }

    /* Read variability data and extract bin data for selected time */

    /* Time-independent parameters */
    sm_param_read(stream, x);
    height = x[0].d;
    sm_param_read(stream, x);
    scale = x[0].d;
    sm_param_read(stream, x);
    cons = x[0].d;
    slope = x[1].d;

    /* Mean value for selected time */
    for (j = 0; j < ntime + 1; j++) {
        sm_param_read(stream, m);
        if (j == modelpar.time) {
            mean = m[modelpar.season].d;
        }
    }

    /* Standard deviation for selected time */
    for (j = 0; j < ntime + 1; j++) {
        sm_param_read(stream, m);
        if (j == modelpar.time) {
            sig = m[modelpar.season].d;
        }
    }

    /* Emission line brightness depending on airmass and layer height
       (van Rhijn 1921) */
    z = (90. - modelpar.alt) * CPL_MATH_RAD_DEG;
    cvr = 1 / sqrt(1 - pow((SM_ERAD / (SM_ERAD + height)) * sin(z), 2));

    /* Influence of solar radio flux [sfu] */
    csol = slope * modelpar.msolflux + cons;

    /* Correction factors for airglow continuum */
    corr = scale * cvr * csol * mean;
    rsig = sig / mean;

    /* Read relative airglow continuum fluxes and uncertainties
       and compute absolute values */
    for (i = 0; i < acont.n; i++) {
        sm_param_read(stream, x);
        acont.dat[i].lam = x[0].d;
        acont.dat[i].flux = x[1].d * corr;
        if (x[1].d == 0) {
            rx = 0;
        } else {
            rx = x[2].d / x[1].d;
        }
        acont.dat[i].dflux1 = acont.dat[i].flux *
                              sqrt(pow(rx, 2) + pow(rsig, 2));
        acont.dat[i].dflux2 = acont.dat[i].dflux1;
    }

    /* Close file */
    fclose(stream);

    /* Linear interpolation of airglow continuum */
    sm_spec_interpol(spec, &acont);
    sm_spec_free(&acont);

    /* Correction for molecular absorption of airglow emission in the lower
       atmosphere */
    sm_spec_copy(&atrans, molabstrans);
    sm_spec_scale(&atrans, '^', modelpar.airmass);
    sm_spec_calc(spec, '*', &atrans);
    sm_spec_free(&atrans);

    /* Extinction correction considering reduced loss of scattered light
       for extended emission compared to point source */

    /* Calculate airglow airmass (based on van Rhijn formula for 90 km) */
    xag = 1 / sqrt(1.0 - 0.972 * pow(sin(z), 2));

    /* Rayleigh scattering (coefficients from 3D scattering calculations) */
    lgx = log10(xag);
    fr = mr * lgx + cr;
    cxr = fr * xag;
    sm_spec_copy(&rtrans, rscattrans);
    sm_spec_scale(&rtrans, '^', cxr);
    sm_spec_calc(spec, '*', &rtrans);
    sm_spec_free(&rtrans);

    /* Mie scattering/aerosol extinction (coeff. from 3D scattering calc.) */
    fm = mm * lgx + cm;
    cxm = fm * xag;
    sm_spec_copy(&mtrans, mscattrans);
    sm_spec_scale(&mtrans, '^', cxm);
    sm_spec_calc(spec, '*', &mtrans);
    sm_spec_free(&mtrans);

    return CPL_ERROR_NONE;
}

/**@}*/
