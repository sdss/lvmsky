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

/**@{*/

/*!
 * skycorr.c
 *
 * Authors:     Stefan Noll & ESO In-Kind Team Innsbruck
 * Created:     01 Aug 2011
 * Last update: 28 Aug 2013
 *
 * Stand-alone programme for sky correction code
 *
 * The programme has one mandatory parameter:
 * (1) path + name of driver file
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <sc_skycorr.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

int main(int argc, char *argv[])
{
    cpl_errorstate errstate = cpl_errorstate_get();
    char parfile[SC_MAXLEN], errtxt[SC_MAXLEN];
    cpl_boolean err = CPL_FALSE;
    double cs, ce;

    cpl_init(CPL_INIT_DEFAULT);
    cpl_msg_info(cpl_func, "Skycorr %s", VERSION);

    /* Read driver file name */
    if (argc > 1) {
        strcpy(parfile, (char *) argv[1]);
    } else {
        sprintf(errtxt, "%s: no parameter file", SC_ERROR_IIP_TXT);
        cpl_error_set_message(cpl_func, SC_ERROR_IIP, "%s", errtxt);
        err = CPL_TRUE;
    }

    /* Print name of parameter file */
    if (err == CPL_FALSE) {
         cpl_msg_info(cpl_func, "Using parameter file: %s", parfile);
    }

    /* Perform sky correction and get run time */
    if (err == CPL_FALSE) {
        cs = cpl_test_get_walltime();
        sc_skycorr(parfile);
        ce = cpl_test_get_walltime();
        cpl_msg_info(cpl_func, "Code run time: %g s", ce - cs);
    }

    /* Show errors */
    if (cpl_errorstate_is_equal(errstate)) {
        cpl_msg_info(cpl_func, "No errors occurred");
    } else {
        cpl_errorstate_dump(errstate, CPL_TRUE, cpl_errorstate_dump_one);
    }

    if (cpl_errorstate_is_equal(errstate)) {
        cpl_end();
        return EXIT_SUCCESS;
    } else {
        cpl_end();
        return EXIT_FAILURE;
    }
}

/**@}*/
