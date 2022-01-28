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
 * \ingroup sm_module1a
 */

/**@{*/

/*!
 * \file preplinetrans.c
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  06 Dec 2011
 * \date   06 Oct 2015
 *
 * Programme for preparing a library of airglow line transmission files
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <sm_general.h>
#include <sm_skyemcomp.h>
#include <sm_linetrans.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

int main(int argc, char *argv[])
{
    /*!
     * Main routine for preparing a library of airglow transmission files.
     * It has to be started from the base directory of the sky model code by
     * typing "bin/preplinetrans". The programme has one optional parameter
     * for the path to the file ::SM_FILENAMELIST. The default input is
     * "data".
     */

    cpl_errorstate errstate = cpl_errorstate_get();
    char path[FILENAME_MAX], filenamelist[FILENAME_MAX];
    double cs, ce;

    cpl_init(CPL_INIT_DEFAULT);

    /* Read driver file name */
    if (argc > 1) {
        strcpy(path, (char *) argv[1]);
    } else {
        strcpy(path, "data");
    }

    /* Print path to file for file names */
    cpl_msg_info(cpl_func, "Path to %s: %s", SM_FILENAMELIST, path);

    /* Check existence of path */
    if (sm_basic_access(path, F_OK) != CPL_ERROR_NONE) {
        cpl_errorstate_dump(errstate, CPL_FALSE, cpl_errorstate_dump_one);
        cpl_end();
        return EXIT_FAILURE;
    }

    /* Build full name of file for file names */
    sprintf(filenamelist, "%s/%s", path, SM_FILENAMELIST);

    /* Check existence of file for file names */
    if (sm_basic_access(filenamelist, F_OK) != CPL_ERROR_NONE) {
        cpl_errorstate_dump(errstate, CPL_FALSE, cpl_errorstate_dump_one);
        cpl_end();
        return EXIT_FAILURE;
    }

    /* Prepare lists of line transmission and get run time */
    cs = cpl_test_get_walltime();
    sm_linetrans(path);
    ce = cpl_test_get_walltime();
    cpl_msg_info(cpl_func, "Code run time: %.2f min", (ce - cs) / 60.);

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
