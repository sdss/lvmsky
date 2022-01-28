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

/*
 * test_sc_skycorr.c
 *
 * Authors:     Stefan Noll & ESO In-Kind Team Innsbruck
 * Created:     22 Feb 2011
 * Last update: 28 Aug 2013
 *
 * Test programme for sky correction code
 *
 * The programme has one optional parameter:
 * (1) path + name of driver file
 *     (default: config/sctest_sinfo_H.par)
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
    cpl_errorstate prev_state = cpl_errorstate_get();
    char parfile[SC_MAXLEN] = "config/sctest_sinfo_H.par";
    double cs, ce;

    cpl_init(CPL_INIT_DEFAULT);

    /* Read driver file name if given */
    if (argc > 1) {
        strcpy(parfile, (char *) argv[1]);
    }

    /* Print name of parameter file */
    cpl_msg_info(cpl_func, "Driver file: %s", parfile);

    /* Perform sky correction and get run time */
    cs = cpl_test_get_walltime();
    sc_skycorr(parfile);
    ce = cpl_test_get_walltime();
    cpl_msg_info(cpl_func, "Code run time: %g s", ce - cs);

    /* Show errors */
    cpl_errorstate_dump(prev_state, CPL_TRUE, cpl_errorstate_dump_one);

    if (cpl_errorstate_is_equal(prev_state)) {
        cpl_end();
        return EXIT_SUCCESS;
    } else {
        cpl_end();
        return EXIT_FAILURE;
    }
}

/**@}*/
