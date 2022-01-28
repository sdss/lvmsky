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
 * \callgraph
 *
 * \file test_sc_contsub.c
 *
 * \brief test routine for continuum subtraction library
 *
 * \author Marco Barden, Stefan Noll, & ESO In-Kind Team Innsbruck
 *
 */


/*****************************************************************************
 *                                 DEFINES                                   *
 ****************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <sc_basic.h>
#include <sc_par.h>
#include <sc_contsub.h>
#include <sc_specdiss.h>
#include <sc_readspec.h>
#include <sc_lines.h>
#include <sc_conv.h>
#include <sc_plot.h>

/*****************************************************************************
 *                                    CODE                                   *
 ****************************************************************************/

int main(void)
{
    cpl_errorstate err_state;            // CPL error state
    cpl_parameterlist *parlist;

    cpl_table *spec, *linetab, *groups;

    char parfile[SC_MAXLEN] = "config/sctest_sinfo_H.par";

    double ts = 0;                       // start timer

    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    err_state = cpl_errorstate_get();

    ts = cpl_test_get_cputime();

    // Read parameter file into CPL parameter list
    parlist = cpl_parameterlist_new();
    sc_par_readfile(parlist, parfile);

    // Create line table
    linetab = cpl_table_new(0);
    sc_specdiss_init_linetab(linetab);

    // Create test science spectrum
    spec = cpl_table_new(0);
    groups = cpl_table_new(0);
    if (sc_readspec(spec, parlist) != CPL_ERROR_NONE)
    {
        cpl_msg_error(cpl_func, "sc_readspec() failed");
        cpl_errorstate_dump(err_state, CPL_FALSE, cpl_errorstate_dump_one);
    } else if (sc_lines(groups, parlist) != CPL_ERROR_NONE) {
        cpl_msg_error(cpl_func, "sc_lines() failed");
        cpl_errorstate_dump(err_state, CPL_FALSE, cpl_errorstate_dump_one);
    } else if (sc_specdiss_find_emissionlines(spec, linetab, parlist,
                                              groups) != CPL_ERROR_NONE) {
        cpl_msg_error(cpl_func, "sc_specdiss_find_emissionlines() failed");
        cpl_errorstate_dump(err_state, CPL_FALSE, cpl_errorstate_dump_one);
    } else if (sc_contsub(spec) != CPL_ERROR_NONE) {
        cpl_msg_error(cpl_func, "sc_contsub() failed");
        cpl_errorstate_dump(err_state, CPL_FALSE, cpl_errorstate_dump_one);
    } else {
        cpl_test_rel(cpl_table_get_double(spec, "cflux", 10, NULL)+
                     cpl_table_get_double(spec, "lflux", 10, NULL),
                     cpl_table_get_double(spec, "flux", 10, NULL), 1e-10);
        cpl_test_rel(cpl_table_get_double(spec, "cflux", 38, NULL)+
                     cpl_table_get_double(spec, "lflux", 38, NULL),
                     cpl_table_get_double(spec, "flux", 38, NULL), 1e-10);
    }

//    cpl_table_dump(spec, 0, 50, NULL);

    // cleanup
    cpl_table_delete(linetab);
    cpl_table_delete(spec);
    cpl_table_delete(groups);
    cpl_parameterlist_delete(parlist);

    // CPU-time
    printf("run time: %g s\n", cpl_test_get_cputime() - ts);

    cpl_memory_dump();

//    return cpl_test_end(0);
//    if (cpl_errorstate_is_equal(err_state)) {
//        cpl_end();
//        return EXIT_SUCCESS;
//    } else {
//        cpl_end();
//        return EXIT_FAILURE;
//    }
    return (0);
}

/**@}*/
