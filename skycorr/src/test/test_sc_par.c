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
 * \file test_sc_par.c
 *
 * test routines for handling parameter files
 *
 * \author Wolfgang Kausch & ESO In-Kind Team Innsbruck
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

#include <sc_par.h>
#include <sc_basic.h>


/*****************************************************************************
 *                                    CODE                                   *
 ****************************************************************************/

int main(void)
{
    cpl_errorstate error_state;         /* CPL error state                  */

    cpl_parameterlist *parlist;         /* Parameter list to be filled      */
    cpl_parameter *p;

    char parfile_in[SC_MAXLEN];         /* Name of paramfile to be read     */
    char parfile_out[SC_MAXLEN];        /* Name of paramfile to be written  */
    char sys[SC_MAXLEN];                /* Copy command                     */

/*--------------------------------------------------------------------------*/

    /* Initialising */
    cpl_init(CPL_INIT_DEFAULT);
    error_state = cpl_errorstate_get();
    cpl_msg_set_time_off();
    cpl_msg_set_component_on();
    cpl_msg_set_domain("[test_sc_par]");
    cpl_msg_set_domain_on();
    cpl_msg_set_log_name("test_sc_par.log");
    cpl_msg_set_level(CPL_MSG_INFO);
    cpl_msg_set_log_level(CPL_MSG_DEBUG);

    /* Prepare strings */
    sprintf(parfile_in, "%s", "config/sctest_sinfo_H.par");
    sprintf(parfile_out, "%s", "../../output/test.par");
    sprintf(sys,"cp %s %s", parfile_in, parfile_out);

    /* Read parameter file */
    parlist = cpl_parameterlist_new();
    sc_par_readfile(parlist, parfile_in);

    /* Change content of parameter 'base_dir' */
    p = cpl_parameterlist_find(parlist, "inst_dir");
    cpl_parameter_set_string(p, "/nirvana");

    /* Change content of parameter 'col_dflux' */
    p = cpl_parameterlist_find(parlist, "col_dflux");
    cpl_parameter_set_string(p, "dflux");

    /* Change content of parameter 'fwhm' */
    p = cpl_parameterlist_find(parlist, "fwhm");
    cpl_parameter_set_double(p, 4.0);

    /* Change content of parameter 'rebintype' */
    p = cpl_parameterlist_find(parlist, "rebintype");
    cpl_parameter_set_int(p, 0);

    /* Change content of parameter 'plot_type' */
    p = cpl_parameterlist_find(parlist, "plot_type");
    cpl_parameter_set_string(p, "X");

    /* Write output parameter file */
    if (system(sys)) {};
    sc_par_writefile(parlist, parfile_out);

    /* Errors and run time */
    cpl_parameterlist_delete(parlist);
    cpl_errorstate_dump(error_state, CPL_FALSE, cpl_errorstate_dump_one);
    cpl_memory_dump();

    if (cpl_errorstate_is_equal(error_state)) {
        cpl_end();
        return EXIT_SUCCESS;
    } else {
        cpl_end();
        return EXIT_FAILURE;
    }
}

/**@}*/
