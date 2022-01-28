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
 * \file test_sc_specdiss.c
 *
 * test routine for spectral analyse library
 *
 * \author Wolfgang Kausch, Stefan Noll, & ESO In-Kind Team Innsbruck
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
#include <sc_lines.h>
#include <sc_specdiss.h>
#include <sc_readspec.h>
#include <sc_conv.h>
#include <sc_plot.h>

/*****************************************************************************
 *                                    CODE                                   *
 ****************************************************************************/

int main(int argc, char *argv[])
{

    double cs=0, ce=0;

    /* CPL test tables */
    cpl_table *spectrum1;
    cpl_table *groups;
    cpl_errorstate error_state=CPL_ERROR_NONE; /* CPL error state           */
    cpl_error_code errstat;

    /* Reading parameter file */
    char parfile[SC_MAXLEN] = "config/sctest_sinfo_H.par";
    cpl_parameterlist *parlist;
    cpl_table *linetab;

/*--------------------------------------------------------------------------*/

    /* Initialising */
    cpl_init(CPL_INIT_DEFAULT);
//     cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    error_state = cpl_errorstate_get();
    cs = cpl_test_get_cputime();
    cpl_msg_set_time_off();
    cpl_msg_set_component_on();
    cpl_msg_set_domain("[test_sc_specdiss]");
    cpl_msg_set_domain_on();
    cpl_msg_set_log_name("test_sc_specdiss.log");
    cpl_msg_set_level(CPL_MSG_INFO);
    cpl_msg_set_log_level(CPL_MSG_DEBUG);

    /* Read driver file name if given */
    if (argc > 1) {
        strcpy(parfile, (char *) argv[1]);
    }

    /* Creating parameter list */
    parlist = cpl_parameterlist_new();

    /* Filling parameter list */
    sc_par_readfile(parlist, parfile);

    /* creating test science spectrum */
    spectrum1 = cpl_table_new(0);
    if ((errstat = sc_readspec(spectrum1, parlist)) != CPL_ERROR_NONE) {
        cpl_parameterlist_delete(parlist);
        cpl_table_delete(spectrum1);
        cpl_errorstate_dump(error_state, CPL_FALSE, cpl_errorstate_dump_one);
        cpl_end();
        return EXIT_FAILURE;
    }

    /* Read line groups from ASCII file and adapt it to observing
       conditions */
    groups = cpl_table_new(0);
    if ((errstat = sc_lines(groups, parlist)) != CPL_ERROR_NONE) {
        cpl_parameterlist_delete(parlist);
        cpl_table_delete(spectrum1);
        cpl_table_delete(groups);
        cpl_end();
        return EXIT_FAILURE;
    }

/*--------------------------------------------------------------------------*/

    /* Creating line property table */
    linetab=cpl_table_new(0);
    sc_specdiss_init_linetab(linetab);

    /* Spectra dissection stuff */
    sc_specdiss_find_emissionlines(spectrum1, linetab, parlist, groups);

//     cpl_table_dump(spectrum1,0,4,NULL);
/*--------------------------------------------------------------------------*/

//   cpl_table_dump(linetab,0,10,NULL);
//   cpl_table_dump_structure(spectrum1,NULL);

    /* Cleaning */
    cpl_parameterlist_delete(parlist);
    cpl_table_delete(spectrum1);
    cpl_table_delete(groups);
    cpl_table_delete(linetab);

    printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    printf("'%s' routine end memory dump:\n",cpl_func);
    cpl_memory_dump();
    printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");


    ce = cpl_test_get_cputime();
    printf("run time: %g s\n", ce - cs);
    cpl_errorstate_dump(error_state, CPL_FALSE, cpl_errorstate_dump_one);

// return cpl_test_end(0);

    if (cpl_errorstate_is_equal(error_state)) {
        cpl_end();
        return EXIT_SUCCESS;
    } else {
        cpl_end();
        return EXIT_FAILURE;
    }

}
