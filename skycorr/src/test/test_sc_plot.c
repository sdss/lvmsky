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
 * \ingroup molecfit
 */

/**@{*/

/*!
 * \callgraph
 *
 * \file test_sc_plot.c
 *
 * test routine for plotting library
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

#include <sc_basic.h>
#include <sc_par.h>
#include <sc_readspec.h>
#include <sc_conv.h>
#include <sc_plot.h>

/*****************************************************************************
 *                                    CODE                                   *
 ****************************************************************************/

int main(void)
{
    cpl_errorstate error_state=CPL_ERROR_NONE;  /* CPL error state          */
    cpl_error_code errstat;

    /* CPL test tables */
    cpl_table *spectrum1, *spectrum2;

    /* Reading parameter file */
    char parfile[SC_MAXLEN];
    cpl_parameter *p;
    cpl_parameterlist *parlist;
    cpl_parameterlist *plottags;

/*--------------------------------------------------------------------------*/
    /* Initialising */
    cpl_init(CPL_INIT_DEFAULT);
    error_state = cpl_errorstate_get();
    cpl_msg_set_time_off();
    cpl_msg_set_component_on();
    cpl_msg_set_domain("[test_sc_plot]");
    cpl_msg_set_domain_on();
    cpl_msg_set_log_name("test_sc_plot.log");
    cpl_msg_set_level(CPL_MSG_INFO);
    cpl_msg_set_log_level(CPL_MSG_DEBUG);

    /* Reading two test spectra */
    sprintf(parfile,"%s","config/sctest_sinfo_H.par");

    /* Creating parameter list */
    parlist = cpl_parameterlist_new();
    /* Filling parameter list */
    sc_par_readfile(parlist,parfile);

    /* creating test OBJ/SKY spectra */
    p = cpl_parameterlist_find(parlist, "spectype");
    cpl_parameter_set_string(p, "SCI");
    spectrum1=cpl_table_new(0);
    if ((errstat = sc_readspec(spectrum1, parlist)) !=
        CPL_ERROR_NONE)
    {
        cpl_parameterlist_delete(parlist);
        cpl_table_delete(spectrum1);
    }

    p = cpl_parameterlist_find(parlist, "spectype");
    cpl_parameter_set_string(p, "SKY");
    spectrum2=cpl_table_new(0);
    if ((errstat = sc_readspec(spectrum2, parlist)) !=
        CPL_ERROR_NONE)
    {
        cpl_parameterlist_delete(parlist);
        cpl_table_delete(spectrum2);
    }

/* ------------------------------------------------------------------------ */

    /* Single plot */
    plottags=cpl_parameterlist_new();
    sc_setplottags_single_spec(plottags,"lambda","flux","Title",
                               "X-label","Y-label",parlist);
    sc_plot_single_spec(spectrum1,plottags);
    cpl_parameterlist_delete(plottags);

    /* Double plot */
    plottags=cpl_parameterlist_new();
    sc_setplottags_double_spec(plottags,"lambda","flux","Title1",
                               "X-label1","Y-label1","lambda","flux",
                               "Title2","X-label2","Y-label2",parlist);
    sc_plot_double_spec(spectrum1,spectrum2,plottags);
    cpl_parameterlist_delete(plottags);

    /* over plot */
    plottags=cpl_parameterlist_new();
    sc_setplottags_overplot_spec(plottags,"lambda","flux","flux2",
                                 "Spectrum1","Spectrum2",
                                 "Title1","X-label","Y-label1",parlist);
    cpl_table_duplicate_column(spectrum1,"flux2",spectrum2,"flux");
    sc_overplot_spec(spectrum1,plottags);
    cpl_parameterlist_delete(plottags);

/* ------------------------------------------------------------------------ */

    /* Cleaning */
    cpl_parameterlist_delete(parlist);
    cpl_table_delete(spectrum1);
    cpl_table_delete(spectrum2);

    cpl_memory_dump();

    if (cpl_errorstate_is_equal(error_state)) {
        cpl_end();
        return EXIT_SUCCESS;
    } else {
        cpl_end();
        return EXIT_FAILURE;
    }

}
