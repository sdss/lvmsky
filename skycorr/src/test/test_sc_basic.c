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
 * \file test_sc_basic.c
 *
 * \brief test routine for basic library
 *
 * \author Marco Barden & ESO In-Kind Team Innsbruck
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

/*****************************************************************************
 *                                    CODE                                   *
 ****************************************************************************/

int main(void)
{
    cpl_errorstate err_state;            // CPL error state
    cpl_array *arr;
    cpl_table *tab;

    double mean = 0, rms = 0;
    double data[] = { 4.5, 4.9, 5.6, 4.2, 6.2, 5.2, 9.9, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    int ndata = 20;

    char targetcol[SC_MAXLEN];

    int i;

    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    err_state = cpl_errorstate_get();

    arr = cpl_array_new(ndata, CPL_TYPE_DOUBLE);
    for (i = 0; i < ndata; i++) {
        cpl_array_set_double(arr, i, data[i]);
        if (data[i] == 0.0) {
            cpl_array_set_invalid(arr, i);
        }
    }

    /* Test of sc_basic_clipmean */
    if (sc_basic_clipmean(&mean, &rms, arr, CPL_FALSE) != CPL_ERROR_NONE) {
        cpl_msg_error(cpl_func, "sc_basic_clipmean() failed");
        cpl_errorstate_dump(err_state, CPL_FALSE, cpl_errorstate_dump_one);
    } else {
//        printf("mean: %lf +- %lf\n", mean, rms);
        cpl_test_rel(mean, 5.34, 0.02);
        cpl_test_rel(rms, 1.04, 0.02);
    }

    /* Test of sc_basic_copy_array2tablecol                                 */
    tab = cpl_table_new(ndata);
    sprintf(targetcol,"%s","target_col");
    if (sc_basic_arr2col(tab,targetcol,arr) != CPL_ERROR_NONE) {
        cpl_msg_error(cpl_func, "sc_basic_copy_array2tablecol() failed");
        cpl_errorstate_dump(err_state, CPL_FALSE, cpl_errorstate_dump_one);
    }
//     cpl_table_dump_structure(tab,NULL);
//     cpl_table_dump(tab,0,19,NULL);

    /* Test of sc_basic_copy_tablecol2array                                 */
    cpl_array_delete(arr);
    arr = cpl_array_new(0, CPL_TYPE_DOUBLE);
    if (sc_basic_col2arr(arr,tab,targetcol) != CPL_ERROR_NONE) {
        cpl_msg_error(cpl_func, "sc_basic_copy_tablecol2array() failed");
        cpl_errorstate_dump(err_state, CPL_FALSE, cpl_errorstate_dump_one);
    }
//     cpl_array_dump_structure(arr,NULL);
//     cpl_array_dump(arr,0,19,NULL);

    // cleanup
    cpl_array_delete(arr);
    cpl_table_delete(tab);

    // test sc_basic_sort
    arr = cpl_array_new(5, CPL_TYPE_DOUBLE);
    cpl_array_set(arr, 0, 3.);
    cpl_array_set(arr, 1, 7.);
    cpl_array_set(arr, 2, 4.);
    cpl_array_set(arr, 3, 9.);
    cpl_array_set(arr, 4, 4.);

    sc_basic_sortarr(arr, CPL_TYPE_DOUBLE);

//    cpl_array_dump(arr, 0, 5, NULL);

    // cleanup
    cpl_array_delete(arr);

    cpl_memory_dump();

    return cpl_test_end(0);
}

/**@}*/
