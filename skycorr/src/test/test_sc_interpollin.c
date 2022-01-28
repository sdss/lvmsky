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
 * \file test_sc_interpollin.c
 *
 * \brief test routine for linear interpolation library
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
    //cpl_errorstate err_state;            // CPL error state

    cpl_error_code err_code;

    double xr[2] = {       2,    4 },
           yr[2] = {       1,    2 },
           xo[7] = { 0, 1, 2, 3, 4, 5, 6 },
           yo[7];

    double ts = 0;                       // start timer

    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    //err_state = cpl_errorstate_get();

    ts = cpl_test_get_cputime();

    err_code = sc_basic_interpollin(xo, yo, 7, xr, yr, 2, 0);

    printf("x(out): %5lf %5lf %5lf %5lf %5lf %5lf %5lf\n",
           xo[0], xo[1], xo[2], xo[3], xo[4], xo[5], xo[6]);
    printf("y(out): %5lf %5lf %5lf %5lf %5lf %5lf %5lf\n",
           yo[0], yo[1], yo[2], yo[3], yo[4], yo[5], yo[6]);
    printf("error_code: %i\n", err_code);

    // CPU-time
    printf("run time: %g s\n", cpl_test_get_cputime() - ts);

    cpl_memory_dump();

    return cpl_test_end(0);
//    if (cpl_errorstate_is_equal(err_state)) {
//        cpl_end();
//        return EXIT_SUCCESS;
//    } else {
//        cpl_end();
//        return EXIT_FAILURE;
//    }
//    return (0);
}
/**@}*/
