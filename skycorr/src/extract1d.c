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
 * extract1d.c
 *
 * Authors:     Stefan Noll & ESO In-Kind Team Innsbruck
 * Created:     28 Apr 2013
 * Last update: 28 Aug 2013
 *
 * Programme for extracting a 1D science and sky reference spectrum from two
 * 2D FITS images
 *
 * The programme has five input parameters:
 * (1) path + name of driver file (provides path + name of resulting sky 1D
 *     FITS image and the FITS extension names)
 * (2) path + name of input 2D FITS image for science spectrum
 * (3) path + name of input 2D FITS image for reference sky spectrum
 *     (can be the same as (2))
 * (4) minimum flux in spatial profile function of science spectrum relative
 *     to profile maximum to select pixel for extraction (range: [0,1];
 *     default: 0.01)
 * (5) fraction of pixels with flux below selected pixel in sky spectrum
 *     (range: [0,1]; default: 0.5 = median); pixel selection for median
 *     along spatial direction for each wavelength separately
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <sc_conv.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

int main(int argc, char *argv[])
{
    cpl_errorstate errstate = cpl_errorstate_get();
    char parfile[SC_MAXLEN], inscifile[SC_MAXLEN];
    char inskyfile[SC_MAXLEN], errtxt[SC_MAXLEN];
    cpl_boolean err = CPL_FALSE;
    double minrelprofflux = 0.01; // default value
    double lowpixfrac = 0.5; // default value

    cpl_init(CPL_INIT_DEFAULT);

    /* Read input parameters */
    if (argc >= 4) {
        strcpy(parfile, (char *) argv[1]);
        strcpy(inscifile, (char *) argv[2]);
        strcpy(inskyfile, (char *) argv[3]);
        if (argc == 5) {
            if (sc_basic_isnumber((char *) argv[4]) == CPL_FALSE) {
                sprintf(errtxt, "%s: minrelprofflux (no number)",
                        SC_ERROR_IIP_TXT);
                cpl_error_set_message(cpl_func, SC_ERROR_IIP, "%s", errtxt);
                err = CPL_TRUE;
            }
            minrelprofflux = atof(argv[4]);
            if (minrelprofflux < 0. || minrelprofflux > 1.) {
                sprintf(errtxt, "%s: minrelprofflux (value outside [0,1])",
                        SC_ERROR_IIP_TXT);
                cpl_error_set_message(cpl_func, SC_ERROR_IIP, "%s", errtxt);
                err = CPL_TRUE;
            }
        }
        if (argc >= 6) {
            if (sc_basic_isnumber((char *) argv[5]) == CPL_FALSE) {
                sprintf(errtxt, "%s: lowpixfrac (no number)",
                        SC_ERROR_IIP_TXT);
                cpl_error_set_message(cpl_func, SC_ERROR_IIP, "%s", errtxt);
                err = CPL_TRUE;
            }
            lowpixfrac = atof(argv[5]);
            if (lowpixfrac < 0. || lowpixfrac > 1.) {
                sprintf(errtxt, "%s: lowpixfrac (value outside [0,1])",
                        SC_ERROR_IIP_TXT);
                cpl_error_set_message(cpl_func, SC_ERROR_IIP, "%s", errtxt);
                err = CPL_TRUE;
            }
        }
    } else {
        if (argc <= 1) {
            sprintf(errtxt, "%s: no parameter file", SC_ERROR_IIP_TXT);
            cpl_error_set_message(cpl_func, SC_ERROR_IIP, "%s", errtxt);
            err = CPL_TRUE;
        }
        if (argc <= 2) {
            sprintf(errtxt, "%s: no input 2D FITS images", SC_ERROR_IIP_TXT);
            cpl_error_set_message(cpl_func, SC_ERROR_IIP, "%s", errtxt);
            cpl_end();
            err = CPL_TRUE;
        } else {
            sprintf(errtxt, "%s: no input 2D FITS sky image",
                    SC_ERROR_IIP_TXT);
            cpl_error_set_message(cpl_func, SC_ERROR_IIP, "%s", errtxt);
            cpl_end();
            err = CPL_TRUE;
        }
    }

    /* Print name of parameter file */
    if (err == CPL_FALSE) {
        cpl_msg_info(cpl_func, "Using parameter file: %s", parfile);
    }

    /* Extract science and sky reference 1D FITS images from 2D FITS images */
    if (err == CPL_FALSE) {
        sc_conv_extract1d(parfile, inscifile, inskyfile, minrelprofflux,
                          lowpixfrac);
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
