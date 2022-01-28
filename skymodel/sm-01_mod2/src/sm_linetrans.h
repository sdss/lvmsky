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
 * \defgroup sm_module1a Sky Model Module 1a
 *
 * This module provides functions for preparing a library of airglow line
 * transmission files
 */

/*!
 * \ingroup sm_module1a
 */

/**@{*/

/*!
 * \file sm_linetrans.h
 *
 * Header for routines concerning the preparing a library of airglow line
 * transmission files
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  06 Dec 2011
 * \date   06 Oct 2015
 */

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************************
 *                                INCLUDES                                   *
 ****************************************************************************/

/* Config header */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* Sky model headers */

#include <sm_general.h>
#include <sm_skyemcomp.h>

/*****************************************************************************
 *                                DEFINES                                    *
 ****************************************************************************/

#ifndef SM_LINETRANS_H
#define SM_LINETRANS_H

/* Definition of constants */

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code sm_linetrans(const char *path);
cpl_error_code sm_linetrans_getclassdat(cpl_table *classdat,
                                        const smparmodel *modelpar);
cpl_error_code sm_linetrans_getlibstruct(cpl_table *filelist,
                                         const smparmodel *modelpar,
                                         const int libflag);
cpl_error_code sm_linetrans_calclib(const cpl_table *filelist,
                                    const cpl_table *classdat,
                                    const smspec *linetab);

#endif /* SM_LINETRANS_H */

#ifdef __cplusplus
}
#endif

/**@}*/

