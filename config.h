/*****************************************************************************
 *   Copyright (C) 2008, 2009 Advanced Concepts Team (European Space Agency) *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#ifndef PAGMO_CONFIG_H
#define PAGMO_CONFIG_H

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 100000 \
	+ __GNUC_MINOR__ * 1000 \
	+ __GNUC_PATCHLEVEL__ * 10)
#else
#error "The GCC compiler is needed to compile PaGMO."
#endif

#ifdef PAGMO_WIN32
#ifdef PAGMO_DLL_EXPORT_API
#define __PAGMO_VISIBLE __declspec(dllexport)
#elif defined ( PAGMO_DLL_IMPORT_API )
#define __PAGMO_VISIBLE __declspec(dllimport)
#else
#define __PAGMO_VISIBLE
#endif
#define __PAGMO_VISIBLE_FUNC __PAGMO_VISIBLE
#else
#define __PAGMO_VISIBLE __attribute__ ((visibility("default")))
#define __PAGMO_VISIBLE_FUNC
#endif

#endif
