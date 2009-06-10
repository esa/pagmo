/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
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

// 26/12/08 Created by Francesco Biscani.

#ifndef PAGMO_ATOMIC_COUNTERS_H
#define PAGMO_ATOMIC_COUNTERS_H

#if defined ( _PAGMO_GCC_ATOMIC_BUILTINS )

#include "atomic_counter_gcc_41.h"

typedef atomic_counter_gcc_41<size_t> atomic_counter_size_t;

#elif defined ( _PAGMO_MSVC_ATOMIC_BUILTINS )

#include "atomic_counter_msvc_long.h"

// TODO: here for win64 bit we probably need another counter altogether and another #ifdef, since
// MSVC's atomic builtins operate on 32bit and 64bit with different naming conventions.
typedef atomic_counter_msvc_long atomic_counter_size_t;

#else

#include "atomic_counter_generic.h"

typedef atomic_counter_generic<size_t> atomic_counter_size_t;

#endif

#endif
