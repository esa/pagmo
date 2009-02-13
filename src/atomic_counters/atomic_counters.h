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

// 26/12/08 Created by Francesco Biscani.

#ifndef PAGMO_ATOMIC_COUNTERS_H
#define PAGMO_ATOMIC_COUNTERS_H

#include "../config.h"

#if defined( __GNUC__ ) && GCC_VERSION >= 401000

#include "atomic_counter_gcc_41.h"

namespace PaGMO
{
	typedef atomic_counter_gcc_41<int> atomic_counter_int;
	typedef atomic_counter_gcc_41<size_t> atomic_counter_size_t;
}

#else // Not GCC or GCC < 4.1.

// TODO: for MSVC, use its atomic builtins instead of the generic counter.
#include "atomic_counter_generic.h"

namespace PaGMO
{
	typedef atomic_counter_generic<int> atomic_counter_int;
	typedef atomic_counter_generic<size_t> atomic_counter_size_t;
}

#endif // Compiler selection in case of MT.

#endif
