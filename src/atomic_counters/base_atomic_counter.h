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

#ifndef PAGMO_BASE_ATOMIC_COUNTER_H
#define PAGMO_BASE_ATOMIC_COUNTER_H

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace PaGMO
{
	template <class IntType, class Derived>
	class base_atomic_counter
	{
		public:
			base_atomic_counter():m_value(0) {}
			template <class IntType2>
			base_atomic_counter(const IntType2 &n):m_value(static_cast<IntType2>(n)) {}
			operator IntType() const {
				return m_value;
			}
		protected:
			IntType m_value;
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
