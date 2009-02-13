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

#ifndef PAGMO_ATOMIC_COUNTER_GCC_41_H
#define PAGMO_ATOMIC_COUNTER_GCC_41_H

#include "base_atomic_counter.h"

namespace PaGMO
{
	template <class IntType>
	class atomic_counter_gcc_41: public base_atomic_counter<IntType, atomic_counter_gcc_41<IntType> >
	{
			typedef base_atomic_counter<IntType, atomic_counter_gcc_41<IntType> > ancestor;
		public:
			atomic_counter_gcc_41():ancestor::base_atomic_counter() {}
			template <class IntType2>
			atomic_counter_gcc_41(const IntType2 &n):ancestor::base_atomic_counter(n) {}
			template <class IntType2>
			atomic_counter_gcc_41 &operator+=(const IntType2 &n) {
				__sync_add_and_fetch(&(this->m_value),static_cast<IntType>(n));
				return *this;
			}
			template <class IntType2>
			atomic_counter_gcc_41 &operator-=(const IntType2 &n) {
				__sync_sub_and_fetch(&(this->m_value),static_cast<IntType>(n));
				return *this;
			}
			atomic_counter_gcc_41 &operator++() {
				return operator+=(static_cast<IntType>(1));
			}
			atomic_counter_gcc_41 &operator--() {
				return operator-=(static_cast<IntType>(1));
			}
			atomic_counter_gcc_41 operator++(int) {
				return atomic_counter_gcc_41(__sync_fetch_and_add(&(this->m_value),1));
			}
			atomic_counter_gcc_41 operator--(int) {
				return atomic_counter_gcc_41(__sync_fetch_and_sub(&(this->m_value),1));
			}
			template <class IntType2>
			bool compare_and_swap(const IntType2 &oldval, const IntType2 &newval) {
				return __sync_bool_compare_and_swap(&(this->m_value),oldval,newval);
			}
	};
}

#endif
