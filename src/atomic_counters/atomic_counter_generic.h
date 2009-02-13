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

#ifndef PAGMO_ATOMIC_COUNTER_GENERIC_H
#define PAGMO_ATOMIC_COUNTER_GENERIC_H

#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>

#include "base_atomic_counter.h"

namespace PaGMO
{
	template <class IntType>
	class atomic_counter_generic: public base_atomic_counter<IntType,atomic_counter_generic<IntType> >
	{
			typedef base_atomic_counter<IntType,atomic_counter_generic<IntType> > ancestor;
			typedef boost::lock_guard<boost::mutex> lock_type;
		public:
			atomic_counter_generic():ancestor::base_atomic_counter(),m_mutex() {}
			template <class IntType2>
			atomic_counter_generic(const IntType2 &n):ancestor::base_atomic_counter(n),m_mutex() {}
			atomic_counter_generic(const atomic_counter_generic &a):ancestor::base_atomic_counter(a.m_value) {}
			template <class IntType2>
			atomic_counter_generic &operator+=(const IntType2 &n) {
				lock_type lock(m_mutex);
				this->m_value += n;
				return *this;
			}
			template <class IntType2>
			atomic_counter_generic &operator-=(const IntType2 &n) {
				lock_type lock(m_mutex);
				this->m_value -= n;
				return *this;
			}
			atomic_counter_generic &operator++() {
				return operator+=(static_cast<IntType>(1));
			}
			atomic_counter_generic &operator--() {
				return operator-=(static_cast<IntType>(1));
			}
			atomic_counter_generic operator++(int) {
				lock_type lock(m_mutex);
				atomic_counter_generic retval(this->m_value);
				++this->m_value;
				return retval;
			}
			atomic_counter_generic operator--(int) {
				lock_type lock(m_mutex);
				atomic_counter_generic retval(this->m_value);
				--this->m_value;
				return retval;
			}
			template <class IntType2>
			bool compare_and_swap(const IntType2 &oldval, const IntType2 &newval) {
				lock_type lock(m_mutex);
				if (this->m_value == oldval) {
					this->m_value = newval;
					return true;
				}
				return false;
			}
		private:
			boost::mutex m_mutex;
	};
}

#endif
