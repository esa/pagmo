/***************************************************************************
 *   Copyright (C) 2009 by Francesco Biscani   *
 *   bluescarni@gmail.com   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef ATOMIC_COUNTER_MSVC_LONG_H
#define ATOMIC_COUNTER_MSVC_LONG_H

#include <boost/numeric/conversion/cast.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_same.hpp>
#include <windows.h>

/// Atomic counter class for Visual Studio C++.
/**
 * Will use MSVC's atomic builtins.
 */
class atomic_counter_msvc_long
{
	public:
		/// Default constructor.
		/**
		 * It will initialise the internal value to zero.
		 */
		atomic_counter_msvc_long():m_value(0) {}
		/// Constructor from generic integer type.
		/**
		 * IntType2 must be able to be used to construct long.
		 * @throws boost::numeric::bad_numeric_cast if long and IntType2 are not the same type
		 * and conversion from IntType2 to long results in loss of information.
		 */
		template <class IntType2>
		atomic_counter_msvc_long(const IntType2 &n):m_value(n)
		{
			BOOST_STATIC_ASSERT(boost::is_integral<IntType2>::value);
			if (!boost::is_same<long,IntType2>::value) {
				long tmp(boost::numeric_cast<long>(n));
				(void)tmp;
			}
		}
		/// Copy constructor.
		/**
		 * Operation is atomic with respect to input argument a.
		 */
		atomic_counter_msvc_long(const atomic_counter_msvc_long &a):m_value(a.get_value()) {}
		/// Assignment operator.
		/**
		 * This operation acquires atomically the internal value of a, but the assignment
		 * to this is not atomic.
		 */
		atomic_counter_msvc_long &operator=(const atomic_counter_msvc_long &a)
		{
			if (this != &a) {
				m_value = a.get_value();
			}
			return *this;
		}
		/// In-place addition.
		atomic_counter_msvc_long &operator+=(const long &n)
		{
			InterlockedExchangeAdd(&m_value,n);
			return *this;
		}
		/// In-place subtraction.
		atomic_counter_msvc_long &operator-=(const long &n)
		{
			InterlockedExchangeAdd(&m_value,-n);
			return *this;
		}
		/// Prefix increment.
		atomic_counter_msvc_long &operator++()
		{
			InterlockedIncrement(&m_value);
			return *this;
		}
		/// Prefix decrement.
		atomic_counter_msvc_long &operator--()
		{
			InterlockedDecrement(&m_value);
			return *this;
		}
		/// Postfix increment.
		atomic_counter_msvc_long operator++(int)
		{
			return atomic_counter_msvc_long(InterlockedExchangeAdd(&m_value,1));
		}
		/// Postfix decrement.
		atomic_counter_msvc_long operator--(int)
		{
			return atomic_counter_msvc_long(InterlockedExchangeAdd(&m_value,-1));
		}
		/// Get copy of internal value.
		long get_value() const
		{
			return InterlockedExchangeAdd(&m_value,0);
		}
		/// Fast type-trait for increment operation.
		static const bool is_increment_fast = true;
		/// Fast type-trait for arithmetics.
		static const bool is_arithmetics_fast = true;
	private:
		/// Internal value.
		/**
		 * Declared mutable because atomic_counter_msvc_long::get_value needs to perform the operation
		 * this + 0 in order to fetch safely the current m_value with MSVC's atomic builtins.
		 */
		mutable long m_value;
};

#endif
 