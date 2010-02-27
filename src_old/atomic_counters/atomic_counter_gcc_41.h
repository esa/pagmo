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

#ifndef ATOMIC_COUNTER_GCC_41_H
#define ATOMIC_COUNTER_GCC_41_H

#include <boost/numeric/conversion/cast.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_same.hpp>

/// Atomic counter class for GCC >= 4.1.
/**
 * Will use GCC's atomic builtins, available from version 4.1.
 *
 * If IntType is not an integral type, a compile-time error will be produced upon instantiation.
 */
template <class IntType>
class atomic_counter_gcc_41
{
		BOOST_STATIC_ASSERT(boost::is_integral<IntType>::value);
	public:
		/// Default constructor.
		/**
		 * It will initialise the internal value to zero.
		 */
		atomic_counter_gcc_41():m_value(0) {}
		/// Constructor from generic integer type.
		/**
		 * IntType must be constructible from IntType2.
		 * @throws boost::numeric::bad_numeric_cast if IntType and IntType2 are not the same type
		 * and conversion from IntType2 to IntType results in loss of information.
		 */
		template <class IntType2>
		atomic_counter_gcc_41(const IntType2 &n):m_value(n)
		{
			BOOST_STATIC_ASSERT(boost::is_integral<IntType2>::value);
			if (!boost::is_same<IntType,IntType2>::value) {
				IntType tmp(boost::numeric_cast<IntType>(n));
				(void)tmp;
			}
		}
		/// Copy constructor.
		/**
		 * Operation is atomic with respect to input argument a.
		 */
		atomic_counter_gcc_41(const atomic_counter_gcc_41 &a):m_value(a.get_value()) {}
		/// Assignment operator.
		/**
		 * This operation acquires atomically the internal value of a, but the assignment
		 * to this is not atomic.
		 */
		atomic_counter_gcc_41 &operator=(const atomic_counter_gcc_41 &a)
		{
			if (this != &a) {
				m_value = a.get_value();
			}
			return *this;
		}
		/// In-place addition.
		/**
		 * The operation will be forwarded to the internal value.
		 */
		template <class IntType2>
		atomic_counter_gcc_41 &operator+=(const IntType2 &n)
		{
			__sync_add_and_fetch(&m_value,static_cast<IntType>(n));
			return *this;
		}
		/// In-place subtraction.
		/**
		 * The operation will be forwarded to the internal value.
		 */
		template <class IntType2>
		atomic_counter_gcc_41 &operator-=(const IntType2 &n)
		{
			__sync_sub_and_fetch(&m_value,static_cast<IntType>(n));
			return *this;
		}
		/// Prefix increment.
		atomic_counter_gcc_41 &operator++()
		{
			return operator+=(static_cast<IntType>(1));
		}
		/// Prefix decrement.
		atomic_counter_gcc_41 &operator--()
		{
			return operator-=(static_cast<IntType>(1));
		}
		/// Postfix increment.
		atomic_counter_gcc_41 operator++(int)
		{
			return atomic_counter_gcc_41(__sync_fetch_and_add(&m_value,static_cast<IntType>(1)));
		}
		/// Postfix decrement.
		atomic_counter_gcc_41 operator--(int)
		{
			return atomic_counter_gcc_41(__sync_fetch_and_sub(&m_value,static_cast<IntType>(1)));
		}
		/// Get copy of internal value.
		IntType get_value() const
		{
			return __sync_fetch_and_add(&m_value,static_cast<IntType>(0));
		}
		/// Fast type-trait for increment operation.
		static const bool is_increment_fast = true;
		/// Fast type-trait for arithmetics.
		static const bool is_arithmetics_fast = true;
	private:
		/// Internal value.
		/**
		 * Declared mutable because atomic_counter_gcc_41::get_value needs to perform the operation
		 * this + 0 in order to fetch safely the current m_value with GCC's atomic builtins.
		 */
		mutable IntType m_value;
};

#endif
 