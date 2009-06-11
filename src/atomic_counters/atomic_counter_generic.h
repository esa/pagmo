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

#ifndef ATOMIC_COUNTER_GENERIC_H
#define ATOMIC_COUNTER_GENERIC_H

#include <boost/numeric/conversion/cast.hpp>
#include <boost/static_assert.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_same.hpp>

/// Generic atomic counter class.
/**
 * Will use a mutex internally to provide atomic operations on a generic integer type IntType.
 * Performance will be poor due to the locking mechanism.
 *
 * If IntType is not an integral type, a compile-time error will be produced upon instantiation.
 */
template <class IntType>
class atomic_counter_generic
{
		BOOST_STATIC_ASSERT(boost::is_integral<IntType>::value);
		typedef boost::lock_guard<boost::mutex> lock_type;
	public:
		/// Default constructor.
		/**
		 * It will initialise the internal value to zero.
		 */
		atomic_counter_generic():m_value(0),m_mutex() {}
		/// Constructor from generic integer type.
		/**
		 * IntType must be constructible from IntType2.
		 * @throws boost::numeric::bad_numeric_cast if IntType and IntType2 are not the same type
		 * and conversion from IntType2 to IntType results in loss of information.
		 */
		template <class IntType2>
		atomic_counter_generic(const IntType2 &n):m_value(n),m_mutex()
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
		atomic_counter_generic(const atomic_counter_generic &a):m_value(a.get_value()),m_mutex() {}
		/// Assignment operator.
		/**
		 * This operation acquires atomically the internal value of a, but the assignment
		 * to this is not atomic.
		 */
		atomic_counter_generic &operator=(const atomic_counter_generic &a)
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
		atomic_counter_generic &operator+=(const IntType2 &n)
		{
			BOOST_STATIC_ASSERT(boost::is_integral<IntType2>::value);
			lock_type lock(m_mutex);
			m_value += n;
			return *this;
		}
		/// In-place subtraction.
		/**
		 * The operation will be forwarded to the internal value.
		 */
		template <class IntType2>
		atomic_counter_generic &operator-=(const IntType2 &n)
		{
			BOOST_STATIC_ASSERT(boost::is_integral<IntType2>::value);
			lock_type lock(m_mutex);
			m_value -= n;
			return *this;
		}
		/// Prefix increment.
		atomic_counter_generic &operator++()
		{
			return operator+=(static_cast<IntType>(1));
		}
		/// Prefix decrement.
		atomic_counter_generic &operator--()
		{
			return operator-=(static_cast<IntType>(1));
		}
		/// Postfix increment.
		atomic_counter_generic operator++(int)
		{
			lock_type lock(m_mutex);
			atomic_counter_generic retval(m_value);
			++m_value;
			return retval;
		}
		/// Postfix decrement.
		atomic_counter_generic operator--(int)
		{
			lock_type lock(m_mutex);
			atomic_counter_generic retval(m_value);
			--m_value;
			return retval;
		}
		/// Get copy of internal value.
		IntType get_value() const
		{
			lock_type lock(m_mutex);
			return m_value;
		}
		/// Fast type-trait for increment operation.
		static const bool is_increment_fast = false;
		/// Fast type-trait for arithmetics.
		static const bool is_arithmetics_fast = false;
	private:
		/// Internal value.
		IntType					m_value;
		/// Mutex.
		mutable boost::mutex	m_mutex;
};

#endif
