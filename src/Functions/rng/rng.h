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

// 27/12/08 Created by Francesco Biscani.

#ifndef PAGMO_RNG_H
#define PAGMO_RNG_H

#include <boost/cstdint.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>

// This rng returns an unsigned integer in the [0,2**32-1] range.
typedef boost::mt19937 rng_uint32;
// This rng returns a double in the [0,1[ range.
typedef boost::lagged_fibonacci607 rng_double;

// Generic thread-safe wrapper class around a Boost-like random number generator.
template <class Rng>
class static_rng {
	public:
		typedef typename Rng::result_type result_type;
		result_type operator()() {
			boost::lock_guard<boost::mutex> lock(m_mutex);
			return m_rng();
		}
	private:
		static boost::mutex	m_mutex;
		static Rng		m_rng;
};

template <class Rng>
boost::mutex static_rng<Rng>::m_mutex;

// Use as initial seed the number of microseconds elapsed since 01/01/1970, cast to
// uint32_t.
template <class Rng>
Rng static_rng<Rng>::m_rng(boost::uint32_t((boost::posix_time::microsec_clock::local_time() -
	boost::posix_time::ptime(boost::gregorian::date(1970,1,1))).total_microseconds()));

// Thread-safe uint32 rng.
typedef static_rng<rng_uint32> static_rng_uint32;

// Thread-safe double rng.
typedef static_rng<rng_double> static_rng_double;

#endif
