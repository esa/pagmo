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

// 27/12/08 Created by Francesco Biscani.

#ifndef PAGMO_RNG_H
#define PAGMO_RNG_H

#include <boost/cstdint.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>

namespace pagmo
{
/// This rng returns an unsigned integer in the [0,2**32-1] range.
/**
 * @see http://www.boost.org/doc/libs/release/libs/random/random-generators.html
 */
typedef boost::mt19937 rng_uint32;
/// This rng returns a double in the [0,1[ range.
/**
 * @see http://www.boost.org/doc/libs/release/libs/random/random-generators.html
 */
typedef boost::lagged_fibonacci607 rng_double;

/// Generic thread-safe wrapper class around a Boost-like pseudo-random number generator.
/**
 * To use, construct and call operator() to get a pseudo-random number. Type Rng must be a Boost-like
 * pseudo-random number generator.
 *
 * Implementation internally uses a mutex, so that this generator can
 * be safely called concurrently from multiple threads. The initial seed used
 * is the number of microseconds elapsed since 01/01/1970, cast to uint32_t.
 * Please note that the initial seed is set once at program startup and shared among all
 * instances for a given Rng type.
 * @see http://www.boost.org/doc/libs/release/libs/random/index.html
 */
template <class Rng>
class static_rng {
	public:
		/// Result type.
		typedef typename Rng::result_type result_type;
		/// Return random number.
		/**
		 * Return next pseudo-random number in the sequence.
		 */
		result_type operator()()
		{
			boost::lock_guard<boost::mutex> lock(m_mutex);
			return m_rng();
		}
		/// Set seed.
		/**
		 * Set seed to n. Thread-safe. It is assumed that the underlying Rng type
		 * can be successfully seeded with an int value.
		 */
		static void set_seed(int n)
		{
			boost::lock_guard<boost::mutex> lock(m_mutex);
			m_rng.seed(n);
		}
	private:
		static boost::mutex	m_mutex;
		static Rng		m_rng;
};

template <class Rng>
boost::mutex static_rng<Rng>::m_mutex;

// Use as initial seed the number of microseconds elapsed since 01/01/1970, cast to uint32_t.
template <class Rng>
Rng static_rng<Rng>::m_rng(boost::uint32_t((boost::posix_time::microsec_clock::local_time() -
	boost::posix_time::ptime(boost::gregorian::date(1970,1,1))).total_microseconds()));

/// Thread-safe version of pagmo::rng_uint32.
typedef static_rng<rng_uint32> static_rng_uint32;

/// Thread-safe version of pagmo::rng_double.
typedef static_rng<rng_double> static_rng_double;

}

#endif
