/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
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
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <sstream>
#include <string>

#include "serialization.h"
#include "config.h"

namespace pagmo
{
/// This rng returns an unsigned integer in the [0,2**32-1] range.
/**
 * @see http://www.boost.org/doc/libs/release/libs/random/random-generators.html
 */
class __PAGMO_VISIBLE rng_uint32: public boost::mt19937 {
		friend class boost::serialization::access;
	public:
		/// Return value of the generator.
		typedef boost::mt19937::result_type result_type;
		/// Default constructor.
		/**
		 * Will invoke the base default constructor.
		 */
		rng_uint32():boost::mt19937() {}
		/// Constructor from unsigned integer.
		/**
		 * Will invoke the corresponding base constructor.
		 */
		rng_uint32(const result_type &n):boost::mt19937(n) {}
		// Default generated copy ctor and assignment are fine.
	private:
		// Serialization exploits the fact that the state of Boost RNGs
		// can be sent/received to/from standard streams.
		template <class Archive>
		void save(Archive &ar, const unsigned int) const
		{
			std::stringstream ss;
			ss << *static_cast<boost::mt19937 const *>(this);
			std::string tmp(ss.str());
			ar << tmp;
		}
		template <class Archive>
		void load(Archive &ar, const unsigned int)
		{
			std::string tmp;
			ar >> tmp;
			std::stringstream ss(tmp);
			ss >> *static_cast<boost::mt19937 *>(this);
		}
		BOOST_SERIALIZATION_SPLIT_MEMBER()
};

/// This rng returns a double in the [0,1[ range.
/**
 * @see http://www.boost.org/doc/libs/release/libs/random/random-generators.html
 */
class __PAGMO_VISIBLE rng_double: public boost::lagged_fibonacci607 {
		friend class boost::serialization::access;
	public:
		/// Default constructor.
		/**
		 * Will invoke the base default constructor.
		 */
		rng_double():boost::lagged_fibonacci607() {}
		/// Constructor from unsigned integer.
		/**
		 * Will invoke the corresponding base constructor.
		 */
		rng_double(const boost::uint32_t &n):boost::lagged_fibonacci607(n) {}
		// Default generated copy ctor and assignment are fine.
	private:
		template <class Archive>
		void save(Archive &ar, const unsigned int) const
		{
			std::stringstream ss;
			ss << *static_cast<boost::lagged_fibonacci607 const *>(this);
			std::string tmp(ss.str());
			ar << tmp;
		}
		template <class Archive>
		void load(Archive &ar, const unsigned int)
		{
			std::string tmp;
			ar >> tmp;
			std::stringstream ss(tmp);
			ss >> *static_cast<boost::lagged_fibonacci607 *>(this);
		}
		BOOST_SERIALIZATION_SPLIT_MEMBER()
};

/// Generic thread-safe generator of pseudo-random number generators.
/**
 * To use, call the static member get() to get a pseudo-random number generator seeded with an initial pseudo-random value.
 *
 * Implementation internally uses a mutex, so that this generator can
 * be safely called concurrently from multiple threads. The initial seed used
 * is the number of microseconds elapsed since 01/01/1970, cast to uint32_t.
 *
 * @see http://www.boost.org/doc/libs/release/libs/random/index.html
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */



class __PAGMO_VISIBLE rng_generator {
	public:
		/// Return pseudo-random number generator.
		/**
		 * Type Rng must be a Boost-like pseudo-random number generator initialisable
		 * with a boost::uint32_t. Return value is seeded with an internal
		 * static pagmo::rng_uint32.
		 *
		 * @return pseudo-random number generator seeded with pseudo-random value.
		 */

		template <class Rng> 
		static Rng get();
		static void set_seed(int);

	private:
		static  boost::mutex  m_mutex;
		static  rng_uint32  m_seeder;
};

}

#endif
