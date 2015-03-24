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

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>

#include "rng.h"

namespace pagmo
{

boost::mutex rng_generator::m_mutex;

// Use as initial seed the number of microseconds elapsed since 01/01/1970, cast to uint32_t.
rng_uint32 rng_generator::m_seeder(boost::uint32_t((boost::posix_time::microsec_clock::local_time() -
	boost::posix_time::ptime(boost::gregorian::date(1970,1,1))).total_microseconds()));

/// Set seed.
/**
 * Set the seed of the internal generator to n. Thread-safe. Note that input integer n will be
 * cast to uint32_t.
 *
 * @param[in] n seed for the generator of pseudo-random number generators.
 */
void rng_generator::set_seed(int n)
{
	boost::lock_guard<boost::mutex> lock(m_mutex);
	m_seeder.seed(n);
}


template <class Rng> 
Rng rng_generator::get()
{
	boost::lock_guard<boost::mutex> lock(m_mutex);
	return Rng(m_seeder());
}

template __PAGMO_VISIBLE rng_double rng_generator::get<rng_double>();
template __PAGMO_VISIBLE rng_uint32 rng_generator::get<rng_uint32>();

}
