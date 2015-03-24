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

#include <algorithm>
#include <boost/numeric/conversion/cast.hpp>
#include <utility>
#include <vector>

#include "../population.h"
#include "../rng.h"
#include "base.h"
#include "base_r_policy.h"
#include "random_r_policy.h"

namespace pagmo { namespace migration {

/// Constructor.
/**
 * @param[in] rate migration rate.
 * @param[in] type migration rate type.
 *
 * @see base_r_policy::base_r_policy.
 */
random_r_policy::random_r_policy(const double &rate, rate_type type):base_r_policy(rate,type),m_urng(rng_generator::get<rng_uint32>()) {}

base_r_policy_ptr random_r_policy::clone() const
{
	return base_r_policy_ptr(new random_r_policy(*this));
}

std::vector<std::pair<population::size_type,std::vector<population::individual_type>::size_type> >
	random_r_policy::select(const std::vector<population::individual_type> &immigrants, const population &dest) const
{
	const population::size_type rate_limit = std::min<population::size_type>(get_n_individuals(dest),boost::numeric_cast<population::size_type>(immigrants.size()));
	// Temporary vectors to store sorted indices of the populations.
	std::vector<population::size_type> immigrants_idx(boost::numeric_cast<std::vector<population::size_type>::size_type>(immigrants.size()));
	std::vector<population::size_type> dest_idx(boost::numeric_cast<std::vector<population::size_type>::size_type>(dest.size()));
	// Fill in the arrays of indices.
	iota(immigrants_idx.begin(),immigrants_idx.end(),population::size_type(0));
	iota(dest_idx.begin(),dest_idx.end(),population::size_type(0));
	// Permute the indices (immigrants).
	for (population::size_type i = 0; i < rate_limit; ++i) {
		population::size_type next_idx = i + (m_urng() % (rate_limit - i));
		if (next_idx != i) {
			std::swap(immigrants_idx[i], immigrants_idx[next_idx]);
		}
	}
	// Permute the indices (destination).
	for (population::size_type i = 0; i < rate_limit; ++i) {
		population::size_type next_idx = i + (m_urng() % (dest.size() - i));
		if (next_idx != i) {
			std::swap(dest_idx[i], dest_idx[next_idx]);
		}
	}
	// Return value.
	std::vector<std::pair<population::size_type,std::vector<population::individual_type>::size_type> >
		retval;
	for (population::size_type i = 0; i < rate_limit; ++i) {
		retval.push_back(std::make_pair(dest_idx[i],immigrants_idx[i]));
	}
	return retval;
}

} }

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::migration::random_r_policy)
