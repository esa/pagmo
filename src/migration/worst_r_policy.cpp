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
#include "base.h"
#include "base_r_policy.h"
#include "worst_r_policy.h"

namespace pagmo { namespace migration {

/// Constructor from rate and rate type.
/**
 * @param[in] rate migration rate.
 * @param[in] type migration rate type.
 *
 * @see base_r_policy::base_r_policy.
 */
worst_r_policy::worst_r_policy(const double &rate, rate_type type):base_r_policy(rate,type) {}

base_r_policy_ptr worst_r_policy::clone() const
{
	return base_r_policy_ptr(new worst_r_policy(*this));
}

// Helper object used to sort arrays of indices of object placed in another container.
template <class Container>
struct indirect_individual_sorter
{
	indirect_individual_sorter(const Container &container, const population &pop):
		m_container(container),m_pop(pop) {}
	template <class Idx>
	bool operator()(const Idx &idx1, const Idx &idx2) const
	{
		typedef typename Container::const_iterator::difference_type diff_type;
		return m_pop.n_dominated(*(m_container.begin() + boost::numeric_cast<diff_type>(idx1))) >
			m_pop.n_dominated(*(m_container.begin() + boost::numeric_cast<diff_type>(idx2)));
	}
	// The original container.
	const Container		&m_container;
	const population	&m_pop;
};

// Selection implementation.
std::vector<std::pair<population::size_type,std::vector<population::individual_type>::size_type> >
	worst_r_policy::select(const std::vector<population::individual_type> &immigrants, const population &dest) const
{
	const population::size_type rate_limit = std::min<population::size_type>(get_n_individuals(dest),boost::numeric_cast<population::size_type>(immigrants.size()));
	// Temporary vectors to store sorted indices of the populations.
	std::vector<population::size_type> immigrants_idx(boost::numeric_cast<std::vector<population::size_type>::size_type>(immigrants.size()));
	std::vector<population::size_type> dest_idx(boost::numeric_cast<std::vector<population::size_type>::size_type>(dest.size()));
	// Fill in the arrays of indices.
	iota(immigrants_idx.begin(),immigrants_idx.end(),population::size_type(0));
	iota(dest_idx.begin(),dest_idx.end(),population::size_type(0));
	// Sort the arrays of indices.
	// From best to worst.
	std::sort(immigrants_idx.begin(),immigrants_idx.end(),indirect_individual_sorter<std::vector<population::individual_type> >(immigrants,dest));
	// From worst to best.
	std::sort(dest_idx.begin(),dest_idx.end(),indirect_individual_sorter<population>(dest,dest));
	std::reverse(dest_idx.begin(),dest_idx.end());
	// Create the result.
	std::vector<std::pair<population::size_type,std::vector<population::individual_type>::size_type> > result;
	for (population::size_type i = 0; i < rate_limit; ++i) {
		// Similar to fair policy, but replace unconditionally, without checking if the incoming individuals are better.
		result.push_back(std::make_pair(dest_idx[i],immigrants_idx[i]));
	}
	return result;
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::migration::worst_r_policy)
