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

#ifndef PAGMO_MIGRATION_HV_FAIR_R_POLICY_H
#define PAGMO_MIGRATION_HV_FAIR_R_POLICY_H

#include <utility>
#include <vector>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"
#include "base_r_policy.h"

namespace pagmo { namespace migration {

/// Replace the 'n' least contributors with the incoming set.
/**
 * Policy replaces a set of 'n' least contributors with the incoming set of immigrants in a 'fair' way.
 * That is: an original inhabitant of the island is only substituted if its exclusive hypervolume contribution is less than the 
 * the worst contribution of a replaced immigrant.
 * Least contributors are computed at once for single population (per-front).
 *
 * Note: If the problem is a single objective one, we fall back to the pagmo::migration::fair_r_policy instead.
 *
 * @author Marcus Maertens (mmarcusx@gmail.com)
 * @author Krzysztof Nowak (kn@kiryx.net)
 */
class __PAGMO_VISIBLE hv_fair_r_policy: public base_r_policy
{
public:
	hv_fair_r_policy(const double &rate = 1, rate_type type = absolute, double nadir_eps = 1.0);
	base_r_policy_ptr clone() const;
	std::vector<std::pair<population::size_type,std::vector<population::individual_type>::size_type> >
		select(const std::vector<population::individual_type> &, const population &) const;
private:
	static bool ind_cmp(const std::pair<unsigned int, double> &, const std::pair<unsigned int, double> &);
	static bool sort_point_pairs_asc(const std::pair<unsigned int, double> &a, const std::pair<unsigned int, double> &b);
	const double m_nadir_eps;

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base_r_policy>(*this);
		ar & const_cast<double &>(m_nadir_eps);
	}
};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::migration::hv_fair_r_policy)

#endif
