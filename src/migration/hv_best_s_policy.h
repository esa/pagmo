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

#ifndef PAGMO_MIGRATION_HV_BEST_S_POLICY_H
#define PAGMO_MIGRATION_HV_BEST_S_POLICY_H

#include <vector>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"
#include "base_s_policy.h"

namespace pagmo { namespace migration {

///  Choose 'n' greatest contributors migration policy
/**
 * This policy chooses the individuals that contribute the greatest amount of volume to the total hypervolume.
 * Contributions of all individuals are computed, after which a set of 'n' best candidates is determined.
 *
 * Note: If the problem is a single objective one, we fall back to the pagmo::migration::best_s_policy instead.
 *
 * @author Marcus Maertens (mmarcusx@gmail.com)
 * @author Krzysztof Nowak (kn@kiryx.net)
 */
class __PAGMO_VISIBLE hv_best_s_policy: public base_s_policy
{
public:
	hv_best_s_policy(const double &rate = 1, rate_type type = absolute, double nadir_eps = 1.0);
	base_s_policy_ptr clone() const;
	std::vector<population::individual_type> select(population &) const;
private:
	const double m_nadir_eps;
	static bool sort_point_pairs_desc(const std::pair<unsigned int, double> &, const std::pair<unsigned int, double> &);

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base_s_policy>(*this);
		ar & const_cast<double &>(m_nadir_eps);
	}
};

} }

BOOST_CLASS_EXPORT_KEY(pagmo::migration::hv_best_s_policy)

#endif
