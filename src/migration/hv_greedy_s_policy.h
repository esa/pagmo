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

#ifndef PAGMO_MIGRATION_HV_GREEDY_S_POLICY_H
#define PAGMO_MIGRATION_HV_GREEDY_S_POLICY_H

#include <vector>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"
#include "base_s_policy.h"

namespace pagmo { namespace migration {

///  Choose 'n' successive greatest contributors migration policy
/**
 * This policy revolves around choosing the individuals that contribute the greatest amount of volume to the total hypervolume.
 * Individuals are chosen iteratively, thus it is regarded as a greedy strategy.
 *
 * Note: If the problem is a single objective one, we fall back to the pagmo::migration::best_s_policy instead.
 *
 * @see Karl Bringmann, Tobias Friedrich, "An efficient algorithm for computing hypervolume contributions", Evolutionary Computation, v.18 n.3, p.383-402, Fall 2010
 * @see Karl Bringmann, Tobias Friedrich, "Don't be greedy when calculating hypervolume contributions",
 * Proceedings of the 10th International Workshop on Foundations of Genetic Algorithms (FOGA), pp. 103-112, ACM Press, 2009.
 *
 * @author Marcus Maertens (mmarcusx@gmail.com)
 * @author Krzysztof Nowak (kn@kiryx.net)
 */
class __PAGMO_VISIBLE hv_greedy_s_policy: public base_s_policy
{
public:
	hv_greedy_s_policy(const double &rate = 1, rate_type type = absolute, double nadir_eps = 1.0);
	base_s_policy_ptr clone() const;
	std::vector<population::individual_type> select(population &) const;
private:
	const double m_nadir_eps;

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base_s_policy>(*this);
		ar & const_cast<double &>(m_nadir_eps);
	}
};

} }

BOOST_CLASS_EXPORT_KEY(pagmo::migration::hv_greedy_s_policy)

#endif
