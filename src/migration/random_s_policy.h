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

#ifndef PAGMO_MIGRATION_RANDOM_S_POLICY_H
#define PAGMO_MIGRATION_RANDOM_S_POLICY_H

#include <vector>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"
#include "base_s_policy.h"

namespace pagmo { namespace migration {

/// Random migration selection policy.
/**
 * This policy chooses the the individuals according to uniform random distribution
 *
 * @author Marcus Maertens (mmarcusx@gmail.com)
 * 
 */
class __PAGMO_VISIBLE random_s_policy: public base_s_policy
{
	public:
		random_s_policy(const double &rate = 1, rate_type type = absolute);
		base_s_policy_ptr clone() const;
		std::vector<population::individual_type> select(population &) const;
	private:
		struct dom_comp;
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base_s_policy>(*this);
		}
};

} }

BOOST_CLASS_EXPORT_KEY(pagmo::migration::random_s_policy)

#endif
