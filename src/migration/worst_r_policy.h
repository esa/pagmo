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

#ifndef PAGMO_MIGRATION_WORST_R_POLICY_H
#define PAGMO_MIGRATION_WORST_R_POLICY_H

#include <utility>
#include <vector>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"
#include "base_r_policy.h"

namespace pagmo { namespace migration {

/// Worst replacement policy.
/**
 * Best individuals from the incoming population will forcefully replace the worst individuals in the
 * target population.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE worst_r_policy: public base_r_policy
{
	public:
		worst_r_policy(const double &rate = 1, rate_type type = absolute);
		base_r_policy_ptr clone() const;
		std::vector<std::pair<population::size_type,std::vector<population::individual_type>::size_type> >
			select(const std::vector<population::individual_type> &, const population &) const;
	private:	
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base_r_policy>(*this);
		}
};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::migration::worst_r_policy)

#endif
