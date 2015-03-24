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

#ifndef PAGMO_MIGRATION_BASE_S_POLICY_H
#define PAGMO_MIGRATION_BASE_S_POLICY_H

#include <boost/shared_ptr.hpp>
#include <vector>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"

namespace pagmo { namespace migration {

// Base class for selection policies for migration.
class base_s_policy;

/// Shared pointer to base selection policy.
typedef boost::shared_ptr<base_s_policy> base_s_policy_ptr;

/// Base class for migration selection policies.
/**
 * The task of a migration selection policy is to select in a population the individuals that will emigrate. The selection
 * is performed by the pure virtual select() method.
 *
 * The base::get_n_individuals() method for this class is meant to represent the number of individuals emigrating from the population.
 *
 * @author Marek Ruciński (marek.rucinski@gmail.com)
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE base_s_policy: public base
{
	public:
		base_s_policy(const double &rate = 1, rate_type type = absolute);
		virtual ~base_s_policy();
		/// Clone method.
		/**
		 * Provided that the derived policy implements properly the copy constructor, virtually all implementations of this method will
		 * look like this:
@verbatim
return base_ptr(new derived_policy(*this));
@endverbatim
		 *
		 * @return migration::base_s_policy_ptr to a copy of this.
		 */
		virtual base_s_policy_ptr clone() const = 0;
		/// Select individuals to emigrate from the given population.
		/**
		 * This is the method that actually implements the policy.
		 * Output vector should contain copies of selected individuals.
		 *
		 * \param[in,out] pop source population. In some cases (see best_kill_s-policy)
		 * it can be modified (e.g. killing the selected individual)
		 *
		 * \return a vector containing selected individuals.
		 */
		virtual std::vector<population::individual_type> select(population &pop) const = 0;
	private:	
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
		}
};

}}

BOOST_SERIALIZATION_ASSUME_ABSTRACT(pagmo::migration::base_s_policy)

#endif
