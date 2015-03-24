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

#ifndef PAGMO_MIGRATION_BASE_R_POLICY_H
#define PAGMO_MIGRATION_BASE_R_POLICY_H

#include <boost/shared_ptr.hpp>
#include <utility>
#include <vector>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"

namespace pagmo { namespace migration {

// Fwd declaration.
class base_r_policy;

/// Shared pointer to base replacement policy.
typedef boost::shared_ptr<base_r_policy> base_r_policy_ptr;

/// Base class for migration replacement policies.
/**
 * The task of a migration replacement policy is to select in a population the individuals that will replaced by immigrating individuals. The selection
 * is performed by the pure virtual select() method.
 *
 * The base::get_n_individuals() method for this class is meant to represent the maximum number of individuals in the target population
 * that can be replaced by the immigrants. How many of these will actually be replaced will depend on the specific policy implementation.
 *
 * @author Marek Rucinski (marek.rucinski@gmail.com)
 * @author Francesco Biscani (bluescarni@gmail.com)
 * @author Dario Izzo (dario.izzo@gmail.com)
 */
class __PAGMO_VISIBLE base_r_policy: public base
{
	public:
		base_r_policy(const double &rate = 1, rate_type type = fractional);
		virtual ~base_r_policy();
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
		virtual base_r_policy_ptr clone() const = 0;
		/// Assign pairs of individuals for replacement during migration.
		/**
		 * Note, that this method does not alter the target population, it just provides the replacement choice.
		 * The actual replacement is done in the archipelago class.
		 * The first element of a pair should be the index of the one from the destination population.
		 * The second one - the index of the individual from the immigrants vector which is to replace the first one.
		 *
		 * \param[in] immigrants vector of incoming individuals.
		 * \param[in] destination population into which the immigrants will be replaced.
		 *
		 * \return replacement assignment.
		 */
		virtual std::vector<std::pair<population::size_type,std::vector<population::individual_type>::size_type> >
			select(const std::vector<population::individual_type> &immigrants, const population &destination) const = 0;
	protected:

	private:	
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
		}
};

}}

BOOST_SERIALIZATION_ASSUME_ABSTRACT(pagmo::migration::base_r_policy)

#endif
