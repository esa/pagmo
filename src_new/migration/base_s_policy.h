/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
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
#include <iostream>
#include <string>
#include <vector>

#include "../config.h"
#include "../population.h"

namespace pagmo {

/// Migration policies namespace.
/**
 * This namespace contains selection/replacement policies used during migration in the archipelago class.
 */
namespace migration {

// Base class for selection policies for migration.
class base_s_policy;

/// Shared pointer to base selection policy.
typedef boost::shared_ptr<base_s_policy> base_s_policy_ptr;

/// Base class for selection policies for migration.
/**
 * This class provides its subclasses with means for specifying the migration rate.
 * The migration rate can be specified either as an absolute value (the number of individuals to migrate)
 * or as a fracion of the population size. The type of migration rate is selected upon construction.
 *
 * @author Marek Ruciński (marek.rucinski@gmail.com)
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE base_s_policy
{
	public:
		/// Type of migration rate.
		enum migration_rate_type {
			/// Migration rate is interpreted as the absolute number of individuals to migrate.
			absolute = 0,
			/// Migration rate is interpreted as the fraction of individuals to migrate.
			fractional = 1
		};
		base_s_policy(const double &rate,  migration_rate_type = fractional);
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
		 * \param[in] pop source population.
		 *
		 * \return a vector containing selected individuals.
		 */
		virtual std::vector<population::individual_type> select(const population &pop) const = 0;
		population::size_type get_n_individuals(const population &) const;
		std::string human_readable() const;
	protected:
		virtual std::string human_readable_extra() const;
	protected:
		/// Migration rate.
		/**
		 * It will be interpreted as an integer in case of absolute rate migration type, as a floating-point value
		 * in case of fractional migration type.
		 */
		double			m_rate;
		/// Migration rate type.
		migration_rate_type	m_rate_type;
};

std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const base_s_policy &);

}}

#endif
