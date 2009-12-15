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

// 09/03/2009: Initial version by Marek Rucinski.

#ifndef PAGMO_RANDOM_MIGRATION_SELECTION_POLICY_H
#define PAGMO_RANDOM_MIGRATION_SELECTION_POLICY_H

#include <boost/cstdint.hpp>

#include "../../../config.h"
#include "../../../Functions/rng/rng.h"
#include "MigrationSelectionPolicy.h"

// TODO: fix missing headers.

namespace pagmo
{

/// Random migration selection policy.
/**
 * This policy is to choose migrating individuals randomly.
 */
class __PAGMO_VISIBLE RandomMigrationSelectionPolicy: public MigrationSelectionPolicy
{
	public:
		/// Constructor.
		/**
		 * \see MigrationSelectionPolicy::MigrationSelectionPolicy().
		 * \param[in] seed initial random seed for the internal rng. If it is not specified, it is generated using the static RNG.
		 */
		RandomMigrationSelectionPolicy(const boost::uint32_t seed = static_rng_uint32()()):MigrationSelectionPolicy(), rng(seed) { }

		/// Constructor.
		/**
		 * Allows migration rate specification (absolute).
		 * \see MigrationSelectionPolicy::MigrationSelectionPolicy(const int&)
		 * \param[in] _migrationRate desired migration rate (absolute).
		 * \param[in] seed initial random seed for the internal rng. If it is not specified, it is generated using the static RNG.
		 */
		RandomMigrationSelectionPolicy(const int& _migrationRate, const boost::uint32_t seed = static_rng_uint32()()):MigrationSelectionPolicy(_migrationRate), rng(seed) { }

		/// Constructor.
		/**
		 * Allows migration rate specification (fractional).
		 * \see MigrationSelectionPolicy::MigrationSelectionPolicy(const double&)
		 * \param[in] _migrationRate desired migration rate (fractional).
		 * \param[in] seed initial random seed for the internal rng. If it is not specified, it is generated using the static RNG.
		 */
		RandomMigrationSelectionPolicy(const double& _migrationRate, const boost::uint32_t seed = static_rng_uint32()()):MigrationSelectionPolicy(_migrationRate), rng(seed) { }

		/// Copy constructor.
		/**
		 * The copy inherits migration rate parameters, but it is not strictly identical. This is because the RNG
		 * of the copy is seeded with the rng of the source, and thus future behaviour of the original and the copy
		 * will be different.
		 * \todo figure out if it is possible to create a truly identical copy.
		 */
		RandomMigrationSelectionPolicy(const RandomMigrationSelectionPolicy& rmsp):MigrationSelectionPolicy(rmsp),rng(rmsp.rng()) { }

		/// Virtual destructor.
		virtual ~RandomMigrationSelectionPolicy() { }

		/// \see MigrationSelectionPolicy::selectForMigration
		virtual std::vector<individual> selectForMigration(const population& population);

		/// \see MigrationSelectionPolicy::clone
		virtual RandomMigrationSelectionPolicy* clone() const {
			return new RandomMigrationSelectionPolicy(*this);
		}

	private:
		mutable rng_uint32 rng; ///< Random Number Generator
};

}

#endif
