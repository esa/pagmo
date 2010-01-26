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

#ifndef PAGMO_RANDOM_MIGRATION_REPLACEMENT_POLICY_H
#define PAGMO_RANDOM_MIGRATION_REPLACEMENT_POLICY_H

#include <boost/cstdint.hpp>

#include "../../../config.h"
#include "../../../rng.h"
#include "MigrationReplacementPolicy.h"

namespace pagmo
{

/// Random migration replacement policy.
/**
 * Here randomly selected individuals in the destination population are replaced by the randomly selected
 * incoming individuals. The number of individuals replaced is limited by the maximum incoming rate parameter.
 */
class __PAGMO_VISIBLE RandomMigrationReplacementPolicy: public MigrationReplacementPolicy
{
	public:
		/// Constructor.
		/**
		 * \see MigrationReplacementPolicy::MigrationReplacementPolicy().
		 * \param[in] seed initial random seed for the internal rng. If it is not specified, it is generated using the static RNG.
		 */
		RandomMigrationReplacementPolicy(const boost::uint32_t seed = static_rng_uint32()()):MigrationReplacementPolicy(), rng(seed) { }

		/// Constructor.
		/**
		 * Allows maximum migration rate specification (absolute).
		 * \see MigrationReplacementPolicy::MigrationReplacementPolicy(const int&)
		 * \param[in] _maxMigrationRate desired maximum migration rate (absolute).
		 * \param[in] seed initial random seed for the internal rng. If it is not specified, it is generated using the static RNG.
		 */
		RandomMigrationReplacementPolicy(const int& _maxMigrationRate, const boost::uint32_t seed = static_rng_uint32()()):MigrationReplacementPolicy(_maxMigrationRate), rng(seed) { }

		/// Constructor.
		/**
		 * Allows maximum migration rate specification (fractional).
		 * \see MigrationReplacementPolicy::MigrationReplacementPolicy(const double&)
		 * \param[in] _maxMigrationRate desired migration rate (fractional).
		 * \param[in] seed initial random seed for the internal rng. If it is not specified, it is generated using the static RNG.
		 */
		RandomMigrationReplacementPolicy(const double& _maxMigrationRate, const boost::uint32_t seed = static_rng_uint32()()):MigrationReplacementPolicy(_maxMigrationRate), rng(seed) { }

		/// Copy constructor.
		/**
		 * The copy inherits migration rate parameters, but it is not strictly identical. This is because the RNG
		 * of the copy is seeded with the rng of the source, and thus future behaviour of the original and the copy
		 * will be different.
		 * \todo figure out if it is possible to create a truly identical copy.
		 */
		RandomMigrationReplacementPolicy(const RandomMigrationReplacementPolicy& rmrp):MigrationReplacementPolicy(rmrp),rng(rmrp.rng()) { }

		/// \see MigrationReplacementPolicy::selectForReplacement
		virtual std::list<std::pair<int, int> > selectForReplacement(const std::vector<individual>& incomingPopulation, const population& destinationPopulation);

		/// \see MigrationReplacementPolicy::clone
		virtual RandomMigrationReplacementPolicy* clone() const {
			return new RandomMigrationReplacementPolicy(*this);
		};

	private:
		/// Object used to generate sequences of number from 0 to n using the std::generate function.
		struct IndexGenerator {
			int currentValue; ///< Current counter value
			/// Counter constructor; Initialises the counter to 0.
			IndexGenerator() {
				currentValue = 0;
			}
			/// Call operator - give the next index value.
			int operator()() {
				return currentValue++;
			}
		};

		mutable rng_uint32 rng; ///< Random Number Generator
};

}

#endif
