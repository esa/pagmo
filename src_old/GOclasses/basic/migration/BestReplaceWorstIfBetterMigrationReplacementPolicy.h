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

#ifndef PAGMO_BEST_REPLACE_WORST_IF_BETTER_MIGRATION_REPLACEMENT_POLICY_H
#define PAGMO_BEST_REPLACE_WORST_IF_BETTER_MIGRATION_REPLACEMENT_POLICY_H

#include "../../../config.h"
#include "../../../rng.h"
#include "MigrationReplacementPolicy.h"

namespace pagmo
{

/// Best-replace-worst-if-better migration replacement policy.
/**
 * Best individuals from the incoming population are matched with the worst ones from the destination
 * population. The replacement is moreover condition to that the new individual is better than the old one.
 */
class __PAGMO_VISIBLE BestReplaceWorstIfBetterMigrationReplacementPolicy: public MigrationReplacementPolicy
{
	public:
		/// Constructor.
		/** \see MigrationReplacementPolicy::MigrationReplacementPolicy(). */
		BestReplaceWorstIfBetterMigrationReplacementPolicy():MigrationReplacementPolicy() { }

		/// Constructor.
		/** \see MigrationReplacementPolicy::MigrationReplacementPolicy(const int&) */
		BestReplaceWorstIfBetterMigrationReplacementPolicy(const int& _maxMigrationRate):MigrationReplacementPolicy(_maxMigrationRate) { }

		/// Constructor.
		/** \see MigrationReplacementPolicy::MigrationReplacementPolicy(const double&) */
		BestReplaceWorstIfBetterMigrationReplacementPolicy(const double& _maxMigrationRate):MigrationReplacementPolicy(_maxMigrationRate) { }

		/// Copy constructor.
		BestReplaceWorstIfBetterMigrationReplacementPolicy(const BestReplaceWorstIfBetterMigrationReplacementPolicy& rmrp):MigrationReplacementPolicy(rmrp) { }

		/// Virtual destructor.
		virtual ~BestReplaceWorstIfBetterMigrationReplacementPolicy() { }

		/// \see MigrationReplacementPolicy::selectForReplacement
		virtual std::list<std::pair<int, int> > selectForReplacement(const std::vector<individual>& incomingpopulation, const population& destinationpopulation);

		/// \see MigrationReplacementPolicy::clone
		virtual BestReplaceWorstIfBetterMigrationReplacementPolicy* clone() const {
			return new BestReplaceWorstIfBetterMigrationReplacementPolicy(*this);
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

		/// Helper object used to sort arrays of indices of object placed in another arrays
		struct IndirectIndividualSorter {
			const std::vector<individual>& actualArray;

			IndirectIndividualSorter(const std::vector<individual>& _actualArray)
					:actualArray(_actualArray) { }

			int operator() (int i, int j) {
				return individual::compare_by_fitness(actualArray[i], actualArray[j]);
			}
		};
};

}

#endif
