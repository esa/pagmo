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

#ifndef PAGMO_MIGRATION_REPLACEMENT_POLICY_H
#define PAGMO_MIGRATION_REPLACEMENT_POLICY_H

#include <boost/shared_ptr.hpp>
#include <list>
#include <utility>

#include "../../../config.h"
#include "../population.h"

namespace pagmo
{

/// Base class for replacement policies for migration.
/**
 * The class provides it's subclasses with the basic means for determining the maximum incoming migration rate
 * (i.e. the maximum number of individuals which are replaced each time migrating individuals arrive)
 * in both absolute and fractional way.
 */
class __PAGMO_VISIBLE MigrationReplacementPolicy
{
	public:
		/// Default constructor.
		/** Maximum migration rate is assumed to be 1.0 (whole destination population can be replaced). */
		MigrationReplacementPolicy():maxMigrationRateAbs(-1), maxMigrationRateFrac(1.0) { }

		/// Constructor.
		/** Creates a policy with specified absolute maximum migration rate.
		 * \param[in] _maxMigrationRateAbs desired absolute migration rate.
		 */
		MigrationReplacementPolicy(const int& _maxMigrationRateAbs):maxMigrationRateAbs(_maxMigrationRateAbs), maxMigrationRateFrac(0.0) { }

		/// Constructor.
		/** Creates a policy with specified fractional migration rate.
		 * \param[in] _maxMigrationRateFrac desired fractional migration rate.
		 */
		MigrationReplacementPolicy(const double& _maxMigrationRateFrac):maxMigrationRateAbs(-1), maxMigrationRateFrac(_maxMigrationRateFrac) { }

		/// Copy constructor.
		MigrationReplacementPolicy(const MigrationReplacementPolicy& mrp):maxMigrationRateAbs(mrp.maxMigrationRateAbs), maxMigrationRateFrac(mrp.maxMigrationRateFrac) { }

		/// Assign pairs of individuals for replacement during migration.
		/**
		 * Note, that this method does not alter the target population, it just provides the replacement choice.
		 * The actual replacement is done in another place.
		 * The first element of a pair should be the index of the one from the destination population.
		 * The second one - the index of the individual from the incoming population which is to replace the first one.
		 * \param[in] incomingPopulation A vector of incoming individuals.
		 * \param[in] destinationPopulation Destination population.
		 * \return Replacement assignment.
		 */
		virtual std::list<std::pair<int, int> > selectForReplacement(const std::vector<individual>& incomingPopulation, const population& destinationPopulation) = 0;

		/// Clone object.
		/**
		 * Cloned object should be the exact copy of the original, but should be safe to use in a multithreaded environment.
		 * \todo Determine if the state variables should also be cloned and a separate method for resetting them be provided, or not.
		 */
		virtual MigrationReplacementPolicy* clone() const = 0;

		/// Calculates the maximum number of individuals that could be replaced in the specified population according to the policy.
		/**
		 * \param[in] population destination population.
		 */
		int getMaxMigrationRate(const population& population);

		/// Set the absolute maximum migration rate.
		/** Fractional maximum migration rate parameter is reset to 0.0 */
		void setMaxMigrationRate(const int& _maxMigrationRateAbs) {
			maxMigrationRateAbs = _maxMigrationRateAbs;
			maxMigrationRateFrac = 0.0;
		}

		/// Set the fractional maximum migration rate.
		/** Absolute maximum migration rate parameter is reset to -1 (the fractional one will be used). */
		void setMigrationRate(const double & _maxMigrationRateFrac) {
			maxMigrationRateAbs = -1;
			maxMigrationRateFrac = _maxMigrationRateFrac;
		}

		virtual ~MigrationReplacementPolicy() {}
	protected:

		int maxMigrationRateAbs; ///< Maximum incoming migration rate (absolute value), -1 means: use fraction.

		double maxMigrationRateFrac; ///< Maximum incoming migration rate (fractional value).

	private:
		/// Stream output operator.
		friend __PAGMO_VISIBLE_FUNC std::ostream &operator<<(std::ostream &s, const MigrationReplacementPolicy& ms);
};

/// Stream output operator.
__PAGMO_VISIBLE_FUNC std::ostream &operator<<(std::ostream &s, const MigrationReplacementPolicy& ms);

}

#endif
