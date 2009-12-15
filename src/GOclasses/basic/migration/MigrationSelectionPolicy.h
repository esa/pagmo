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

#ifndef PAGMO_MIGRATION_SELECTION_POLICY_H
#define PAGMO_MIGRATION_SELECTION_POLICY_H

#include <vector>

#include "../../../config.h"
#include "../individual.h"
#include "../population.h"

namespace pagmo
{

/// Base class for selection policies for migration.
/**
 * This class provides it's subclasses with means for specifying the migration rate.
 * The migration rate can be specified either as an absolute value (the number of individuals to migrate)
 * or as a fracion of the population size.
 * The meaningful value is determined in the following way: if the absolute rate is negative, fractional
 * value is assumed to be correct. Otherwise, the absolute value is taken.
 */
class __PAGMO_VISIBLE MigrationSelectionPolicy
{
	public:
		/// Default constructor.
		/** Migration rate is assumed to be 1 individual. */
		MigrationSelectionPolicy():migrationRateAbs(1), migrationRateFrac(0.0) { }

		/// Constructor.
		/** Creates a policy with specified absolute migration rate.
		 * \param[in] _migrationRateAbs desired absolute migration rate.
		 */
		MigrationSelectionPolicy(const int& _migrationRateAbs):migrationRateAbs(_migrationRateAbs), migrationRateFrac(0.0) { }

		/// Constructor.
		/** Creates a policy with specified fractional migration rate.
		 * \param[in] _migrationRateFrac desired fractional migration rate.
		 */
		MigrationSelectionPolicy(const double& _migrationRateFrac):migrationRateAbs(-1), migrationRateFrac(_migrationRateFrac) { }

		/// Copy constructor.
		MigrationSelectionPolicy(const MigrationSelectionPolicy& msp):migrationRateAbs(msp.migrationRateAbs), migrationRateFrac(msp.migrationRateFrac) { }

		/// Select individuals to migrate out of the given population.
		/**
		 * This is the method that actually implements the policy.
		 * Output vector should contain <b>copies</b> of selected individuals.
		 * \param[in] population Source population.
		 * \return A vector containing selected individuals.
		 */
		virtual std::vector<individual> selectForMigration(const population& pop) = 0;

		/// Clone object.
		/**
		 * Cloned object should be the exact copy of the original, but should be safe to use in a multithreaded environment.
		 * \todo Determine if the state variables should also be cloned and a separate method for resetting them be provided, or not.
		 */
		virtual MigrationSelectionPolicy* clone() const = 0;

		/// Get the migration rate for the given population.
		int getNumberOfIndividualsToMigrate(const population& population);

		/// Set the absolute migration rate.
		/** Fractional migration rate parameter is reset to 0.0 */
		void setMigrationRate(const int& _migrationRateAbs) {
			migrationRateAbs = _migrationRateAbs;
			migrationRateFrac = 0.0;
		}

		/// Set the fractional migration rate.
		/** Absolute migration rate parameter is reset to -1 (the fractional one will be used). */
		void setMigrationRate(const double & _migrationRateFrac) {
			migrationRateAbs = -1;
			migrationRateFrac = _migrationRateFrac;
		}

	protected:
		/// Migration Rate (absolute value)
		/** This variable specifies the number of individuals to migrate. */
		int migrationRateAbs;
		/// Migration rate (fractional value)
		/** This variable specifies the fraction of the population to migrate. */
		double migrationRateFrac;

	private:
		/// Stream output operator.
		friend __PAGMO_VISIBLE_FUNC std::ostream &operator<<(std::ostream &s, const MigrationSelectionPolicy& ms);
};

/// Stream output operator.
__PAGMO_VISIBLE_FUNC std::ostream &operator<<(std::ostream &s, const MigrationSelectionPolicy& ms);

}

#endif
