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

#ifndef PAGMO_CHOOSE_BEST_MIGRATION_SELECTION_POLICY_H
#define PAGMO_CHOOSE_BEST_MIGRATION_SELECTION_POLICY_H

#include "../../../config.h"
#include "MigrationSelectionPolicy.h"

namespace pagmo
{

/// Choose best migration selection policy.
/**
 * This policy is to choose best individuals from the population as migrating individuals.
 */
class __PAGMO_VISIBLE ChooseBestMigrationSelectionPolicy: public MigrationSelectionPolicy
{
	public:
		/// Constructor.
		/**
		 * \see MigrationSelectionPolicy::MigrationSelectionPolicy().
		 */
		ChooseBestMigrationSelectionPolicy():MigrationSelectionPolicy() { }

		/// Constructor.
		/** \see MigrationSelectionPolicy::MigrationSelectionPolicy(const int&) */
		ChooseBestMigrationSelectionPolicy(const int& _migrationRate):MigrationSelectionPolicy(_migrationRate) { }

		/// Constructor.
		/** \see MigrationSelectionPolicy::MigrationSelectionPolicy(const double&) */
		ChooseBestMigrationSelectionPolicy(const double& _migrationRate):MigrationSelectionPolicy(_migrationRate) { }

		/// Copy constructor.
		ChooseBestMigrationSelectionPolicy(const ChooseBestMigrationSelectionPolicy& rmsp):MigrationSelectionPolicy(rmsp) { }

		/// Virtual destructor.
		virtual ~ChooseBestMigrationSelectionPolicy() { }

		/// \see MigrationSelectionPolicy::selectForMigration
		virtual std::vector<individual> selectForMigration(const population& population);

		/// \see MigrationSelectionPolicy::clone
		virtual ChooseBestMigrationSelectionPolicy* clone() const {
			return new ChooseBestMigrationSelectionPolicy(*this);
		}
};

}

#endif
