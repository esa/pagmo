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

#ifndef PAGMO_MIGRATION_POLICY_H
#define PAGMO_MIGRATION_POLICY_H

#include <boost/scoped_ptr.hpp>
#include <iostream>

#include "../../../config.h"
#include "ChooseBestMigrationSelectionPolicy.h"
#include "RandomMigrationReplacementPolicy.h"

namespace pagmo
{

/// This object gathers island-specific migration-related parameters.
/**
 * Introduction of this object has been lobbied by Dario because in his oppinion it's presence will increase
 * the user-friendliness of the python interface.
 */
class __PAGMO_VISIBLE MigrationPolicy
{
	public:
		/// Default constructor. All class fields are initialised to the default values imposed by Dario.
		/**
		 * These default values are:
		 * <ul>
		 *   <li>Migration probability: 1</li>
		 *   <li>Migration selection policy: ChooseBestMigrationSelectionPolicy with the migration rate of 1 individual</li>
		 *   <li>Migration selection policy: RandomMigrationReplacementPolicy with the maximum migration rate limit of 1.0 (the whole population)</li>
		 * </ul>
		 */
		MigrationPolicy()
				:migrationProbability(1.0),
				selectionPolicy(ChooseBestMigrationSelectionPolicy(1).clone()),
				replacementPolicy(RandomMigrationReplacementPolicy(1.0).clone()) { }

		/// Constructor with migration probability.
		/** This constructor allows specifying migration probability without touching the policies. For default values of policies see @see MigrationPolicy::MigrationPolicy(). */
		MigrationPolicy(const double _migrationProbability)
				:selectionPolicy(ChooseBestMigrationSelectionPolicy(1).clone()),
				replacementPolicy(RandomMigrationReplacementPolicy(1.0).clone()) {
			setMigrationProbability(_migrationProbability); // does range checking
		}

		/// Constructor.
		/** This constructor allows specyfiyng values of all migration parameters. */
		MigrationPolicy(const double _migrationProbability, const MigrationSelectionPolicy& _selectionPolicy, const MigrationReplacementPolicy& _replacementPolicy)
				:selectionPolicy(_selectionPolicy.clone()),
				replacementPolicy(_replacementPolicy.clone()) {
			setMigrationProbability(_migrationProbability); // does range checking
		}

		/// Copy constructor.
		/** Creates a deep copy of the given migration policy. */
		MigrationPolicy(const MigrationPolicy& mp)
				:migrationProbability(mp.migrationProbability),
				selectionPolicy(mp.selectionPolicy->clone()),
				replacementPolicy(mp.replacementPolicy->clone()) { }

		// Getters and setters.
		double getMigrationProbability() const {
			return migrationProbability;
		}

		void setMigrationProbability(const double _migrationProbability);


		MigrationSelectionPolicy& getMigrationSelectionPolicy() const {
			return *selectionPolicy;
		}

		void setMigrationSelectionPolicy(const MigrationSelectionPolicy& msp) {
			selectionPolicy.reset(msp.clone());
		}


		MigrationReplacementPolicy& getMigrationReplacementPolicy() const {
			return *replacementPolicy;
		}

		void setMigrationReplacementPolicy(const MigrationReplacementPolicy& mrp) {
			replacementPolicy.reset(mrp.clone());
		}


	private:

		/// Migration probability.
		/** This number specifies the probability of <b>inserting migrated individuals to the island's population</b>. */
		double migrationProbability;

		/// The selection policy. Must not be null.
		boost::scoped_ptr<MigrationSelectionPolicy> selectionPolicy;

		/// The replacement policy. Must not be null.
		boost::scoped_ptr<MigrationReplacementPolicy> replacementPolicy;

		/// Stream output operator.
		friend __PAGMO_VISIBLE_FUNC std::ostream &operator<<(std::ostream &s, const MigrationPolicy& msp);
};

/// The stream output operator.
__PAGMO_VISIBLE_FUNC std::ostream &operator<<(std::ostream &s, const MigrationPolicy& msp);

}

#endif
