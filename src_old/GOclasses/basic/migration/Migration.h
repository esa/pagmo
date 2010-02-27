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

#ifndef PAGMO_MIGRATION_H
#define PAGMO_MIGRATION_H

#include <boost/scoped_ptr.hpp>
#include <iostream>

#include "../../../config.h"
#include "MigrationPolicy.h"
#include "MigrationScheme.h"

namespace pagmo
{

/// This object gathers migration-related parameters necessary to construct an archipelago.
/**
 * Introduction of this object has been lobbied by Dario because in his opinion its presence will increase
 * the user-friendliness of the Python interface.
 */
class __PAGMO_VISIBLE Migration
{
	public:
		/// Default constructor. All class fields are initialised to the default values imposed by Dario.
		/**
		 * These default values are:
		 * <ul>
		 *   <li>Migration scheme: see default values in the MigrationScheme::MigrationScheme() constructor</li>
		 *   <li>Migration policy: see default values in the MigrationPolicy::MigrationPolicy() constructor</li>
		 * </ul>
		 * Cool, huh?
		 */
		Migration()
				:migrationScheme(MigrationScheme().clone()),
				migrationPolicy(new MigrationPolicy()) { }

		/// Constructor.
		/** This constructor allows specifying migration scheme and policy. */
		Migration(const MigrationScheme& ms, const MigrationPolicy& mp)
				:migrationScheme(ms.clone()),
				migrationPolicy(new MigrationPolicy(mp)) { }

		/// Copy constructor.
		/** Creates a deep copy of the given migration descriptor. */
		Migration(const Migration& migration)
				:migrationScheme((*(migration.migrationScheme)).clone()),
				migrationPolicy(new MigrationPolicy(*(migration.migrationPolicy))) { }

		// Getters and setters.
		const MigrationScheme& getMigrationScheme() const {
			return *migrationScheme;
		}

		void setMigrationScheme(const MigrationScheme& ms) {
			migrationScheme.reset(ms.clone());
		}


		MigrationPolicy& getMigrationPolicy() const {
			return *migrationPolicy;
		}

		void setMigrationPolicy(const MigrationPolicy& mp) {
			migrationPolicy.reset(new MigrationPolicy(mp));
		}

	private:

		/// The archipelago's migration scheme. Must not be null.
		boost::scoped_ptr<MigrationScheme> migrationScheme;

		/// The islands' migration policy. Must not be null.
		boost::scoped_ptr<MigrationPolicy> migrationPolicy;

		/// Stream output operator.
		friend __PAGMO_VISIBLE_FUNC std::ostream &operator<<(std::ostream &s, const Migration& migration);
};

/// The stream output operator.
__PAGMO_VISIBLE_FUNC std::ostream &operator<<(std::ostream &s, const Migration& migration);

}

#endif
