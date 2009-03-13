/*****************************************************************************
 *   Copyright (C) 2008, 2009 Advanced Concepts Team (European Space Agency) *
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

#include "population.h"

/// Base class for selection policies for migration.
class MigrationSelectionPolicy
{
	public:
		/// Virtual destructor.
		virtual ~MigrationSelectionPolicy() { }
		
		/// Select individuals to migrate out of the given population.
		/**
		 * This is the method that actually implements the policy.
		 * Output population should contain <b>copies</b> of selected individuals.
		 * \param[in] population Source population.
		 * \return Population containing selected individuals.
		 */
		virtual Population selectForMigration(const Population& population) = 0;
		
		///Clone object.
		/**
		 * Cloned object should be the exact copy of the original, but should be safe to use in a multithreaded environment.
		 * \todo Determine if the state variables should also be cloned and a separate method for resetting them be provided, or not.
		 */
		virtual MigrationSelectionPolicy* clone() const = 0;
};

#endif
