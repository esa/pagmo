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

#ifndef PAGMO_MIGRATION_REPLACEMENT_POLICY_H
#define PAGMO_MIGRATION_REPLACEMENT_POLICY_H

#include <boost/shared_ptr.hpp>
#include <list>
#include <utility>
#include "population.h"

/// Base class for replacement policies for migration.
class MigrationReplacementPolicy
{
	public:
		/// Virtual destructor.
		virtual ~MigrationReplacementPolicy() { }
		
		/// Assign pairs of individuals for replacement during migration. 
		/**
		 * Note, that this method does not alter any of the populations, it just provides the replacement choice.
		 * The actual replacement is done in another place.
		 * The first element of a pair should be the one from the destination population.
		 * The second one - an individual from the incoming population which is to replace the first one.
		 * \param[in] incomingPopulation Population of incoming individuals.
		 * \param[in] destinationPopulation Destination population.
		 * \return Replacement assignment.
		 */
		virtual std::list<std::pair<int, int> > selectForReplacement(const Population& incomingPopulation, const Population& destinationPopulation) = 0;
		
		///Clone object.
		/**
		 * Cloned object should be the exact copy of the original, but should be safe to use in a multithreaded environment.
		 * \todo Determine if the state variables should also be cloned and a separate method for resetting them be provided, or not.
		 */
		virtual MigrationReplacementPolicy* clone() const = 0;
};

#endif
