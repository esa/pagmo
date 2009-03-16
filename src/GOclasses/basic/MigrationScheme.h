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

#ifndef PAGMO_MIGRATION_SCHEME_H
#define PAGMO_MIGRATION_SCHEME_H

#include "island.h"
#include <boost/scoped_ptr.hpp>
#include "base_topology.h"
#include "../../config.h"

//In order to prevent recursive inclusion... TODO: move MigrationScheme to a separate file?
class island;

/// Base class for the migration schemes.
/**
 * Migration schemes are the actual implementation of migration. 
 */
class __PAGMO_VISIBLE MigrationScheme
{
	public:
		/// Constructor.
		/**
		 * Creates a migration shceme associated with the given topology.
		 * A deep copy of topology is created and stored.
		 */
		MigrationScheme(const base_topology& _topology)
				:topology(_topology.clone())
		{
		}
		
		/// Virtual destructor.
		virtual ~MigrationScheme() { }; 
	
		/// Pre-evolution callback.
		/**
		 * This method is called by islands just before the actual evolution starts.
		 * The method must be thread-safe, but it may assume that the island being the argument of the function is not
		 * evolving at the moment.
		 * \param[in,out] _island island calling the function.
		 */		
		virtual void preEvolutionCallback(island& _island) = 0;
		
		/// Post-evolution callback.
		/**
		 * This method is called by islands right after the actual evolution finishes.
		 * The method must be thread-safe, but it may assume that the island being the argument of the function is not
		 * evolving at the moment.
		 * \param[in,out] _island island calling the function.
		 */		
		virtual void postEvolutionCallback(island& _island) = 0;

		///Clone object.
		/**
		 * Cloned object should be the exact copy of the original, but should be safe to use in a multithreaded environment.
		 * \todo Determine if the state variables should also be cloned and a separate method for resetting them be provided, or not.
		 */
		virtual MigrationScheme* clone() const = 0;
		
	private:
		
		// Class fields.
		boost::scoped_ptr<base_topology>	topology; ///< Migration topology. \todo I'm not so sure if all possible migration schemes require a topology...
};

#endif
