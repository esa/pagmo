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

#ifndef PAGMO_MIGRATION_H
#define PAGMO_MIGRATION_H

#include "population.h"
#include "individual.h"
#include <list>
#include <utility>
#include <boost/shared_ptr.hpp>

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
		virtual boost::shared_ptr<Population> selectForMigration(const Population& population) = 0;
};

/// Dummy selection policy.
/**
 * This policy does no selection. Use it as the defualt value, for example for islands that are not connected
 * to any archipelago.
 */
class DummySelectionPolicy : public MigrationSelectionPolicy
{
	public:
		/// \see MigrationSelectionPolicy::selectForMigration
		virtual boost::shared_ptr<Population> selectForMigration(const Population& population);
};




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
		virtual boost::shared_ptr<std::list<std::pair<int, int> > > selectForReplacement(const Population& incomingPopulation, const Population& destinationPopulation) = 0;
};

/// Dummy replacement policy.
/**
 * This policy does nothing. Use it as the defualt vlaue, for example for islands that are not connected
 * to any archipelago.
 */
class DummyReplacementPolicy : public MigrationReplacementPolicy
{
	public:
		/// \see MigrationReplacementPolicy::selectForReplacement
		virtual boost::shared_ptr<std::list<std::pair<int, int> > > selectForReplacement(const Population& incomingPopulation, const Population& destinationPopulation);
};



//In order to prevent recursive inclusion... TODO: move MigrationScheme to a separate file?
class island;

/// Base class for the migration schemes.
/**
 * Migration schemes are the actual implementation of migration. 
 */
class MigrationScheme
{
	public:
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
};

/// "No migration" policy.
class NoMigrationScheme : public MigrationScheme
{
	public:
		/// \see MigrationScheme::preEvolutionCallback
		virtual void preEvolutionCallback(island& _island) { };
		
		/// \see MigrationScheme::postEvolutionCallback
		virtual void postEvolutionCallback(island& _island) { };
};

#endif
