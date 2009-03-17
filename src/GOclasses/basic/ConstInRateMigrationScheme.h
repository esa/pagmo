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

#ifndef PAGMO_CONST_IN_RATE_MIGRATION_SCHEME_H
#define PAGMO_CONST_IN_RATE_MIGRATION_SCHEME_H

#include "MigrationScheme.h"
#include <boost/unordered_map.hpp>
#include <vector>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>

/// Migration scheme with constant incoming migration rate.
/**
 * This is the migration model in which the migration is "destination driven".
 */
class __PAGMO_VISIBLE ConstInRateMigrationScheme: public MigrationScheme
{
		/// Lock guard type abbreviation.
		typedef boost::lock_guard<boost::mutex> lock_type;

	public:
		/// Constructor.
		/**
		 * Creates a migration shceme associated with the given topology.
		 */
		ConstInRateMigrationScheme(const base_topology& _topology)
				:MigrationScheme(_topology)
		{
		}
		
		/// Copy constructor.
		/**
		 * Creates a deep copy of the object
		 */
		ConstInRateMigrationScheme(const ConstInRateMigrationScheme& ms):MigrationScheme(ms) { }
		
		/// Virtual destructor.
		virtual ~ConstInRateMigrationScheme() { }; 
	
		/// \see MigrationScheme::preEvolutionCallback
		virtual void preEvolutionCallback(island& _island);
		
		/// \see MigrationScheme::postEvolutionCallback
		virtual void postEvolutionCallback(island& _island);

		/// \see MigrationScheme::clone
		virtual MigrationScheme* clone() const { return new ConstInRateMigrationScheme(*this); }
		
	private:
		/// Database of the migrating individuals.
		/**
		 * Maps island's id to a vector of the most recent migrating individuals from that island.
		 */
		boost::unordered_map<size_t, std::vector<Individual> > immigrants;
		
		/// Mutex for protecting the database.
		boost::mutex							immigrants_mutex;		
};

#endif
