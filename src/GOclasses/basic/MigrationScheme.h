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
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <iostream>
#include <boost/unordered_map.hpp>
#include "../../Functions/rng/rng.h"
#include <vector>

/// Base class for the migration schemes.
/**
 * Migration schemes are the actual implementation of migration. 
 * By the collective decision, the objective programming paradigm was sacrificed here for the sake of the ease of use in python.
 */
class __PAGMO_VISIBLE MigrationScheme
{
	protected:
		/// Lock guard type abbreviation.
		typedef boost::lock_guard<boost::mutex> lock_type;
	
	public:
		/// Constructor.
		/**
		 * Creates a migration scheme not associated with a topology.
		 */
		MigrationScheme(int _distributionType, int _migrationDirection, uint32_t seed = static_rng_uint32()())
			: distributionType(_distributionType), migrationDirection(_migrationDirection), topology(0), rng(seed) { }

		/// Constructor.
		/**
		 * Creates a migration scheme associated with the given topology.
		 * A deep copy of topology is created and stored.
		 */
		MigrationScheme(int _distributionType, int _migrationDirection, const base_topology& _topology, uint32_t seed = static_rng_uint32()())
				:distributionType(_distributionType), migrationDirection(_migrationDirection), topology(_topology.clone()), rng(seed)
		{
		}
		
		/// Copy constructor.
		/**
		 * Creates a deep copy of the object
		 */
		MigrationScheme(const MigrationScheme& ms)				
		{
			distributionType = ms.distributionType;
			migrationDirection = ms.migrationDirection;
			topology.reset(ms.topology ? ms.topology->clone() : 0);
			rng = ms.rng;
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
		void preEvolutionCallback(island& _island);
		
		/// Post-evolution callback.
		/**
		 * This method is called by islands right after the actual evolution finishes.
		 * The method must be thread-safe, but it may assume that the island being the argument of the function is not
		 * evolving at the moment.
		 * \param[in,out] _island island calling the function.
		 */		
		void postEvolutionCallback(island& _island);

		///Clone object.
		/**
		 * Cloned object should be the exact copy of the original, but should be safe to use in a multithreaded environment.
		 * \todo Determine if the state variables should also be cloned and a separate method for resetting them be provided, or not.
		 */
		MigrationScheme* clone() const { return new MigrationScheme(*this); }
		
		/// Register an island with the migration scheme.
		virtual void push_back(const island& _island)
		{
			lock_type lock(topology_mutex); ///\todo just in case... is it really needed?
			if(topology) {
				topology->push_back(_island.id());
			}
		}
		
		/// Reset the migration scheme.
		/**
		 * Associated topology, if present, is cleared. All islands have to be re-registered with the migration scheme.		 
		 */
		virtual void reset()
		{
			immigrantsByIsland.clear();
			
			if(topology) {
				topology->clear();
			}
		}		
		
		//Getters and setters
		/// Public topology getter.
		/**
		 * Note that this function will throw an exception if no topology is associated with the scheme.		 
		 */
		const base_topology& getTopology() const
		{
			if(!topology) {
				pagmo_throw(value_error, "The migration scheme has no associated topology!");
			}
			return *topology;
		}
		
		/// Topology setter.
		/**
		 * A deep copy of the given object is created and stored. Apart from this,the method has no embedded logic;
		 * it just performs a replacement of the associated topology. In consequence, the toplogy object must contain
		 * valid content.
		 * \todo Is this the right behaviour?
		 */
		void setTopology(const base_topology* _topology) { topology.reset(_topology ? _topology->clone() : 0); }		
		
	private:
		/// Immigrants distribution type flag.
		/** 0 = point-to-point, 1 = broadcast among neighbours */
		int distributionType;
		
		/// Migration direction flag.
		/** 0 = immigrants flow is initiated by the source island, 1 = immingrants flow is initiated by the destination island */
		int migrationDirection;
	
		boost::scoped_ptr<base_topology>	topology; ///< Migration topology. \todo I'm not so sure if all possible migration schemes require a topology...
		mutable boost::mutex				topology_mutex; ///< Topology mutex. <b>Access to the topology must be synchronised!!!</b>

		/// Container for migrating individuals.
		/**
		 * The way in which this container is used depends on the type of migration performed.
		 */
		boost::unordered_map<size_t, std::vector<Individual> > immigrantsByIsland;
		
		rng_uint32 rng; ///< Random number generator	
	
		/// Stream output operator.
		friend std::ostream &operator<<(std::ostream &s, const MigrationScheme& ms);
};

/// Stream output operator.
__PAGMO_VISIBLE_FUNC std::ostream &operator<<(std::ostream &s, const MigrationScheme& ms);

#endif
