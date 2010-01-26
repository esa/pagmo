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

#ifndef PAGMO_MIGRATION_SCHEME_H
#define PAGMO_MIGRATION_SCHEME_H

#include <boost/cstdint.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/unordered_map.hpp>
#include <iostream>
#include <vector>

#include "../../../config.h"
#include "../../../rng.h"
#include "../../basic/island.h"
#include "../topology/base_topology.h"
#include "../topology/one_way_ring_topology.h"

namespace pagmo
{

/// Base class for the migration schemes.
/**
 * Migration scheme is the actual implementation of migration.
 * By the collective decision, the objective programming paradigm was sacrificed here for the sake of the ease of use in python.
 *
 * In this class it is decided how migrating individuals are distributed, and by which island (source or destination) the migration is initiated.
 * So far, two ways of immigrants distribution are possible:
 * <ul>
 *   <li>Point-to-point: the migration takes place only between the initiating island and one randomly selected neighbour</li>
 *   <li>Broadcast: the migration takes place between the initiating island and all its neghbours</li>
 * </ul>
 * The type of distribution is set with the first argument of the class' constructor.
 *
 * Asynchronous migration may be initiated in two ways: either by the island being the "source" of migrating individuals in the processs of migration,
 * or by the one being the "destination" of the immigrants. The difference between the two cases may not be evident immidiately, but they
 * are not identical.
 *
 * Consider the case of a star-like (broadcast) topology with one node in the centre (the hub) and many nodes connected only to it (the leaves).
 * Consider the incoming and outgoing flow rates of immigrants for bot the hub and a leaf:
 * <ul>
 *   <li>When the migration is initiated by the source island, outgoing immigrant flow rate of the hub is contant and outgoing rate of the leaves is also constant.
 *   Incoming rates however are low for the leaves, especially in case of the point-to-point migration, as migration out of the hub occurs relatively rarely
 *   and the chance that the particular leaf is the migration target is very low. Symetrically, the incoming migration rate of the hub is very high, as many leaves
 *   migrate to it simultanously.</li>
 *   <li>When the migration is initiated by the destination island, the situation is reversed. Incoming immigrant flow rates are constant for both hub and the
 *   leaves, but the outgoing rate is very high for the hub and very low for the leaves (again, probability of choosing a particular leaf for migration is
 *   relatively low)</li>
 * </ul>
 * The type of migration initiation is governed by the second argument of the class' constructor.
 */
class __PAGMO_VISIBLE MigrationScheme
{
	protected:
		/// Lock guard type abbreviation.
		typedef boost::lock_guard<boost::mutex> lock_type;

	public:
		/// Default constructor.
		/**
		 * Creates a migration scheme with a one-way ring topology, using Dario's favourite defualt parameters.
		 * Their values are:
		 * <ul>
		 * <li>distributionType: point-to-point</li>
		 * <li>migrationDirection: destination</li>
		 * </ul>
		 * @param[in] seed RNG seed.
		 */
		MigrationScheme(boost::uint32_t seed = static_rng_uint32()())
				: distributionType(0), migrationDirection(1), topology(one_way_ring_topology().clone()), rng(seed) { }

		/// Almost default constructor.
		/**
		 * Creates a migration scheme with a given topology, using Dario's favourite defualt parameters.
		 * For values of these parameters see \see MigrationScheme::MigrationScheme()
		 * @param[in] _topology Topology to be used.
		 * @param[in] seed RNG seed.
		 */
		MigrationScheme(const base_topology& _topology, boost::uint32_t seed = static_rng_uint32()())
				: distributionType(0), migrationDirection(1), topology(_topology.clone()), rng(seed) { }


		/// Constructor.
		/**
		 * Creates a migration scheme not associated with a topology.
		 * @param _distributionType Immigrants distribution type: 0 = point-to-point, 1 = broadcast
		 * @param _migrationDirection Initiating island: 0 = source, 1 = destination
		 * @param seed RNG seed.
		 */
		MigrationScheme(int _distributionType, int _migrationDirection, boost::uint32_t seed = static_rng_uint32()())
				: distributionType(_distributionType), migrationDirection(_migrationDirection), topology(0), rng(seed) { }

		/// Constructor.
		/**
		 * Creates a migration scheme associated with the given topology.
		 * A deep copy of topology is created and stored.
		 * @param _distributionType Immigrants distribution type: 0 = point-to-point, 1 = broadcast
		 * @param _migrationDirection Initiating island: 0 = source, 1 = destination
		 * @param _topology Topology to be used for the migration.
		 * @param seed RNG seed.
		 */
		MigrationScheme(int _distributionType, int _migrationDirection, const base_topology& _topology, boost::uint32_t seed = static_rng_uint32()())
				:distributionType(_distributionType), migrationDirection(_migrationDirection), topology(_topology.clone()), rng(seed) {
		}

		/// Copy constructor.
		/**
		 * Creates a deep copy of the object
		 */
		MigrationScheme(const MigrationScheme& ms) {
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
		MigrationScheme* clone() const {
			return new MigrationScheme(*this);
		}

		/// Register an island with the migration scheme.
		virtual void push_back(const island& _island) {
			lock_type lock(topology_mutex); ///\todo just in case... is it really needed?
			if (topology) {
				topology->push_back(_island.id());
			}
		}

		/// Reset the migration scheme.
		/**
		 * Associated topology, if present, is cleared. All islands have to be re-registered with the migration scheme.
		 */
		virtual void reset() {
			immigrantsByIsland.clear();

			if (topology) {
				topology->clear();
			}
		}

		//Getters and setters
		/// Public topology getter.
		/**
		 * Note that this function will throw an exception if no topology is associated with the scheme.
		 */
		const base_topology& getTopology() const {
			if (!topology) {
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
		void setTopology(const base_topology* _topology) {
			topology.reset(_topology ? _topology->clone() : 0);
		}

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
		boost::unordered_map<size_t, std::vector<individual> > immigrantsByIsland;

		rng_double rng; ///< Random number generator

		/// Stream output operator.
		friend __PAGMO_VISIBLE_FUNC std::ostream &operator<<(std::ostream &s, const MigrationScheme& ms);
};

/// Stream output operator.
__PAGMO_VISIBLE_FUNC std::ostream &operator<<(std::ostream &s, const MigrationScheme& ms);

}

#endif
