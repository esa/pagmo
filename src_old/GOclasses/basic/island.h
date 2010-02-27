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

// 04/01/2009: Initial version by Francesco Biscani.

#ifndef PAGMO_ISLAND_H
#define PAGMO_ISLAND_H

#include <boost/scoped_ptr.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <iostream>

#include "../../config.h"
#include "../../atomic_counters/atomic_counters.h"
#include "../algorithms/base.h"
#include "../problems/base.h"
#include "individual.h"
#include "migration/MigrationPolicy.h"
#include "population.h"

namespace pagmo
{

class archipelago;

/// Island class.
/**
 * An island encorporates a population of individuals, evolved using a specific algorithm.
 * It is usually a part of archipelago and possess some specific properties like migration
 * selection/replacement policies.
 * Note, that many methods of this class are synchronised, i.e. their execution waits until the eventual running
 * evolution finishes.
 */
class __PAGMO_VISIBLE island
{
		/// Mutex type abbreviation.
		typedef boost::mutex mutex_type;
		/// Lock guard type abbreviation.
		typedef boost::lock_guard<mutex_type> lock_type;

		/// Friend class... dirty?
		friend class archipelago;

		/// Stream output operator.
		friend __PAGMO_VISIBLE_FUNC std::ostream &operator<<(std::ostream &, const island &);

		/// Stream output operator (for the archipelago class).
		friend __PAGMO_VISIBLE_FUNC std::ostream &operator<<(std::ostream &s, const archipelago &a);

	public:
		/// Constructor.
		/**
		 * Creates an island associated with a specified problem and using the specified algorithm.
		 * The island is not associated to any archipelago, the island's population is empty,
		 * the total evolution time is reset to 0 and the island uses dummy selection and replacement policies.
		 * \param[in] p Problem to be associated with the island.
		 * \param[in] al Algorithm to be used by the island.
		 */
		island(const problem::base& p, const algorithm::base& al);

		/// Constructor.
		/**
		 * Creates an island associated with a specified problem, using the specified algorithm and having
		 * an initial population of given size.
		 * The island is not associated to any archipelago, the total evolution time is reset to 0
		 * and the island uses the dummy selection and replacement policies.
		 * \param[in] p Problem to be associated with the island.
		 * \param[in] al Algorithm to be used by the island.
		 * \param[in] n Size of the island's population.
		 */
		island(const problem::base& p, const algorithm::base& al, int n);

		/// Constructor.
		/**
		 * Creates an island associated with a specified problem, using the specified algorithm, having
		 * an initial population of given size and using specified migration parameters.
		 * Deep copies of the migration policy is created and internally stored (which allows safe multi-thread operation).
		 * \param[in] p Problem to be associated with the island.
		 * \param[in] al Algorithm to be used by the island.
		 * \param[in] n Size of the island's population.
		 * \param[in] mp Migration policy.
		 */
		island(const problem::base& p, const algorithm::base& al, int n, const MigrationPolicy& mp);

		/// Copy constructor.
		/**
		 * Creates a deep copy of an island. The only differences are that the new island is assigned a new id
		 * and is not connected to any archipelago.
		 * \todo Maybe make an automatic insertion to archipelago, if the original island is connected?
		 * \param i island to be duplicated.
		 */
		island(const island& i);

		/// Assignment operator.
		/**
		 * Creates a copy of an island, dropping current contents of the left-side argument.
		 * The only differences are that the left-side island id and the assiciated archipelago are not changed.
		 */
		island &operator=(const island &);

		/// Destructor (<b>synchronised</b>).
		~island();


		//Getters and setters

		/// Population getter (<b>synchronised</b>).
		population get_population() const;

		/// Problem getter (<b>synchronised</b>).
		const problem::base &problem() const;

		/// Algorithm getter (<b>synchronised</b>).
		const algorithm::base &algorithm() const;

		/// Algorithm seter (<b>synchronised</b>).
		void set_algorithm(const algorithm::base &);


		/// Migration selection policy public getter (<b>synchronised</b>).
		/** Note, that this function will throw an exception if no policy is associated with the island. */
		const MigrationSelectionPolicy& getMigrationSelectionPolicy() const;

		/// Migration selection policy public setter (<b>synchronised</b>).
		/**
		 * A deep copy of the given object is created and stored.
		 * \param[in] msp Selection policy to be used with the island. May not be null.
		 */
		void setMigrationSelectionPolicy(const MigrationSelectionPolicy& msp);

		/// Migration replacement policy getter (<b>synchronised</b>).
		/**  Note, that this function will throw an exception if no policy is associated with the island. */
		const MigrationReplacementPolicy& getMigrationReplacementPolicy() const;

		/// Migration replacement policy setter (<b>synchronised</b>).
		/**
		 * A deep copy of the given object is created and stored.
		 * \param[in] msp Replacement policy to be used with the island. May not be null.
		 */
		void setMigrationReplacementPolicy(const MigrationReplacementPolicy& mrp);

		/// Migration policy getter (<b>synchronised</b>).
		const MigrationPolicy& getMigrationPolicy() const;

		/// Migration policy setter (<b>synchronised</b>).
		/**
		 * A deep copy of the given object is created and stored.
		 * \param[in] mp Migration policy to be used with the island. May be null.
		 */
		void setMigrationPolicy(const MigrationPolicy* mp);

		/// Migration probability getter.
		/**
		 * This is provided here for convenience - it allows avoiding bothersome checking for migration policy presence,
		 * which would require introducing friend classes etc. Clear indication that the code starts getting a bit spaghettish...
		 */
		double getMigrationProbability() const {
			return migrationPolicy ? migrationPolicy->getMigrationProbability() : 0.0;
		}

		/// Get the island size (the number of individuals) (<b>synchronised</b>).
		size_t size() const;

		/// Get the island id (<b>synchronised</b>).
		size_t id() const;

		/// Get the total time spent by the island on evolution (<b>synchronised</b>).
		size_t evo_time() const;


		//Collection interface functions

		/// Individual indexed access operator (<b>synchronised</b>).
		individual operator[](int) const;

		/// Individual setter (<b>synchronised</b>).
		/** \see population::setIndividual */
		void set_individual(int, const individual &);

		/// Push back operation for underlying population (<b>synchronised</b>).
		/** \see population::push_back */
		void push_back(const individual &);

		/// Insert operation for underlying population (<b>synchronised</b>).
		/** \see population::insert */
		void insert(int, const individual &);

		/// Erase operation for underlying population (<b>synchronised</b>).
		/** \see population::erase */
		void erase(int);


		//Utility functions

		/// Calculate the mean fitness of the individuals (<b>synchronised</b>).
		double mean() const;

		/// Calculate the mean fitness of the individuals (<b>synchronised</b>).
		double std() const;

		/// Get the best individual (<b>synchronised</b>).
		individual best() const;

		/// Get the best worst (<b>synchronised</b>).
		individual worst() const;


		// Migration Functions

		/// Get individuals migrating out of the population.
		/**
		 * The method returns a population of individuals selected to migrate according to the island's policy.
		 * the individuals are copied, i.e. they are not removed from the population.
		 * \return Population of migrating individuals.
		 */
		std::vector<individual> getMigratingIndividuals();

		/// Put migrating individuals to the population.
		/**
		 * The island replaces individuals in it's own population with incoming individuals, according the it's policy.
		 * \param[in] incomingPopulation incoming individuals.
		 */
		void acceptMigratingIndividuals(const std::vector<individual>& incomingPopulation);


		// Evolution functions

		/// Start an evolution fo the specified number of algorithm steps.
		/**
		 * What exactly a 'step' means depends on the implementation of the algorithm the island is using.
		 * Evolution will be launched in a separate thread. The method will fail with an exception if the
		 * evolution is already running.
		 * \param[in] n the number of algorithm steps to perform
		 */
		void evolve(int n = 1);

		/// Start an evolution fo the specified amount of time.
		/**
		 * The algorithm is run for at least the specified amount of time and for at least one iteration.
		 * \param[in] t time for the algorithm to run, in miliseconds.
		 * Evolution will be launched in a separate thread. The method will fail with an exception if the
		 * evolution is already running.
		 */
		void evolve_t(const size_t& t);


		//Synchronisation functions

		/// Wait until the evolution finishes.
		void join() const;

		/// Check the evolution status.
		/**
		 * \return true, if the evolution is still running, false otherwise.
		 */
		bool busy() const;

		/*
		bool t_substitute_worst(const individual &);
		individual t_best() const;
		*/

	private:

		/// Archipelago setter.
		/**
		 * This method is to be used only by the archipelago class, which is the friend class.
		 */
		void set_archipelago(archipelago *);

		// void t_check() const;

		/// Evolver thread object.
		/**
		 * This is a callable helper object used to launch an evolution for a given number of iterations.
		 */
		struct int_evolver {

			/// Constructor.
			/**
			 * \param[in] i the owner island.
			 * \param[in] n the number of algorithm steps to perform.
			 */
			int_evolver(island *i, int n):m_i(i),m_n(n) { }

			/// Call operator.
			/**
			 * Does the actual evolution, including optional calls to the migration scheme.
			 */
			void operator()();

			island 		*m_i; ///< Owner island.
			const int	m_n; ///< Number of steps to perform.
		};

		/// Time-dependent evolver thread object.
		/**
		 * This is a callable helper object used to launch an evolution for a specified amount of time.
		 */
		struct t_evolver {

			/// Constructor.
			/**
			 * \param[in] i the owner island.
			 * \param[in] t time to evolve, in miliseconds.
			 */
			t_evolver(island *i, const size_t &t):m_i(i),m_t(t) {}

			/// Call operator.
			/**
			 * Does the actual evolution, including optional calls to the migration scheme.
			 */
			void operator()();

			island 			*m_i; ///< Owner island.
			const size_t	m_t; ///< Number of steps to perform.
		};


		/// Generates a new uniqe island id.
		static size_t get_new_id();

		static atomic_counter_size_t			id_counter; ///< Counter used to generate the island ids.

		//Class fields
		size_t										m_id; ///< Island id.
		population									m_pop; ///< Island's population.
		boost::scoped_ptr<const algorithm::base>		m_goa; ///< Island's algorithm.
		archipelago									*m_a; ///< Associated archipelago (may be null).
		size_t										m_evo_time; ///< Counts the total time spent by the island on evolution (in milliseconds).
		mutable mutex_type							m_evo_mutex; ///< Mutex used to control evolution synchronisation.
		boost::scoped_ptr<MigrationPolicy>          migrationPolicy; ///< Migration parameters of the island (may be null).
};

/// Stream output operator.
std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const island &);

}

#endif
