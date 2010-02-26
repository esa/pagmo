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

#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <cstddef>
#include <iostream>
#include <string>

#include "config.h"
#include "algorithm/base.h"
#include "migration/base.h"
#include "migration/base_s_policy.h"
#include "migration/best_s_policy.h"
#include "population.h"
#include "problem/base.h"
#include "types.h"

namespace pagmo
{

// Forward declaration of archipelago class, needed to make friend.
class archipelago;

/// Island class.
/**
 * This class incorporates a pagmo::population and a pagmo::algorithm::base used to evolve the population. Each time the evolve() (or evolve_t()) method is called, a
 * local thread is opened and the method returns immediately, while the population is evolved asynchronously in the background. While evolution is undergoing, the island is locked down
 * and no further operations will be allowed. The method join() can be used to wait until evolution on the island has terminated. The busy() methods can be used to query the state
 * of the island.
 *
 * If the island belongs to an archipelago, it can exchange individuals with other islands in the archipelago. The topology of the archipelago determines
 * the connections between islands, whereas every island can define migration policies to specify how to select and replace individuals during migration.
 * The relevant policy classes are migration::base_s_policy (selection policy) and migration::base_r_policy (replacement policy).
 *
 * The interface of this class mirrors the interface of the population class. It is hence possible to get and set individuals, get the population size,
 * access the population champion, etc. The main difference
 * is that the methods of this class will never return references to internal members, in order to protect the internal state of the island while evolution is undergoing.
 * All getters methods will thus return copies instead of references, and all public methods will wait for an ongoing evolution to terminate before performing any action.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE island
{
		// Lock type alias.
		typedef boost::lock_guard<boost::mutex> lock_type;
	public:
		/// The archipelago class needs access to the internals of the island.
		friend class archipelago;
		island(const island &);
		island(const problem::base &, const algorithm::base &, int n = 0,
			const migration::base_s_policy &s_policy = migration::best_s_policy(1,migration::absolute));
		island &operator=(const island &);
		~island();
		/** @name Input/output.*/
		//@{
		std::string human_readable_terse() const;
		std::string human_readable() const;
		//@}
		/** @name Evolution.
		 * Methods related to island evolution.
		 */
		//@{
		void join() const;
		bool busy() const;
		void evolve(int n = 1);
		void evolve_t(int);
		std::size_t get_evolution_time() const;
		//@}
		/** @name Getters and setters.*/
		//@{
		algorithm::base_ptr get_algorithm() const;
		void set_algorithm(const algorithm::base &);
		problem::base_ptr get_problem() const;
		population::size_type get_size() const;
		//@}
	private:
		// Evolver thread object. This is a callable helper object used to launch an evolution for a given number of iterations.
		struct int_evolver {
			int_evolver(island *i, const std::size_t &n):m_i(i),m_n(n) { }
			void operator()();
			island 			*m_i;
			const std::size_t	m_n;
		};
		// Time-dependent evolver thread object. This is a callable helper object used to launch an evolution for a specified amount of time.
		struct t_evolver {
			t_evolver(island *i, const std::size_t &t):m_i(i),m_t(t) {}
			void operator()();
			island 			*m_i;
			const std::size_t	m_t;
		};
	private:
		// Population.
		population			m_pop;
		// Algorithm.
		algorithm::base_ptr		m_algo;
		// Archipelago that, if not null, contains the island.
		archipelago			*m_archi;
		// Counts the total time spent by the island on evolution (in milliseconds).
		std::size_t			m_evo_time;
		// Mutex used to control evolution synchronisation.
		mutable boost::mutex		m_evo_mutex;
		// Migration selection policy.
		migration::base_s_policy_ptr	m_s_policy;
};

std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const island &);

}

#endif
