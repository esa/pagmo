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

#include <boost/numeric/conversion/cast.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/thread/thread.hpp>
#include <cstddef>
#include <exception>
#include <stdexcept>

#include "algorithm/base.h"
#include "base_island.h"
#include "exceptions.h"
#include "island.h"
#include "migration/base_r_policy.h"
#include "migration/base_s_policy.h"
#include "population.h"
#include "problem/base.h"

namespace pagmo
{
/// Constructor from problem::base, algorithm::base, number of individuals, migration probability and selection/replacement policies.
/**
 * @see pagmo::base_island constructors.
 */
island::island(const problem::base &p, const algorithm::base &a, int n, const double &migr_prob,
	const migration::base_s_policy &s_policy, const migration::base_r_policy &r_policy):
	base_island(p,a,n,migr_prob,s_policy,r_policy)
{}

/// Copy constructor.
/**
 * @see pagmo::base_island constructors.
 */
island::island(const island &isl):base_island(isl)
{}

/// Constructor from population.
/**
 * @see pagmo::base_island constructors.
 */
island::island(const population &pop, const algorithm::base &a, const double &migr_prob,
	const migration::base_s_policy &s_policy, const migration::base_r_policy &r_policy):
	base_island(pop,a,migr_prob,s_policy,r_policy)
{}

/// Assignment operator.
island &island::operator=(const island &isl)
{
	base_island::operator=(isl);
	if (this != &isl) {
		m_evo_thread.reset(0);
	}
	return *this;
}

/// Destructor.
/**
 * Will call island::join() before returning. No other side effects.
 */
island::~island()
{
	// TODO: Verify this.
	island::join();
}

void island::join() const
{
	if (m_evo_thread) {
		m_evo_thread->join();
	}
}

// Evolver thread object. This is a callable helper object used to launch an evolution for a given number of iterations.
struct island::int_evolver {
	int_evolver(island *i, const std::size_t &n, bool blocking):m_i(i),m_n(n),m_blocking(blocking) { }
	void operator()();
	void juice_impl(boost::posix_time::ptime &);
	island 			*m_i;
	const std::size_t	m_n;
	const bool		m_blocking;
};

void island::int_evolver::juice_impl(boost::posix_time::ptime &start)
{
	start = boost::posix_time::microsec_clock::local_time();
	// Synchronise start with all other threads if we are not blocking and we are in an archi.
	if (m_i->m_archi && !m_blocking) {
		// TODO: remove this sync_island_start?
		m_i->m_archi->sync_island_start();
	}
	for (std::size_t i = 0; i < m_n; ++i) {
		// Call pre-evolve hooks.
		if (m_i->m_archi) {
			m_i->m_archi->pre_evolution(*m_i);
		}
		m_i->m_pop.problem().pre_evolution(m_i->m_pop);
		// Call the evolution.
		m_i->m_algo->evolve(m_i->m_pop);
		// Post-evolve hooks.
		if (m_i->m_archi) {
			m_i->m_archi->post_evolution(*m_i);
		}
		m_i->m_pop.problem().post_evolution(m_i->m_pop);
		// If we are running in a separate thread, set the interruption point.
		if (!m_blocking) {
			boost::this_thread::interruption_point();
		}
	}
}

// Implementation of int_evolver's juice.
void island::int_evolver::operator()()
{
	boost::posix_time::ptime start;
	if (m_blocking) {
		juice_impl(start);
	} else {
		try {
			juice_impl(start);
		} catch (const boost::thread_interrupted &) {
			// In case of interruption, don't do anything special.
		} catch (const std::exception &e) {
			std::cout << "Error during island evolution: " << e.what() << '\n';
		} catch (...) {
			std::cout << "Error during island evolution, unknown exception caught. :(\n";
		}
	}
	// Try to compute the evolution time before exiting. In case something goes wrong, do not do anything.
	try {
		// We must take care of potentially low-accuracy clocks, where the time difference could be negative for
		// _really_ short evolution times. In that case do not add anything to the total evolution time.
		const boost::posix_time::time_duration diff = boost::posix_time::microsec_clock::local_time() - start;
		if (diff.total_milliseconds() >= 0) {
			m_i->m_evo_time += boost::numeric_cast<std::size_t>(diff.total_milliseconds());
		}
	} catch (...) {
		std::cout << "Error calculating evolution time.\n";
	}
}

void island::evolve(int n)
{
	join();
	const std::size_t n_evo = boost::numeric_cast<std::size_t>(n);
	if (is_blocking_impl() || (m_archi && m_archi->is_blocking_impl())) {
		int_evolver ev(this,n_evo,true);
		ev();
	} else {
		try {
			m_evo_thread.reset(new boost::thread(int_evolver(this,n_evo,false)));
		} catch (...) {
			pagmo_throw(std::runtime_error,"failed to launch the thread");
		}
	}
}

// Time-dependent evolver thread object. This is a callable helper object used to launch an evolution for a specified amount of time.
struct island::t_evolver {
	t_evolver(island *i, const std::size_t &t, bool blocking):m_i(i),m_t(t),m_blocking(blocking) {}
	void operator()();
	void juice_impl(boost::posix_time::ptime &);
	island 			*m_i;
	const std::size_t	m_t;
	const bool		m_blocking;
};

void island::t_evolver::juice_impl(boost::posix_time::ptime &start)
{
	boost::posix_time::time_duration diff;
	start = boost::posix_time::microsec_clock::local_time();
	// Synchronise start.
	if (m_i->m_archi && !m_blocking) {
		m_i->m_archi->sync_island_start();
	}
	do {
		if (m_i->m_archi) {
			m_i->m_archi->pre_evolution(*m_i);
		}
		m_i->m_pop.problem().pre_evolution(m_i->m_pop);
		m_i->m_algo->evolve(m_i->m_pop);
		if (m_i->m_archi) {
			m_i->m_archi->post_evolution(*m_i);
		}
		m_i->m_pop.problem().post_evolution(m_i->m_pop);
		// If we are running in a separate thread, set the interruption point.
		if (!m_blocking) {
			boost::this_thread::interruption_point();
		}
		diff = boost::posix_time::microsec_clock::local_time() - start;
		// Take care of negative timings.
	} while (diff.total_milliseconds() < 0 || boost::numeric_cast<std::size_t>(diff.total_milliseconds()) < m_t);
}

// Perform at least one evolution, and continue evolving until at least a certain amount of time has passed.
void island::t_evolver::operator()()
{
	boost::posix_time::ptime start;
	if (m_blocking) {
		juice_impl(start);
	} else {
		try {
			juice_impl(start);
		} catch (const boost::thread_interrupted &) {
			// In case of interruption, don't do anything special.
		} catch (const std::exception &e) {
			std::cout << "Error during evolution: " << e.what() << '\n';
		} catch (...) {
			std::cout << "Unknown exception caught. :(\n";
		}
	}
	// Try to compute the evolution time before exiting. In case something goes wrong, do not do anything.
	try {
		// We must take care of potentially low-accuracy clocks, where the time difference could be negative for
		// _really_ short evolution times. In that case do not add anything to the total evolution time.
		const boost::posix_time::time_duration diff = boost::posix_time::microsec_clock::local_time() - start;
		if (diff.total_milliseconds() >= 0) {
			m_i->m_evo_time += boost::numeric_cast<std::size_t>(diff.total_milliseconds());
		}
	} catch (...) {
		std::cout << "Error calculating evolution time.\n";
	}
}

/// Evolve island for a specified minimum amount of time.
/**
 * Open a thread and call the internal algorithm's algorithm::base::evolve() method on the population at least once, and keep calling it until at least t milliseconds
 * ("wall clock" time) have elapsed. Will fail if t is negative.
 *
 * This method will return as soon as it has started the first evolution run. During evolution, the island is locked down and no actions on it are possible,
 * but the flow of the rest of the program can continue without waiting for all evolutions to finish. To explicitly block the program until all evolution runs
 * have been performed on the island, call the join() method.
 *
 * @param[in] t minimum evolution time in milliseconds.
 */
void island::evolve_t(int t)
{
	join();
	const std::size_t t_evo = boost::numeric_cast<std::size_t>(t);
	if (is_blocking_impl() || (m_archi && m_archi->is_blocking_impl())) {
		t_evolver ev(this,t_evo,true);
		ev();
	} else {
		try {
			m_evo_thread.reset(new boost::thread(t_evolver(this,t_evo,false)));
		} catch (...) {
			pagmo_throw(std::runtime_error,"failed to launch the thread");
		}
	}
}

/// Interrupt evolution.
/**
 * If an evolution is undergoing, it will be stopped the first time it reaches one of the internal interruption points.
 */
void island::interrupt()
{
	if (m_evo_thread) {
		m_evo_thread->interrupt();
		pagmo_throw(std::runtime_error,"evolution interrupted");
	}
}

/// Query the status of the island.
/**
 * @return true if the island is evolving, false otherwise.
 */
bool island::busy() const
{
	if (!m_evo_thread) {
		return false;
	}
	return m_evo_thread->joinable();
}

// Blocking attribute implementation.
bool island::is_blocking_impl() const
{
	return (m_pop.problem().is_blocking() || m_algo->is_blocking());
}

/// Check if the island is blocking.
/**
 * @return true if either the algorihm or the problem are blocking.
 */
bool island::is_blocking() const
{
	join();
	return is_blocking_impl();
}

}
