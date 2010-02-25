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
#include <sstream>
#include <string>

#include "algorithm/base.h"
#include "archipelago.h"
#include "problem/base.h"
#include "exceptions.h"
#include "island.h"
#include "population.h"
#include "types.h"

namespace pagmo
{
/// Constructor from problem::base, algorithm::base and number of individuals.
/**
 * Will store a copy of the problem and of the algorithm internally, will initialise internal population to n individuals
 * and evolution time to zero. Will fail if n is negative.
 *
 * @param[in] p problem::base to which the internal population will be associated.
 * @param[in] a algorithm::base which will be associated to the island.
 * @param[in] n number of individuals in the internal population.
 */
island::island(const problem::base &p, const algorithm::base &a, int n):m_pop(p,n),m_algo(a.clone()),m_archi(0),m_evo_time(0) {}

/// Copy constructor.
/**
 * Will perform a deep copy of all the elements of island isl, which will be synchronised before any operation takes place.
 *
 * @param[in] isl island to be copied.
 */
island::island(const island &isl)
{
	// Do it like this so that we can synchronise isl before poking into its internals.
	operator=(isl);
}

/// Assignment operator.
/**
 * Performs a deep copy of all the elements of isl into this island. Both island will be synchronised before assignment.
 *
 * @param[in] isl island used for assignment.
 *
 * @return reference to this.
 */
island &island::operator=(const island &isl)
{
	if (this != &isl) {
		// Make sure both islands are in a known state.
		join();
		isl.join();
		// Copy over content.
		m_pop = isl.m_pop;
		m_algo = isl.m_algo->clone();
		m_archi = isl.m_archi;
		m_evo_time = isl.m_evo_time;
	}
	return *this;
}

/// Destructor.
/**
 * Will call island::join() before returning. No other side effects.
 */
island::~island()
{
	join();
}

/// Synchronise the island.
/**
 * After this method returns, any pending evolution has been completed. This method is called internally
 * by all public methods.
 */
void island::join() const
{
	lock_type lock(m_evo_mutex);
}

/// Return terse human readable representation of the island.
/**
 * Will return a formatted string containing:
 * - description of the algorithm,
 * - the output of population::human_readable_terse().
 *
 * @return string containing terse human readable representation of the island.
 */
std::string island::human_readable_terse() const
{
	join();
	std::ostringstream oss;
	oss << *m_algo << '\n';
	oss << m_pop.human_readable_terse() << '\n';
	return oss.str();
}

/// Return human readable representation of the island.
/**
 * Will return a formatted string containing:
 * - description of the algorithm,
 * - the output of population::human_readable().
 *
 * @return string containing complete human readable representation of the island.
 */
std::string island::human_readable() const
{
	join();
	std::ostringstream oss;
	oss << *m_algo << '\n';
	oss << m_pop;
	return oss.str();
}

/// Return the total evolution time in milliseconds.
/**
 * Note that on many 32-bit machines this counter will wrap after roughly 49 days. On most 64-bit machines, the wrapping time
 * will be around 584 million years.
 *
 * @return number of milliseconds spent evolving the population.
 */
std::size_t island::evolution_time() const
{
	join();
	return m_evo_time;
}

/// Return copy of the internal algorithm.
/**
 * @return algorithm::base_ptr to the cloned algorithm.
 */
algorithm::base_ptr island::get_algorithm() const
{
	join();
	return m_algo->clone();
}

/// Query the state of the island.
/**
 * @return true if island is evolving, false otherwise.
 */
bool island::busy() const
{
	if (!m_evo_mutex.try_lock()) {
		return true;
	}
	m_evo_mutex.unlock();
	return false;
}

/// Algorithm setter.
/**
 * The input algorithm will be cloned.
 *
 * @param[in] algo new algorithm::base for the island.
 */
void island::set_algorithm(const algorithm::base &algo)
{
	join();
	m_algo = algo.clone();
}

/// Return copy of the internal problem.
/**
 * @return problem::base_ptr to the cloned problem.
 */
problem::base_ptr island::get_problem() const
{
	join();
	return m_pop.problem().clone();
}

/// Size of the internal population.
/**
 * @return size of the internal population.
 */
population::size_type island::get_size() const
{
	join();
	return m_pop.size();
}

// Implementation of int_evolver's juice.
void island::int_evolver::operator()()
{
	const boost::posix_time::ptime start = boost::posix_time::microsec_clock::local_time();
	try {
		// Synchronise start with all other threads
		if (m_i->m_archi) {
			// TODO: restore.
			// m_i->m_a->sync_island_start();
		}
		for (std::size_t i = 0; i < m_n; ++i) {
			if (m_i->m_archi) {
				// TODO: restore.
				m_i->m_archi->pre_evolution(*m_i);
				//m_i->m_pop.problem().pre_evolution(m_i->m_pop);
			}
			// Call the evolution.
			m_i->m_algo->evolve(m_i->m_pop);
			if (m_i->m_archi) {
				// TODO: restore.
				m_i->m_archi->post_evolution(*m_i);
				//m_i->m_pop.problem().post_evolution(m_i->m_pop);
			}
		}
	} catch (const std::exception &e) {
		std::cout << "Error during island evolution: " << e.what() << '\n';
	} catch (...) {
		std::cout << "Error during island evolution, unknown exception caught. :(\n";
	}
	// We must take care of potentially low-accuracy clocks, where the time difference could be negative for
	// _really_ short evolution times. In that case do not add anything to the total evolution time.
	const boost::posix_time::time_duration diff = boost::posix_time::microsec_clock::local_time() - start;
	if (diff.total_milliseconds() >= 0) {
		m_i->m_evo_time += boost::numeric_cast<std::size_t>(diff.total_milliseconds());
	}
	m_i->m_evo_mutex.unlock();
}

/// Evolve island n times.
/**
 * Open a thread and call the internal algorithm's algorithm::base::evolve() method n times on the internal population. Will fail if n is negative.
 *
 * This method will return as soon as it has started the first evolution run. During evolution, the island is locked down and no actions on it are possible,
 * but the flow of the rest of the program can continue without waiting for all evolutions to finish. To explicitly block the program until all evolution runs
 * have been performed on the island, call the join() method.
 *
 * @param[in] n number of algorithm::base::evolve() calls that will be performed by the internal algorithm on the population.
 */
void island::evolve(int n)
{
	const std::size_t n_evo = boost::numeric_cast<std::size_t>(n);
	if (m_evo_mutex.try_lock()) {
		try {
			boost::thread(int_evolver(this,n_evo));
		} catch (...) {
			pagmo_throw(std::runtime_error,"failed to launch the thread");
		}
	} else {
		pagmo_throw(std::runtime_error,"cannot evolve while still evolving");
	}
}

// Perform at least one evolution, and continue evolving until at least a certain amount of time has passed.
void island::t_evolver::operator()()
{
	const boost::posix_time::ptime start = boost::posix_time::microsec_clock::local_time();
	boost::posix_time::time_duration diff;
	try {
		// Synchronise start
		if (m_i->m_archi) {
			// TODO: restore stuff here.
			//m_i->m_a->sync_island_start();
		}
		do {
			if (m_i->m_archi) {
				m_i->m_archi->pre_evolution(*m_i);
				//m_i->m_pop.problem().pre_evolution(m_i->m_pop);
			}
			m_i->m_algo->evolve(m_i->m_pop);
			if (m_i->m_archi) {
				m_i->m_archi->post_evolution(*m_i);
				//m_i->m_pop.problem().post_evolution(m_i->m_pop);
			}
			diff = boost::posix_time::microsec_clock::local_time() - start;
			// Take care of negative timings.
		} while (diff.total_milliseconds() < 0 || boost::numeric_cast<std::size_t>(diff.total_milliseconds()) < m_t);
	} catch (const std::exception &e) {
		std::cout << "Error during evolution: " << e.what() << '\n';
	} catch (...) {
		std::cout << "Unknown exception caught. :(\n";
	}
	m_i->m_evo_time += boost::numeric_cast<std::size_t>(diff.total_milliseconds());
	m_i->m_evo_mutex.unlock();
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
	const std::size_t t_evo = boost::numeric_cast<std::size_t>(t);
	if (m_evo_mutex.try_lock()) {
		try {
			boost::thread(t_evolver(this,t_evo));
		} catch (...) {
			pagmo_throw(std::runtime_error,"failed to launch the thread");
		}
	} else {
		pagmo_throw(std::runtime_error,"cannot evolve while still evolving");
	}
}

/// Overload stream operator for pagmo::island.
/**
 * Equivalent to printing island::human_readable() to stream.
 *
 * @param[in] s stream to which the island will be sent.
 * @param[in] isl island to be sent to stream.
 *
 * @return reference to s.
 */
std::ostream &operator<<(std::ostream &s, const island &isl)
{
	s << isl.human_readable();
	return s;
}

}
