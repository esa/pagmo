/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
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

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/thread/thread.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <cstddef>
#include <exception>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

#include "algorithm/base.h"
#include "archipelago.h"
#include "base_island.h"
#include "exceptions.h"
#include "migration/base_r_policy.h"
#include "migration/base_s_policy.h"
#include "population.h"
#include "problem/base.h"
#include "types.h"

namespace pagmo
{
/// Constructor from problem::base, algorithm::base, number of individuals, migration probability and selection/replacement policies.
/**
 * Will store a copy of the problem, of the algorithm and of the policies internally, will initialise internal population to n individuals
 * and evolution time to zero. Will fail if n is negative or if migration probability is not in the [0,1] range.
 *
 * @param[in] a algorithm::base which will be associated to the island.
 * @param[in] p problem::base to which the internal population will be associated.
 * @param[in] n number of individuals in the internal population.
 * @param[in] s_policy migration::base_s_policy for the island.
 * @param[in] r_policy migration::base_r_policy for the island.
 *
 * @throws pagmo::value_error if migration probability is outside the [0,1] range.
 */
base_island::base_island(const algorithm::base &a, const problem::base &p, int n,
	const migration::base_s_policy &s_policy, const migration::base_r_policy &r_policy):
	m_algo(a.clone()),m_pop(p,n),m_archi(0),m_evo_time(0),m_s_policy(s_policy.clone()),m_r_policy(r_policy.clone()) { }

/// Copy constructor.
/**
 * Will perform a deep copy of all the elements of island isl, which will be synchronised before any operation takes place.
 *
 * @param[in] isl island to be copied.
 */
base_island::base_island(const base_island &isl):m_pop(isl.get_population())
{
	// Population has already been done and get_population() above already called join().
	m_algo = isl.m_algo->clone();
	m_archi = isl.m_archi;
	m_evo_time = isl.m_evo_time;
	m_s_policy = isl.m_s_policy->clone();
	m_r_policy = isl.m_r_policy->clone();
	m_evo_thread.reset(0);
}

/// Constructor from population.
/**
 * Will construct an island containing the given population and algorithm.
 *
 * @param[in] a algorithm::base which will be associated to the island.
 * @param[in] pop population that will be contained in the island.
 * @param[in] s_policy migration::base_s_policy for the island.
 * @param[in] r_policy migration::base_r_policy for the island.
 *
 * @throws pagmo::value_error if migration probability is outside the [0,1] range.
 */
base_island::base_island(const algorithm::base &a, const population &pop,
	const migration::base_s_policy &s_policy, const migration::base_r_policy &r_policy):
	m_algo(a.clone()),m_pop(pop),m_archi(0),m_evo_time(0),m_s_policy(s_policy.clone()),m_r_policy(r_policy.clone()) { }

/// Assignment operator.
/**
 * Performs a deep copy of all the elements of isl into this island. Both island will be synchronised before assignment.
 *
 * @param[in] isl island used for assignment.
 *
 * @return reference to this.
 */
base_island &base_island::operator=(const base_island &isl)
{
	if (this != &isl) {
		// Make sure both islands are in a known state.
		join();
		isl.join();
		// Copy over content.
		m_algo = isl.m_algo->clone();
		m_pop = isl.m_pop;
		m_archi = isl.m_archi;
		m_evo_time = isl.m_evo_time;
		m_s_policy = isl.m_s_policy->clone();
		m_r_policy = isl.m_r_policy->clone();
		m_evo_thread.reset(0);
	}
	return *this;
}

/// Destructor.
/**
 * Will call base_island::join() (the default implementation) internally. No other side effects.
 */
base_island::~base_island()
{
	base_island::join();
}

/// Return a string identifying the island's type.
/**
 * Default implementation will return the island's C++ mangled name.
 *
 * @return a string identifying the island's type.
 */
std::string base_island::get_name() const
{
	join();
	return typeid(*this).name();
}

/// Return terse human readable representation of the island.
/**
 * Will return a formatted string containing:
 * - island type from base_island::get_name(),
 * - description of the algorithm,
 * - the output of population::human_readable_terse().
 *
 * @return string containing terse human readable representation of the island.
 */
std::string base_island::human_readable_terse() const
{
	join();
	std::ostringstream oss;
	oss << "Island type: " << get_name() << "\n\n";
	oss << *m_algo << "\n\n";
	oss << "Evolution time: " << m_evo_time << " milliseconds\n\n";
	oss << *m_s_policy << '\n';
	oss << *m_r_policy << '\n';
	oss << m_pop.human_readable_terse() << '\n';
	return oss.str();
}

/// Return human readable representation of the island.
/**
 * Will return a formatted string containing:
 * - island type from base_island::get_name(),
 * - description of the algorithm,
 * - the output of population::human_readable().
 *
 * @return string containing complete human readable representation of the island.
 */
std::string base_island::human_readable() const
{
	join();
	std::ostringstream oss;
	oss << "Island type: " << get_name() << "\n\n";
	oss << *m_algo << "\n\n";
	oss << "Evolution time: " << m_evo_time << " milliseconds\n\n";
	oss << *m_s_policy << '\n';
	oss << *m_r_policy << '\n';
	oss << m_pop.human_readable();
	return oss.str();
}

/// Join island.
/**
 * This method is intended to block the flow of the program until any ongoing evolution has terminated.
 * The default implementation will join on the internal thread object if an evolution is ongoing,
 * otherwise it will be a no-op.
 * Re-implementation of this method should always call the default implementation.
 */
void base_island::join() const
{
	if (m_evo_thread && m_evo_thread->joinable()) {
		m_evo_thread->join();
	}
}

/// Thread entry hook.
/**
 * This method will be called before any other operation takes place in the threads spawned during
 * evolution. Default implementation is a no-op.
 */
void base_island::thread_entry()
{}

/// Thread exit hook.
/**
 * This method will be called after any other operation has taken place in the threads spawned during
 * evolution. Default implementation is a no-op.
 */
void base_island::thread_exit()
{}

// RAII class to call thread hooks in base_island.
struct base_island::raii_thread_hook
{
	raii_thread_hook(base_island *ptr):m_ptr(ptr)
	{
		m_ptr->thread_entry();
	}
	~raii_thread_hook()
	{
		m_ptr->thread_exit();
	}
	base_island *m_ptr;
};

// Evolver thread object. This is a callable helper object used to launch an evolution for a given number of iterations.
struct base_island::int_evolver {
	int_evolver(base_island *i, const std::size_t &n):m_i(i),m_n(n) {}
	void operator()();
	void juice_impl(boost::posix_time::ptime &);
	base_island 		*m_i;
	const std::size_t	m_n;
};

void base_island::int_evolver::juice_impl(boost::posix_time::ptime &start)
{
	start = boost::posix_time::microsec_clock::local_time();
	// Synchronise start with all other threads if we are in an archi.
	if (m_i->m_archi) {
		m_i->m_archi->sync_island_start();
	}
	const raii_thread_hook hook(m_i);
	for (std::size_t i = 0; i < m_n; ++i) {
		// Call pre-evolve hooks.
		if (m_i->m_archi) {
			m_i->m_archi->pre_evolution(*m_i);
		}
		m_i->m_pop.problem().pre_evolution(m_i->m_pop);
		// Call the evolution.
		m_i->perform_evolution(*m_i->m_algo,m_i->m_pop);
		// Post-evolve hooks.
		if (m_i->m_archi) {
			m_i->m_archi->post_evolution(*m_i);
		}
		m_i->m_pop.problem().post_evolution(m_i->m_pop);
		// Set the interruption point.
		boost::this_thread::interruption_point();
	}
}

// Implementation of int_evolver's juice.
void base_island::int_evolver::operator()()
{
	boost::posix_time::ptime start;
	try {
		juice_impl(start);
	} catch (const boost::thread_interrupted &) {
		// In case of interruption, don't do anything special.
	} catch (const std::exception &e) {
		std::cout << "Error during island evolution using " << m_i->m_algo->get_name() << ": " << e.what() << std::endl;
	} catch (...) {
		std::cout << "Error during island evolution using " << m_i->m_algo->get_name() << ", unknown exception caught. :(" << std::endl;
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

/// Evolve island n times.
/**
 * Call the internal algorithm's algorithm::base::evolve() method n times on the internal population, using an island-specific
 * mechanism for the actual execution of the code.
 *
 * During evolution, the island is locked down and no actions on it are possible,
 * but the flow of the rest of the program might continue without waiting for all evolutions to finish. To explicitly block the program until all evolution runs
 * have been performed on the island, call the join() method.
 *
 * @param[in] n number of algorithm::base::evolve() calls that will be performed by the internal algorithm on the population.
 */
void base_island::evolve(int n)
{
	join();
	const std::size_t n_evo = boost::numeric_cast<std::size_t>(n);
	try {
		m_evo_thread.reset(new boost::thread(int_evolver(this,n_evo)));
	} catch (...) {
		pagmo_throw(std::runtime_error,"failed to launch the thread");
	}
}

// Time-dependent evolver thread object. This is a callable helper object used to launch an evolution for a specified amount of time.
struct base_island::t_evolver {
	t_evolver(base_island *i, const std::size_t &t):m_i(i),m_t(t) {}
	void operator()();
	void juice_impl(boost::posix_time::ptime &);
	base_island 		*m_i;
	const std::size_t	m_t;
};

void base_island::t_evolver::juice_impl(boost::posix_time::ptime &start)
{
	boost::posix_time::time_duration diff;
	start = boost::posix_time::microsec_clock::local_time();
	// Synchronise start.
	if (m_i->m_archi) {
		m_i->m_archi->sync_island_start();
	}
	const raii_thread_hook hook(m_i);
	do {
		if (m_i->m_archi) {
			m_i->m_archi->pre_evolution(*m_i);
		}
		m_i->m_pop.problem().pre_evolution(m_i->m_pop);
		m_i->perform_evolution(*m_i->m_algo,m_i->m_pop);
		if (m_i->m_archi) {
			m_i->m_archi->post_evolution(*m_i);
		}
		m_i->m_pop.problem().post_evolution(m_i->m_pop);
		// Set the interruption point.
		boost::this_thread::interruption_point();
		diff = boost::posix_time::microsec_clock::local_time() - start;
		// Take care of negative timings.
	} while (diff.total_milliseconds() < 0 || boost::numeric_cast<std::size_t>(diff.total_milliseconds()) < m_t);
}

// Perform at least one evolution, and continue evolving until at least a certain amount of time has passed.
void base_island::t_evolver::operator()()
{
	boost::posix_time::ptime start;
	try {
		juice_impl(start);
	} catch (const boost::thread_interrupted &) {
		// In case of interruption, don't do anything special.
	} catch (const std::exception &e) {
		std::cout << "Error during island evolution using " << m_i->m_algo->get_name() << ": " << e.what() << std::endl;
	} catch (...) {
		std::cout << "Error during island evolution using " << m_i->m_algo->get_name() << ", unknown exception caught. :(" << std::endl;
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
 * Call the internal algorithm's algorithm::base::evolve() method on the population at least once, and keep calling it until at least t milliseconds
 * (in "wall clock" time) have elapsed. Will fail if t is negative.
 *
 * During evolution, the island is locked down and no actions on it are possible,
 * but the flow of the rest of the program might continue without waiting for all evolutions to finish. To explicitly block the program until all evolution runs
 * have been performed on the island, call the join() method.
 *
 * @param[in] t minimum evolution time in milliseconds.
 */
void base_island::evolve_t(int t)
{
	join();
	const std::size_t t_evo = boost::numeric_cast<std::size_t>(t);
	try {
		m_evo_thread.reset(new boost::thread(t_evolver(this,t_evo)));
	} catch (...) {
		pagmo_throw(std::runtime_error,"failed to launch the thread");
	}
}

/// Interrupt evolution.
/**
 * If an evolution is undergoing, the evolution will be stopped the first time the flow reaches one of the internal interruption points.
 * The method will block until the interruption point has been reached.
 */
void base_island::interrupt()
{
	if (m_evo_thread) {
		m_evo_thread->interrupt();
		join();
	}
}

/// Query the status of the island.
/**
 * @return true if the island is evolving, false otherwise.
 */
bool base_island::busy() const
{
	if (!m_evo_thread) {
		return false;
	}
	return (!m_evo_thread->timed_join(boost::posix_time::milliseconds(1)));
}

/// Return the total evolution time in milliseconds.
/**
 * Note that on many 32-bit machines this counter will wrap after roughly 49 days. On most 64-bit machines, the wrapping time
 * will be around 584 million years.
 *
 * @return number of milliseconds spent evolving the population.
 */
std::size_t base_island::get_evolution_time() const
{
	join();
	return m_evo_time;
}

/// Return copy of the internal algorithm.
/**
 * @return algorithm::base_ptr to the cloned algorithm.
 */
algorithm::base_ptr base_island::get_algorithm() const
{
	join();
	return m_algo->clone();
}

/// Sets a decision vector
/**
 * Assigns a decision vector to the i-th individual of the island population
 *
 * @param[in] i individual
 * @param[in] x decision_vector
 */
void base_island::set_x(population::size_type i, const decision_vector &x)
{
	join();
	m_pop.set_x(i,x);
}

/// Sets a velocity vector
/**
 * Assigns a velocity vector to the i-th individual of the island population
 *
 * @param[in] i individual
 * @param[in] v velocity_vector
 */
void base_island::set_v(population::size_type i, const decision_vector &v)
{
	join();
	m_pop.set_v(i,v);
}

/// Algorithm setter.
/**
 * The input algorithm will be cloned.
 *
 * @param[in] algo new algorithm::base for the island.
 */
void base_island::set_algorithm(const algorithm::base &algo)
{
	join();
	m_algo = algo.clone();
}

/// Return copy of the internal problem.
/**
 * @return problem::base_ptr to the cloned problem.
 */
problem::base_ptr base_island::get_problem() const
{
	join();
	return m_pop.problem().clone();
}

/// Size of the internal population.
/**
 * @return size of the internal population.
 */
population::size_type base_island::get_size() const
{
	join();
	return m_pop.size();
}

/// Get a copy of the selection policy.
/**
 * @return migration::base_s_policy_ptr to the cloned policy.
 */
migration::base_s_policy_ptr base_island::get_s_policy() const
{
	join();
	return m_s_policy->clone();
}

/// Get a copy of the replacement policy.
/**
 * @return migration::base_r_policy_ptr to the cloned policy.
 */
migration::base_r_policy_ptr base_island::get_r_policy() const
{
	join();
	return m_r_policy->clone();
}

/// Get a copy of the internal population.
/**
 * @return copy of the population contained in the island.
 */
population base_island::get_population() const
{
	join();
	return m_pop;
}

/// Set internal population.
/**
 * @param[in] pop to be copied into the island.
 */
void base_island::set_population(const population &pop)
{
	join();
	m_pop = pop;
}

struct unary_predicate {
	unary_predicate(std::pair<population::size_type, archipelago::size_type> pair) : m_pair(pair) {};
	bool operator()(std::pair<population::size_type, archipelago::size_type> x){
		return (m_pair.second == x.second);
	}
	std::pair<population::size_type, archipelago::size_type> m_pair;
	
};
// Accept individuals incoming from a migration operation. Returns an std::vector containing the number of accepted individuals from each island
std::vector<std::pair<population::size_type, archipelago::size_type> > base_island::accept_immigrants(std::vector<std::pair<population::size_type, population::individual_type> > &immigrant_pairs)
{
	std::vector<std::pair<population::size_type, archipelago::size_type> > retval;
	// Make sure we are in an archipelago.
	pagmo_assert(m_archi);
	// We shuffle the immigrants as to make sure not to give preference to a particular island
	boost::uniform_int<int> pop_idx(0,immigrant_pairs.size());
	boost::variate_generator<boost::mt19937 &, boost::uniform_int<int> > p_idx(m_pop.m_urng,pop_idx);
	std::random_shuffle(immigrant_pairs.begin(),immigrant_pairs.end(), p_idx);
	// We extract the immigrants from the pair
	std::vector<population::individual_type> immigrants;
	immigrants.reserve(immigrant_pairs.size());
	for (size_t i=0;i<immigrant_pairs.size();++i) {
		immigrants.push_back(immigrant_pairs[i].second);
	}
	
	std::vector<std::pair<population::size_type,std::vector<population::individual_type>::size_type> > rep;
	rep = m_r_policy->select(immigrants,m_pop);
	for (std::vector<std::pair<population::size_type,std::vector<population::individual_type>::size_type> >::const_iterator
		rep_it = rep.begin(); rep_it != rep.end(); ++rep_it)
	{
		pagmo_assert((*rep_it).first < m_pop.m_container.size() && (*rep_it).second < immigrants.size());
		m_pop.m_container[(*rep_it).first] = immigrants[(*rep_it).second];
		m_pop.update_champion((*rep_it).first);
		m_pop.update_dom((*rep_it).first);
		std::pair<population::size_type, archipelago::size_type> pair = std::make_pair(1.0, immigrant_pairs[(*rep_it).second].first);
		std::vector<std::pair<population::size_type, archipelago::size_type> >::iterator where;
		where = std::find_if(retval.begin(), retval.end(), unary_predicate(pair));
		if (where == retval.end()) {
			retval.push_back(pair);
		}
		else {
			(*where).first++;
		}
	}
	return(retval);
}

// Get individuals migrating from here.
std::vector<population::individual_type> base_island::get_emigrants() 
{
	return m_s_policy->select(m_pop);
}

/// Overload stream operator for pagmo::base_island.
/**
 * Equivalent to printing base_island::human_readable() to stream.
 *
 * @param[in] s stream to which the island will be sent.
 * @param[in] isl island to be sent to stream.
 *
 * @return reference to s.
 */
std::ostream &operator<<(std::ostream &s, const base_island &isl)
{
	s << isl.human_readable();
	return s;
}

}
