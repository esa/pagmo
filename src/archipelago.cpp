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

#include <algorithm>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/thread/barrier.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "archipelago.h"
#include "algorithm/base.h"
#include "base_island.h"
#include "exceptions.h"
#include "island.h"
#include "population.h"
#include "problem/base.h"
#include "rng.h"
#include "topology/base.h"
#include "topology/unconnected.h"

namespace pagmo {

// Check we are not using bogus values for the enums.
void archipelago::check_migr_attributes() const
{
	if (m_dist_type < 0 || m_dist_type > 1 ||m_migr_dir < 0 || m_migr_dir > 1) {
		pagmo_throw(value_error,"invalid value for migration attribute (destination and/or direction)");
	}
}

/// Default constructor.
/**
 * Will construct an empty archipelago with topology::unconnected topology, with point_to_point distribution_type and destination migration_direction.
 *
 * @param[in] dt distribution type.
 * @param[in] md migration direction.
 */
archipelago::archipelago(distribution_type dt, migration_direction md):m_islands_sync_point(),m_topology(new topology::unconnected()),
	m_dist_type(dt),m_migr_dir(md),
	m_migr_map(),m_drng(rng_generator::get<rng_double>()),m_urng(rng_generator::get<rng_uint32>()),m_migr_mutex()
{
	check_migr_attributes();
}

/// Constructor from topology.
/**
 * Will construct an empty archipelago with provided topology (which will be deep-copied internally), with point_to_point distribution_type and
 * destination migration_direction.
 *
 * @param[in] t topology that will be associated to the archipelago.
 * @param[in] dt distribution type.
 * @param[in] md migration direction.
 */
archipelago::archipelago(const topology::base &t, distribution_type dt, migration_direction md):
	m_islands_sync_point(),m_topology(),m_dist_type(dt),m_migr_dir(md),
	m_migr_map(),m_drng(rng_generator::get<rng_double>()),m_urng(rng_generator::get<rng_uint32>()),m_migr_mutex()
{
	// NOTE: we cannot set the topology in the initialiser list directly,
	// since we do not know if the topology is suitable. Set it here.
	check_migr_attributes();
	set_topology(t);
}

/// Constructor from problem, algorithm, archipelago size, island sizes, topology and migration attributes.
/**
 * Constructs n islands of m individuals each, with assigned problem p and algorithm a, and inserts them with push_back() into the archipelago,
 * whose topology is set to t, with point_to_point distribution_type and destination migration_direction.
 *
 * @param[in] p problem which will be assigned to all islands.
 * @param[in] a algorithm which will be assigned to all islands.
 * @param[in] n number of islands.
 * @param[in] m number of individuals on each island.
 * @param[in] t topology.
 * @param[in] dt distribution type.
 * @param[in] md migration direction.
 */
archipelago::archipelago(const problem::base &p, const algorithm::base &a, int n, int m, const topology::base &t, distribution_type dt, migration_direction md):
	m_islands_sync_point(),m_topology(new topology::unconnected()),m_dist_type(dt),m_migr_dir(md),
	m_migr_map(),m_drng(rng_generator::get<rng_double>()),m_urng(rng_generator::get<rng_uint32>()),m_migr_mutex()
{
	check_migr_attributes();
	for (size_type i = 0; i < boost::numeric_cast<size_type>(n); ++i) {
		push_back(island(p,a,m));
	}
	// Set topology after pushing back, so that it is possible to give an already-built topology to the constructor
	// without everything blowing up.
	set_topology(t);
}

/// Copy constructor.
/**
 * Will synchronise input archipelago before deep-copying all its elements.
 *
 * @param[in] a archipelago to be copied.
 */
archipelago::archipelago(const archipelago &a)
{
	a.join();
	m_container = a.m_container;
	m_topology = a.m_topology->clone();
	m_dist_type = a.m_dist_type;
	m_migr_dir = a.m_migr_dir;
	m_migr_map = a.m_migr_map;
	m_drng = a.m_drng;
	m_urng = a.m_urng;
	m_migr_hist = a.m_migr_hist;
}

/// Assignment operator.
/**
 * Will synchronise this and the input archipelago before deep-copying all the elements into this.
 *
 * @param[in] a archipelago used for assignment.
 *
 * @return reference to this.
 */
archipelago &archipelago::operator=(const archipelago &a)
{
	if (this != &a) {
		join();
		a.join();
		m_container = a.m_container;
		m_topology = a.m_topology->clone();
		m_dist_type = a.m_dist_type;
		m_migr_dir = a.m_migr_dir;
		m_migr_map = a.m_migr_map;
		m_drng = a.m_drng;
		m_urng = a.m_urng;
		m_migr_hist = a.m_migr_hist;
	}
	return *this;
}

/// Destructor.
/**
 * Will call join() before returning. No other side effects.
 */
archipelago::~archipelago()
{
	join();
}

/// Wait until evolution on each island has terminated.
/**
 * Will call iteratively island::join() on all islands of the archipelago.
 */
void archipelago::join() const
{
	const const_iterator it_f = m_container.end();
	for (const_iterator it = m_container.begin(); it != it_f; ++it) {
		(*it)->join();
	}
}

archipelago::size_type archipelago::locate_island(const base_island &isl) const
{
	// TODO: iterate and increase a size type instead of using distance(), which returns a signed int.
	// Also, assert the island is actually found here?
	const_iterator it = m_container.begin();
	for (; it != m_container.end(); ++it) {
		if (&(*(*it)) == &isl) {
			break;
		}
	}
	return std::distance(m_container.begin(),it);
}

/// Add an island to the archipelago.
/**
 * Both the island and the archipelago will be synchronised before any operation takes place. The island will then
 * be appended to the archipelago and connected to the existing islands via topology::base::push_back().
 *
 * Will fail if check_island() returns false on isl.
 *
 * @param[in] isl island to be added to the archipelago.
 */
void archipelago::push_back(const base_island &isl)
{
	// NOTE: the joins are already done in check_island().
	if (!check_island(isl)) {
		pagmo_throw(value_error,"cannot push_back() incompatible island");
	}
	m_container.push_back(isl.clone());
	// Make sure the clone() method is implemented properly in the problem of the island we just inserted. This is important
	// beacuse otherwise we get inconsistent behaviour in migration: we could be inserting a problem which is not really compatible.
	// NOTE: if we ever implement automatic checking on clone(), this can go away.
// 	if (m_container.back().m_pop.problem() != isl.m_pop.problem()) {
		// Remove the island we just inserted.
// 		m_container.pop_back();
// 		pagmo_throw(value_error,"the problem's clone() method implementation does not produce a problem equal to itself: please correct it");
// 	}
	// Tell the island that it is living in an archipelago now.
	m_container.back()->m_archi = this;
	// Insert the island in the topology.
	m_topology->push_back();
}

/// Set island algorithm.
/**
 * Set algorithm of island number idx to a.
 *
 * @param[in] idx island index.
 * @param[in] a algorithm to be set.
 *
 * @throws pagmo::index_error if index is not smaller than archipelago size.
 */
void archipelago::set_algorithm(const size_type &idx, const algorithm::base &a)
{
	join();
	if (idx >= m_container.size()) {
		pagmo_throw(index_error,"invalid island index");
	}
	m_container[idx]->set_algorithm(a);
}

/// Get the size of the archipelago.
/**
 * @return the number of islands contained in the archipelago.
 */
archipelago::size_type archipelago::get_size() const
{
	join();
	return m_container.size();
}

/// Return human readable representation of the archipelago.
/**
 * Will return a formatted string containing:
 * - the number of islands,
 * - the output of topology::base::human_readable_terse() for the topology,
 * - the output of island::human_readable_terse() for each island.
 *
 * @return formatted string containing the human readable representation of the archipelago.
 */
std::string archipelago::human_readable() const
{
	join();
	std::ostringstream oss;
	oss << "-= Archipelago =-\n\n";
	oss << "Number of islands:\t" << m_container.size() << "\n\n";
	oss << m_topology->human_readable_terse() << '\n';
	for (size_type i = 0; i < m_container.size(); ++i) {
		oss << "Island index: #" << i << "\n\n";
		oss << m_container[i]->human_readable_terse() << '\n';
	}
	return oss.str();
}

/// Return a copy of the topology.
/**
 * @return a topology::base_ptr to a clone of the internal topology.
 */
topology::base_ptr archipelago::get_topology() const
{
	join();
	return m_topology->clone();
}

/// Set topology.
/**
 * A valid topology must contain only island indices smaller than the size of the archipelago. If this condition is satisfied, then the incoming topology
 * will become the new archipelago topology (after calling push_back() a number of times necessary to fill in the island indices missing in the topology).
 * Otherwise, a value_error exception will be raised.
 *
 * @param[in] tp new topology for the archipelago.
 */
void archipelago::set_topology(const topology::base &tp)
{
	join();
	topology::base_ptr t = tp.clone();
	if (m_container.size() < boost::numeric_cast<size_type>(t->get_number_of_vertices())) {
		pagmo_throw(value_error,"invalid topology, too many vertices");
	}
	// Push back the missing vertices, if any.
	for (size_type i = boost::numeric_cast<size_type>(t->get_number_of_vertices());
		i < m_container.size(); ++i)
	{
		t->push_back();
	}
	// The topology is ok, assign it.
	m_topology = t;
}

/// Check whether an island is compatible with the archipelago.
/**
 * Will return true if any of these conditions holds:
 * - archipelago is empty,
 * - the problem of the incoming island is compatible with the problems of the islands already
 *   present in the archipelago (via problem::base::is_compatible()).
 *
 * Otherwise, false will be returned.
 *
 * @param[in] isl island which will be checked.
 *
 * @return true if isl is compatible with the archipelago, false otherwise.
 */
bool archipelago::check_island(const base_island &isl) const
{
	join();
	isl.join();
	// We need to perform checks only if there are other islands in the archipelago.
	// Otherwise, any island will be allowed in.
	if (m_container.size() && (!isl.m_pop.problem().is_compatible(m_container[0]->m_pop.problem()))) {
		return false;
	}
	return true;
}

// Reset island synchronisation barrier. This method is intended as a shortcut,
// do NOT use it if the archipelago has not been joined!
void archipelago::reset_barrier(const size_type &size)
{
	if (size) {
		m_islands_sync_point.reset(new boost::barrier(boost::numeric_cast<unsigned int>(size)));
	}
}

// Helper function to insert a list of candidates immigrants into an immigrants vector, given the source and destination island.
void archipelago::build_immigrants_vector(std::vector<individual_type> &immigrants, const base_island &src_isl,
	base_island &dest_isl, const std::vector<individual_type> &candidates, migr_hist_type &h) const
{
	if (dest_isl.m_pop.problem() == src_isl.m_pop.problem()) {
		// If the problem of the originating island is identical to that of destination island, then we can just
		// append the immigrants.
		immigrants.insert(immigrants.end(),candidates.begin(),candidates.end());
	} else {
		// If the problems are not identical, then we need to re-evaluate the incoming individuals according
		// to dest_isl's problem.
		individual_type tmp;
		tmp.cur_v.resize(dest_isl.m_pop.problem().get_dimension());
		tmp.cur_f.resize(dest_isl.m_pop.problem().get_f_dimension());
		tmp.cur_c.resize(dest_isl.m_pop.problem().get_c_dimension());
		for (std::vector<individual_type>::const_iterator ind_it = candidates.begin();
			ind_it != candidates.end(); ++ind_it)
		{
			// Skip individual if it is not within the bounds of the problem
			// in the destination island.
			if (!dest_isl.m_pop.problem().verify_x(ind_it->cur_x)) {
				continue;
			}
			tmp.cur_x = ind_it->cur_x;
			dest_isl.m_pop.problem().objfun(tmp.cur_f,tmp.cur_x);
			dest_isl.m_pop.problem().compute_constraints(tmp.cur_c,tmp.cur_x);
			// Set the best properties to the current ones.
			tmp.best_x = tmp.cur_x;
			tmp.best_f = tmp.cur_f;
			tmp.best_c = tmp.cur_c;
			immigrants.push_back(tmp);
		}
	}
	// Record the migration history.
	h.push_back(boost::make_tuple(
		boost::numeric_cast<population::size_type>(immigrants.size()),
		locate_island(src_isl),
		locate_island(dest_isl)
	));
}

// This method will be called by each island of the archipelago before starting evolution. Its task is
// to select from the other islands, according to the topology and the migration/distribution type and direction,
// the individuals that will migrate into the island.
void archipelago::pre_evolution(base_island &isl)
{
	// Lock the migration logic.
	lock_type lock(m_migr_mutex);
	// Make sure the island belongs to the archipelago.
	pagmo_assert(isl.m_archi == this);
	// Make sure the archipelago is not empty.
	pagmo_assert(m_container.size());
	// Determine the island's index in the archipelago.
	const size_type isl_idx = locate_island(isl);
	pagmo_assert(isl_idx < m_container.size());
	// Migration history for this migration.
	migr_hist_type h;
	//1. Obtain immigrants.
	std::vector<individual_type> immigrants;
	switch (m_migr_dir) {
		case source:
			// For source migration direction, migration map contains islands' "inboxes". Or, in other words, it contains
			// the individuals that are destined to go into the islands. Such inboxes have been assembled previously,
			// during a post_evolution operation.
			// Iterate over all the vectors of individuals provided by the different islands.
			for (boost::unordered_map<size_type,std::vector<individual_type> >::iterator it = m_migr_map[isl_idx].begin();
				it != m_migr_map[isl_idx].end(); ++it)
			{
				pagmo_assert(it->first < m_container.size());
				build_immigrants_vector(immigrants,*m_container[it->first],isl,it->second,h);
			}
			// Delete stuff in the migration map.
			m_migr_map.erase(isl_idx);
			break;
		case destination:
			// For destination migration direction, items in the migration map behave like "outboxes", i.e. each one is a
			// "database of best individuals" seen in the islands of the archipelago.
			// Get neighbours connecting into isl.
			const std::vector<topology::base::vertices_size_type> inv_adj_islands(m_topology->get_v_inv_adjacent_vertices(boost::numeric_cast<topology::base::vertices_size_type>(isl_idx)));
			// Do something only if there are adjacent islands.
			if (inv_adj_islands.size()) {
				switch (m_dist_type) {
					case point_to_point:
					{
						// Get the index of a random island connecting into isl.
						boost::uniform_int<std::vector<topology::base::vertices_size_type>::size_type> u_int(0,inv_adj_islands.size() - 1);
						const size_type rn_isl_idx = boost::numeric_cast<size_type>(inv_adj_islands[u_int(m_urng)]);
						// Get the immigrants from the outbox of the random island. Note the redundant information in the last
						// argument of the function.
						pagmo_assert(m_migr_map[rn_isl_idx].size() <= 1);
						build_immigrants_vector(immigrants,*m_container[rn_isl_idx],isl,m_migr_map[rn_isl_idx][rn_isl_idx],h);
						break;
					}
					case broadcast:
						// For broadcast migration fetch immigrants from all neighbour islands' databases.
						for (std::vector<topology::base::vertices_size_type>::size_type i = 0; i < inv_adj_islands.size(); ++i) {
							const size_type src_isl_idx = boost::numeric_cast<size_type>(inv_adj_islands[i]);
							pagmo_assert(m_migr_map[src_isl_idx].size() <= 1);
							build_immigrants_vector(immigrants,*m_container[src_isl_idx],isl,m_migr_map[src_isl_idx][src_isl_idx],h);
						}
				}
			}
	}
	//2. Insert immigrants into population.
	if (immigrants.size()) {
		if (m_drng() < isl.m_migr_prob) {
			isl.accept_immigrants(immigrants);
			m_migr_hist.insert(m_migr_hist.end(),h.begin(),h.end());
		}
	}
}

// This method will be called after each island isl has completed an evolution. Its purpose is to get individuals
// emigrating from isl and put them at disposal of the other islands of the archipelago, according to the migration attributes
// and the topology.
void archipelago::post_evolution(base_island &isl)
{
	// Lock the migration logic.
	lock_type lock(m_migr_mutex);
	// Make sure the island belongs to the archipelago.
	pagmo_assert(isl.m_archi == this);
	// Make sure the archipelago is not empty.
	pagmo_assert(m_container.size());
	// Determine the island's index in the archipelago.
	const size_type isl_idx = locate_island(isl);
	pagmo_assert(isl_idx < m_container.size());
	// Create the vector of emigrants.
	std::vector<individual_type> emigrants;
	switch (m_migr_dir) {
		case source:
		{
			// Get the islands to which isl connects.
			const std::vector<topology::base::vertices_size_type> adj_islands(m_topology->get_v_adjacent_vertices(boost::numeric_cast<topology::base::vertices_size_type>(isl_idx)));
			if (adj_islands.size()) {
				emigrants = isl.get_emigrants();
				// Do something only if we have emigrants.
				if (emigrants.size()) {
					switch (m_dist_type)
					{
						case point_to_point:
						{
							// For one-to-one migration choose a random neighbour island and put immigrants to its inbox.
							boost::uniform_int<std::vector<topology::base::vertices_size_type>::size_type> u_int(0,adj_islands.size() - 1);
							const size_type chosen_adj = boost::numeric_cast<size_type>(adj_islands[u_int(m_urng)]);
							m_migr_map[chosen_adj][isl_idx].insert(m_migr_map[chosen_adj][isl_idx].end(),emigrants.begin(),emigrants.end());
							break;
						}
						case broadcast:
							// For broadcast migration put immigrants to all neighbour islands' inboxes.
							for (std::vector<topology::base::vertices_size_type>::size_type i = 0; i < adj_islands.size(); ++i) {
								m_migr_map[boost::numeric_cast<size_type>(adj_islands[i])][isl_idx]
								.insert(m_migr_map[boost::numeric_cast<size_type>(adj_islands[i])][isl_idx].end(),
								emigrants.begin(),emigrants.end());
							}
					}
				}
			}
			break;
		}
		case destination:
			// For destination migration direction, migration map behaves like "outboxes", i.e. each is a "database of best individuals" for corresponding island.
			emigrants = isl.get_emigrants();
			pagmo_assert(m_migr_map[isl_idx].size() <= 1);
			m_migr_map[isl_idx][isl_idx].swap(emigrants);
	}
}

// Functor to count the number of blocking islands.
struct archipelago::count_if_blocking
{
	bool operator()(const base_island_ptr &ptr) const
	{
		return ptr->is_blocking();
	}
};

// Implementation of blocking attribute.
bool archipelago::is_blocking_impl() const
{
	// Make sure we do not overflow.
	boost::numeric_cast<std::iterator_traits<const_iterator>::difference_type>(m_container.size());
	return std::count_if(m_container.begin(),m_container.end(),count_if_blocking());
}

/// Archipelago's blocking attribute.
/**
 * @return true if at least one island's base_island::is_blocking() method returns true.
 */
bool archipelago::is_blocking() const
{
	join();
	return is_blocking_impl();
}

// Functor to count the number of thread-safe islands.
struct archipelago::count_if_thread_safe
{
	bool operator()(const base_island_ptr &ptr) const
	{
		return ptr->is_thread_safe();
	}
};

// Implementation of thread-safe attribute.
bool archipelago::is_thread_safe_impl() const
{
	// Make sure we do not overflow.
	boost::numeric_cast<std::iterator_traits<const_iterator>::difference_type>(m_container.size());
	return std::count_if(m_container.begin(),m_container.end(),count_if_thread_safe());
}

/// Archipelago's thread safety attribute.
/**
 * @return true if at least one island's base_island::is_thread_safe() method returns true.
 */
bool archipelago::is_thread_safe() const
{
	join();
	return is_thread_safe();
}

/// Run the evolution for the given number of iterations.
/**
 * Will iteratively call island::evolve(n) on each island of the archipelago and then return.
 *
 * \param[in] n number of time each island will be evolved.
 */
void archipelago::evolve(int n)
{
	join();
	const iterator it_f = m_container.end();
	if (is_thread_safe_impl()) {
		// Reset thread barrier.
		reset_barrier(m_container.size());
		for (iterator it = m_container.begin(); it != it_f; ++it) {
			(*it)->evolve(n);
		}
	} else {
		// Build a vector of iterators to the islands.
		std::vector<iterator> it_vector;
		for (iterator it = m_container.begin(); it != it_f; ++it) {
			it_vector.push_back(it);
		}
		for (int i = 0; i < n; ++i) {
			// Shuffle the vector of island iterators to simulate async operations.
			std::random_shuffle(it_vector.begin(),it_vector.end());
			for (std::vector<iterator>::iterator it = it_vector.begin(); it != it_vector.end(); ++it) {
				(*(*it))->evolve();
			}
		}
	}
	// If the archipelago is blocking, wait for evolutions to finish.
	if (is_blocking_impl()) {
		join();
	}
}

// Helper function to determine if all islands evolved for at least t milliseconds.
template <class Vector>
static bool all_islands_t_evolved(const Vector &v, int t)
{
	for (typename Vector::const_iterator it = v.begin(); it != v.end(); ++it) {
		if (it->second < t) {
			return false;
		}
	}
	return true;
}

/// Run the evolution for a minimum amount of time.
/**
 * Will iteratively call island::evolve_t(n) on each island of the archipelago and then return.
 *
 * \param[in] t amount of time to evolve each island (in milliseconds).
 */
void archipelago::evolve_t(int t)
{
	join();
	const iterator it_f = m_container.end();
	if (is_thread_safe_impl()) {
		reset_barrier(m_container.size());
		for (iterator it = m_container.begin(); it != it_f; ++it) {
			(*it)->evolve_t(t);
		}
	} else {
		// Build a vector of (iterators,evolution time) pairs.
		std::vector<std::pair<iterator,int> > it_t_vector;
		for (iterator it = m_container.begin(); it != m_container.end(); ++it) {
			it_t_vector.push_back(std::make_pair(it,0));
		}
		while (!all_islands_t_evolved(it_t_vector,t)) {
			// Shuffle the vector to simulate async operations.
			std::random_shuffle(it_t_vector.begin(),it_t_vector.end());
			for (std::vector<std::pair<iterator,int> >::iterator it = it_t_vector.begin(); it != it_t_vector.end(); ++it) {
				// Evolve only if island has not evolved already for the desired time.
				if (it->second < t) {
					// Record the initial evolution time for the island.
					const std::size_t initial_time = (*it->first)->m_evo_time;
					(*it->first)->evolve();
					// Here we must be careful. In "normal" conditions everything is fine and dandy and evo_time will now be higher than
					// intial time. However, if the clock is screwed, if the counter is wrapping past the numerical limit or the evolution
					// lasted 0 milliseconds, etc. then we must detect and fix this. Policy: add one second to the total evolution time for
					// the island, so that at least we are sure we don't end up in an endless cycle.
					const std::size_t time_diff = ((*it->first)->m_evo_time > initial_time) ? ((*it->first)->m_evo_time - initial_time) : 1000;
					it->second += boost::numeric_cast<int>(time_diff);
				}
			}
		}
	}
	// If the archipelago is blocking, wait for evolutions to finish.
	if (is_blocking_impl()) {
		join();
	}
}

/// Query the status of the archipelago.
/**
 * @return true if at least one island is evolving, false otherwise.
 */
bool archipelago::busy() const
{
	const const_iterator it_f = m_container.end();
	for (const_iterator it = m_container.begin(); it != it_f; ++it) {
		if ((*it)->busy()) {
			return true;
		}
	}
	return false;
}

/// Interrupt ongoing evolution.
/**
 * Will iteratively called island::interrupt() on all the islands of the archipelago.
 */
void archipelago::interrupt()
{
	const iterator it_f = m_container.end();
	for (iterator it = m_container.begin(); it != it_f; ++it) {
		(*it)->interrupt();
	}
}

/// Island getter.
/**
 * @return a copy of the island at position idx.
 *
 * @throw value_error if idx is not less than the size of the archipelago.
 */
base_island_ptr archipelago::get_island(const size_type &idx) const
{
	join();
	if (idx >= m_container.size()) {
		pagmo_throw(index_error,"invalid island index");
	}
	base_island_ptr retval = m_container[idx]->clone();
	// The island is no more in an archipelago.
	retval->m_archi = 0;
	return retval;
}

// Synchronise the start of evolution in each island so that all threads are created and initialised
// before actually doing any computation.
void archipelago::sync_island_start() const
{
	m_islands_sync_point->wait();
}

// Dump migration history.
std::string archipelago::dump_migr_history() const
{
	join();
	std::ostringstream oss;
	for (migr_hist_type::const_iterator it = m_migr_hist.begin(); it != m_migr_hist.end(); ++it) {
		oss << *it << '\n';
	}
	return oss.str();
}

// Clear migration history.
void archipelago::clear_migr_history()
{
	join();
	m_migr_hist.clear();
}

/// Overload stream operator for pagmo::archipelago.
/**
 * Equivalent to printing archipelago::human_readable() to stream.
 *
 * @param[in] s stream to which the archipelago will be sent.
 * @param[in] a archipelago to be sent to stream.
 *
 * @return reference to s.
 */
std::ostream &operator<<(std::ostream &s, const archipelago &a)
{
	s << a.human_readable();
	return s;
}

}
