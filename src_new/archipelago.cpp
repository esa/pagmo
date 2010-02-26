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
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "archipelago.h"
#include "algorithm/base.h"
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
archipelago::archipelago(distribution_type dt, migration_direction md):m_island_sync_point(),m_topology(new topology::unconnected()),
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
	m_island_sync_point(),m_topology(t.clone()),m_dist_type(dt),m_migr_dir(md),
	m_migr_map(),m_drng(rng_generator::get<rng_double>()),m_urng(rng_generator::get<rng_uint32>()),m_migr_mutex()
{
	check_migr_attributes();
}

/// Constructor from problem, algorithm, archipelago size, island sizes and topology.
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
	m_island_sync_point(),m_topology(t.clone()),m_dist_type(dt),m_migr_dir(md),
	m_migr_map(),m_drng(rng_generator::get<rng_double>()),m_urng(rng_generator::get<rng_uint32>()),m_migr_mutex()
{
	check_migr_attributes();
	for (size_type i = 0; i < boost::numeric_cast<size_type>(n); ++i) {
		push_back(island(p,a,m));
	}
}

/// Copy constructor.
/**
 * Will synchronise a before deep-copying all its elements.
 *
 * @param[in] a archipelago to be copied.
 */
archipelago::archipelago(const archipelago &a)
{
	a.join();
	m_container = a.m_container;
	m_topology = a.m_topology->clone();
	reset_barrier();
	m_dist_type = a.m_dist_type;
	m_drng = a.m_drng;
	m_urng = a.m_urng;
}

/// Assignment operator.
/**
 * Will synchronise this and a before deep-copying all elements from a into this.
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
		reset_barrier();
		m_dist_type = a.m_dist_type;
		m_migr_dir = a.m_migr_dir;
		m_migr_map = a.m_migr_map;
		m_drng = a.m_drng;
		m_urng = a.m_urng;
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
		it->join();
	}
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
void archipelago::push_back(const island &isl)
{
	// NOTE: the joins are already done in check_island().
	if (!check_island(isl)) {
		pagmo_throw(value_error,"cannot push_back() incompatible island");
	}
	m_container.push_back(isl);
	// Tell the island that it is living in an archipelago now.
	m_container.back().m_archi = this;
	// Insert the island in the topology.
	m_topology->push_back(boost::numeric_cast<int>(m_container.size() - 1));
	// Reset the barrier.
	reset_barrier();
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
	oss << "Archipelago\n";
	oss << "===========\n\n";
	oss << "Number of islands:\t" << m_container.size() << "\n\n";
	oss << m_topology->human_readable_terse() << '\n';
	for (size_type i = 0; i < m_container.size(); ++i) {
		oss << "Island index: #" << i << "\n\n";
		oss << m_container[i].human_readable_terse() << '\n';
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
 * A valid topology must contain all and only the island indices of the current archipelago. I.e., if the size of
 * the archipelago is n, then the vertices list of the topology being set must contain all and only the integers between 0 and n-1.
 *
 * If this condition is satisfied, then the incoming topology will become the new archipelago topology. Otherwise, a value_error exception
 * will be raised.
 *
 * @param[in] t new topology for the archipelago.
 */
void archipelago::set_topology(const topology::base &t)
{
	join();
	if (m_container.size() != boost::numeric_cast<size_type>(t.get_number_of_vertices())) {
		pagmo_throw(value_error,"invalid topology, wrong number of vertices");
	}
	// Get the list of vertices and order it in ascending order.
	std::vector<int> vertices_list(t.get_vertices());
	std::sort(vertices_list.begin(),vertices_list.end());
	// t is a valid topology only if all and only island indices are represented in the vertices list.
	for (std::vector<int>::size_type i = 0; i < vertices_list.size(); ++i) {
		if (i != boost::numeric_cast<std::vector<int>::size_type>(vertices_list[i])) {
			pagmo_throw(value_error,"invalid topology, missing island index");
		}
	}
	// The topology is ok, assign it.
	m_topology = t.clone();
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
bool archipelago::check_island(const island &isl) const
{
	join();
	isl.join();
	// We need to perform checks only if there are other islands in the archipelago.
	// Otherwise, any island will be allowed in.
	if (m_container.size() && (!isl.m_pop.problem().is_compatible(m_container[0].m_pop.problem()))) {
		return false;
	}
	return true;
}

// Reset island synchronisation barrier. This method is intended as a shortcut,
// do NOT use it if the archipelago has not been joined!
void archipelago::reset_barrier()
{
	if (m_container.size()) {
		m_island_sync_point.reset(new boost::barrier(boost::numeric_cast<unsigned int>(m_container.size())));
	} else {
		m_island_sync_point.reset(0);
	}
}

// Helper function to insert a list of candidates immigrants into an immigrants vector, given the source and destination island.
void archipelago::build_immigrants_vector(std::vector<individual_type> &immigrants, const island &src_isl,
	island &dest_isl, const std::vector<individual_type> &candidates) const
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
}

// This method will be called by each island of the archipelago before starting evolution. Its task is
// to select from the other islands, according to the topology and the migration/distribution type and direction,
// the individuals that will migrate into the island.
void archipelago::pre_evolution(island &isl)
{
//std::cout << "In pre\n";
	// Lock the migration logic.
	lock_type lock(m_migr_mutex);
//std::cout << "Locked pre\n";
	// Make sure the island belongs to the archipelago.
	pagmo_assert(isl.m_archi == this);
	// Make sure the archipelago is not empty.
	pagmo_assert(m_container.size());
	// Make sure the island is somewhere after the beignning of the container and within it.
	pagmo_assert(&isl >= &m_container[0] && boost::numeric_cast<size_type>(&isl - &m_container[0]) < m_container.size());
	// Determine the island's index in the archipelago.
	const size_type isl_idx = boost::numeric_cast<size_type>(&isl - &m_container[0]);
	//1. Obtain immigrants.
	std::vector<individual_type> immigrants;
	switch (m_migr_dir) {
		case source:
			// For source migration direction, migration map contains islands' "inboxes". Or, in other words, it contains
			// the individuals that are destined to go into the islands. Such inboxes have been assembled previously,
			// during a post_evolution operation.
			// TODO:verify the above.
			// Iterate over all the vectors of individuals provided by the different islands.
			for (boost::unordered_map<size_type,std::vector<individual_type> >::iterator it = m_migr_map[isl_idx].begin();
				it != m_migr_map[isl_idx].end(); ++it)
			{
				pagmo_assert(it->first < m_container.size());
				build_immigrants_vector(immigrants,m_container[it->first],isl,it->second);
			}
			// Delete stuff in the migration map.
			m_migr_map.erase(isl_idx);
			break;
		case destination:
			// For destination migration direction, items in the migration map behave like "outboxes", i.e. each one is a
			// "database of best individuals" seen in the islands of the archipelago.
			// Get neighbours connecting into isl.
			const std::vector<int> inv_adj_islands(m_topology->get_inv_adjacent_vertices(boost::numeric_cast<int>(isl_idx)));
			// Do something only if there are adjacent islands.
			if (inv_adj_islands.size()) {
				switch (m_dist_type) {
					case point_to_point:
					{
						// Get the index of a random island connecting into isl.
						boost::uniform_int<std::vector<int>::size_type> uint(0,inv_adj_islands.size() - 1);
						const size_type rn_isl_idx = boost::numeric_cast<size_type>(inv_adj_islands[uint(m_urng)]);
						// Get the immigrants from the outbox of the random island. Note the redundant information in the last
						// argument of the function.
						pagmo_assert(m_migr_map[rn_isl_idx].size() <= 1);
						build_immigrants_vector(immigrants,m_container[rn_isl_idx],isl,m_migr_map[rn_isl_idx][rn_isl_idx]);
						break;
					}
					case broadcast:
						// For broadcast migration fetch immigrants from all neighbour islands' databases.
						for (std::vector<int>::size_type i = 0; i < inv_adj_islands.size(); ++i) {
							const size_type src_isl_idx = boost::numeric_cast<size_type>(inv_adj_islands[i]);
							pagmo_assert(m_migr_map[src_isl_idx].size() <= 1);
							build_immigrants_vector(immigrants,m_container[src_isl_idx],isl,m_migr_map[src_isl_idx][src_isl_idx]);
						}
				}
			}
	}
	//2. Insert immigrants into population.
	if (immigrants.size()) {
		if (m_drng() < isl.m_migr_prob) {
			isl.accept_immigrants(immigrants);
		}
	}
//std::cout << "Done pre\n";
}

// This method will be called after each island isl has completed an evolution. Its purpose is to get individuals
// emigrating from isl and put them at disposal of the other islands of the archipelago, according to the migration attributes
// and the topology.
void archipelago::post_evolution(island &isl)
{
//std::cout << "In post\n";
	// Lock the migration logic.
	lock_type lock(m_migr_mutex);
//std::cout << "Locked post\n";
	// Make sure the island belongs to the archipelago.
	pagmo_assert(isl.m_archi == this);
	// Make sure the archipelago is not empty.
	pagmo_assert(m_container.size());
	// Make sure the island is somewhere after the beignning of the container and within it.
	pagmo_assert(&isl >= &m_container[0] && boost::numeric_cast<size_type>(&isl - &m_container[0]) < m_container.size());
	// Determine the island's index in the archipelago.
	const size_type isl_idx = boost::numeric_cast<size_type>(&isl - &m_container[0]);
	// Create the vector of emigrants.
	std::vector<individual_type> emigrants;
	switch (m_migr_dir) {
		case source:
		{
			// Get the islands to which isl connects.
			const std::vector<int> adj_islands(m_topology->get_adjacent_vertices(boost::numeric_cast<int>(isl_idx)));
			if (adj_islands.size()) {
				emigrants = isl.get_emigrants();
				// Do something only if we have emigrants.
				if (emigrants.size()) {
					switch (m_dist_type)
					{
						case point_to_point:
						{
							// For one-to-one migration choose a random neighbour island and put immigrants to its inbox.
							boost::uniform_int<std::vector<int>::size_type> uint(0,adj_islands.size() - 1);
							const size_type chosen_adj = boost::numeric_cast<size_type>(adj_islands[uint(m_urng)]);
							m_migr_map[chosen_adj][isl_idx].insert(m_migr_map[chosen_adj][isl_idx].end(),emigrants.begin(),emigrants.end());
							break;
						}
						case broadcast:
							// For broadcast migration put immigrants to all neighbour islands' inboxes.
							for (std::vector<int>::size_type i = 0; i < adj_islands.size(); ++i) {
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
			// For destination migration directio, migration map behaves like "outboxes", i.e. each is a "database of best individuals" for corresponding island.
			emigrants = isl.get_emigrants();
			pagmo_assert(m_migr_map[isl_idx].size() <= 1);
			m_migr_map[isl_idx][isl_idx].swap(emigrants);
	}
//std::cout << "Done post\n";
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
	for (iterator it = m_container.begin(); it != it_f; ++it) {
		it->evolve(n);
	}
}

/// Run the evolution for a minimum amount of time.
/**
 * Will iteratively call island::evolve_t(n) on each island of the archipelago and then return.
 *
 * \param[in] t amount of time to evolve each island (in miliseconds).
 */
void archipelago::evolve_t(int t)
{
	join();
	const iterator it_f = m_container.end();
	for (iterator it = m_container.begin(); it != it_f; ++it) {
		it->evolve_t(t);
	}
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
