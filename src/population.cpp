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

#include <algorithm>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include "problem/base.h"
#include "exceptions.h"
#include "population.h"
#include "rng.h"
#include "types.h"

namespace pagmo
{

problem::base_ptr &population_access::get_problem_ptr(population &pop)
{
	return pop.m_prob;
}

/// Constructor from problem::base and number of individuals.
/**
 * Will store a copy of the problem and will initialise the population to n randomly-generated individuals.
 * Will fail if n is negative.
 *
 * @param[in] p problem::base that will be associated to the population.
 * @param[in] n integer number of individuals in the population.
 *
 * @throw piranha::value_error if n is negative.
 */
population::population(const problem::base &p, int n):m_prob(p.clone()),m_drng(rng_generator::get<rng_double>()),m_urng(rng_generator::get<rng_uint32>())
{
	if (n < 0) {
		pagmo_throw(value_error,"number of individuals cannot be negative");
	}
	// Store sizes temporarily.
	const size_type size = boost::numeric_cast<size_type>(n);
	const fitness_vector::size_type f_size = m_prob->get_f_dimension();
	const constraint_vector::size_type c_size = m_prob->get_c_dimension();
	const decision_vector::size_type p_size = m_prob->get_dimension();
	for (size_type i = 0; i < size; ++i) {
		// Push back an empty individual.
		m_container.push_back(individual_type());
		m_dom_list.push_back(std::vector<size_type>());
		// Resize individual's elements.
		m_container.back().cur_x.resize(p_size);
		m_container.back().cur_v.resize(p_size);
		m_container.back().cur_c.resize(c_size);
		m_container.back().cur_f.resize(f_size);
		m_container.back().best_x.resize(p_size);
		m_container.back().best_c.resize(c_size);
		m_container.back().best_f.resize(f_size);
		// Initialise randomly the individual.
		reinit(i);
	}
}

// Update the domination lists.
void population::update_dom_list(const size_type &n)
{
	const size_type size = m_container.size();
	pagmo_assert(m_dom_list.size() == size && n < size);
	// Empty the domination list of individual at position n.
	m_dom_list[n].clear();
	for (size_type i = 0; i < size; ++i) {
		if (i != n) {
			// Check if individual in position i dominates individual in position n.
			if (m_prob->compare_fc(m_container[i].cur_f,m_container[i].cur_c,m_container[n].cur_f,m_container[n].cur_c)) {
				// Need to update the domination list in i. If n is already present,
				// do nothing, otherwise push_back.
				if (std::find(m_dom_list[i].begin(),m_dom_list[i].end(),n) == m_dom_list[i].end()) {
					m_dom_list[i].push_back(n);
				}
			} else {
				// We need to erase n from the domination list, if present.
				std::vector<size_type>::iterator it = std::find(m_dom_list[i].begin(),m_dom_list[i].end(),n);
				if (it != m_dom_list[i].end()) {
					m_dom_list[i].erase(it);
				}
			}
			// Check if individual in position n dominates individual in position i.
			if (m_prob->compare_fc(m_container[n].cur_f,m_container[n].cur_c,m_container[i].cur_f,m_container[i].cur_c)) {
				m_dom_list[n].push_back(i);
			}
		}
	}
}

// Init randomly the velocity of the individual in position idx.
void population::init_velocity(const size_type &idx)
{
	const decision_vector::size_type p_size = m_prob->get_dimension();
	for (decision_vector::size_type j = 0; j < p_size; ++j) {
		// Initialise velocities so that in one tick the particles travel at most half the bounds distance.
		m_container[idx].cur_v[j] = boost::uniform_real<double>(m_prob->get_lb()[j] / 2,m_prob->get_ub()[j] / 2)(m_drng);
		// Change randomly the sign of the velocity.
		m_container[idx].cur_v[j] *= (m_drng() < .5) ? 1 : -1;
	}
}

/// Re-initialise all individuals
/**
 * @see population::reinit(const size_type &).
 */
void population::reinit()
{
	for (size_type i = 0; i < size(); ++i)
	{
		reinit(i);
	}
}

/// Re-initialise individual at position idx.
/**
 * The continuous and integer parts of the chromosome will be picked randomly within the problem's bounds, the velocities
 * will be initialised randomly so that in one tick the particles travel at most half the bounds distance. Fitness and constraints
 * will be evaluated, and champion updated.
 *
 * @param[in] idx position of the individual to be re-initialised.
 *
 * @throw index_error if idx is not smaller than size().
 */
void population::reinit(const size_type &idx)
{
	if (idx >= size()) {
		pagmo_throw(index_error,"invalid index");
	}
	const decision_vector::size_type p_size = m_prob->get_dimension(), i_size = m_prob->get_i_dimension();
	// Initialise randomly the continuous part of the decision vector.
	for (decision_vector::size_type j = 0; j < p_size - i_size; ++j) {
		m_container[idx].cur_x[j] = boost::uniform_real<double>(m_prob->get_lb()[j],m_prob->get_ub()[j])(m_drng);
	}
	// Initialise randomly the integer part of the decision vector.
	for (decision_vector::size_type j = p_size - i_size; j < p_size; ++j) {
		m_container[idx].cur_x[j] = boost::uniform_int<int>(m_prob->get_lb()[j],m_prob->get_ub()[j])(m_urng);
	}
	// Initialise randomly the velocity vector.
	init_velocity(idx);
	// Fill in the constraints.
	m_prob->compute_constraints(m_container[idx].cur_c,m_container[idx].cur_x);
	// Compute the fitness.
	m_prob->objfun(m_container[idx].cur_f,m_container[idx].cur_x);
	// Best decision vector is current decision vector, best fitness is current fitness, best constraints are current constraints.
	m_container[idx].best_x = m_container[idx].cur_x;
	m_container[idx].best_f = m_container[idx].cur_f;
	m_container[idx].best_c = m_container[idx].cur_c;
	// Update the champion.
	update_champion(idx);
	// Update the domination lists.
	update_dom_list(idx);
}

/// Copy constructor.
/**
 * Will perform a deep copy of all the elements.
 *
 * @param[in] p population used to initialise this.
 */
population::population(const population &p):m_prob(p.m_prob->clone()),m_container(p.m_container),m_dom_list(p.m_dom_list),
	m_champion(p.m_champion),m_drng(p.m_drng),m_urng(p.m_urng) {}

/// Assignment operator.
/**
 * Performs a deep copy of all the elements of p into this.
 *
 * @param[in] p population to be assigned to this.
 *
 * @return reference to this.
 */
population &population::operator=(const population &p)
{
	if (this != &p) {
		pagmo_assert(m_prob && p.m_prob);
		// Perform the copies.
		if (*m_prob != p.problem()) {
			pagmo_throw(value_error,"cannot assign population with different problem");
		}
		m_container = p.m_container;
		m_dom_list = p.m_dom_list;
		m_champion = p.m_champion;
		m_drng = p.m_drng;
		m_urng = p.m_urng;
	}
	return *this;
}

/// Get constant reference to individual at position n.
/**
 * Will fail if idx is greater than size().
 *
 * @param idx positional index of the individual to get.
 *
 * @return const reference to the individual at position idx.
 *
 * @throws index_error if idx is not smaller than size().
 */
const population::individual_type &population::get_individual(const size_type &idx) const
{
	if (idx >= size()) {
		pagmo_throw(index_error,"invalid index");
	}
	return m_container[idx];
}

/// Get domination list.
/**
 * Will return a vector containing the indices of the individuals dominated by the individual in position idx. Will fail if
 * idx is not smaller than size().
 *
 * @param[in] idx position of the individual whose domination list will be retrieved.
 *
 * @return const reference to the vector containing the indices of the individuals dominated by the individual at position idx.
 *
 * @throws index_error if idx is not smaller than size().
 */
const std::vector<population::size_type> &population::get_domination_list(const size_type &idx) const
{
	if (idx >= size()) {
		pagmo_throw(index_error,"invalid index");
	}
	return m_dom_list[idx];
}

/// Get position of worst individual.
/**
 * The worst individual is the one dominating the smallest number of other individuals in the population.
 *
 * @return the positional index of the worst individual.
 */
population::size_type population::get_worst_idx() const
{
	if (!size()) {
		pagmo_throw(value_error,"empty population, cannot compute position of worst individual");
	}
	container_type::const_iterator it = std::max_element(m_container.begin(),m_container.end(),domination_comp(*this));
	return boost::numeric_cast<size_type>(std::distance(m_container.begin(),it));
}

/// Get position of best individual.
/**
* The best individual is the one dominating the highest number of other individuals in the population.
 *
 * @return the positional index of the best individual.
 */
population::size_type population::get_best_idx() const
{
	if (!size()) {
		pagmo_throw(value_error,"empty population, cannot compute position of best individual");
	}
	container_type::const_iterator it = std::min_element(m_container.begin(),m_container.end(),domination_comp(*this));
	return boost::numeric_cast<size_type>(std::distance(m_container.begin(),it));
}

/// Return terse human-readable representation.
/**
 * Will return a formatted string displaying:
 * - description of the problem (output of problem::base::human_readable())
 * - number of individuals.
 *
 * For use in pagmo::island.
 *
 * @return string with terse human-readable representation of the population.
 */
std::string population::human_readable_terse() const
{
	std::ostringstream oss;
	oss << (*m_prob) << '\n';
	oss << "Number of individuals: " << size() << '\n';
	return oss.str();
}

/// Return human-readable representation.
/**
 * Will return a formatted string displaying:
 * - the output of human_readable_terse(),
 * - list of individuals,
 * - champion.
 *
 * @return string with terse human-readable representation of the population.
 */
std::string population::human_readable() const
{
	std::ostringstream oss;
	oss << human_readable_terse();
	if (size()) {
		oss << "\nList of individuals:\n";
		for (size_type i = 0; i < size(); ++i) {
			oss << '#' << i << ":\n";
			oss << m_container[i] << "\tDominates:\t\t\t" << m_dom_list[i] << '\n';
		}
	}
	if (m_champion.x.size()) {
		oss << "Champion:\n";
		oss << m_champion << '\n';
	} else {
		pagmo_assert(!size());
		oss << "No champion yet.\n";
	}
	return oss.str();
}

// Update the champion with individual in position idx, if better or if the champion has not been set yet.
void population::update_champion(const size_type &idx)
{
	pagmo_assert(idx < m_container.size());
	if (!m_champion.x.size() || m_prob->compare_fc(m_container[idx].cur_f,m_container[idx].cur_c,m_champion.f,m_champion.c)) {
		m_champion.x = m_container[idx].cur_x;
		m_champion.f = m_container[idx].cur_f;
		m_champion.c = m_container[idx].cur_c;
	}
}

/// Set the decision vector of individual at position idx to x.
/**
 * Will update best values of individual and champion if needed. Will fail if problem::base::verify_x() on x returns false.
 *
 * @param[in] idx positional index of the individual to be set.
 * @param[in] x decision vector to be set for the individual at position idx.
 */
void population::set_x(const size_type &idx, const decision_vector &x)
{
	if (idx >= size()) {
		pagmo_throw(index_error,"invalid individual position");
	}
	if (!m_prob->verify_x(x)) {
		pagmo_throw(value_error,"decision vector is not compatible with problem");

	}
	// Set decision vector.
	m_container[idx].cur_x = x;
	// Update current fitness vector.
	m_prob->objfun(m_container[idx].cur_f,x);
	// Update current constraints vector.
	m_prob->compute_constraints(m_container[idx].cur_c,x);
	// If needed, update the best decision, fitness and constraint vectors for the individual.
	if (m_prob->compare_fc(m_container[idx].cur_f,m_container[idx].cur_c,m_container[idx].best_f,m_container[idx].best_c)) {
		m_container[idx].best_x = m_container[idx].cur_x;
		m_container[idx].best_f = m_container[idx].cur_f;
		m_container[idx].best_c = m_container[idx].cur_c;
	}
	// Update the champion.
	update_champion(idx);
	// Updated domination lists.
	update_dom_list(idx);
}

/// Append individual with given decision vector.
/**
 * A new individual with decision vector x will be appended at the end of the population.
 * Velocity will be initialised randomly (as described in reinit()).
 *
 * @param[in] x decision vector of the individual to be appended.
 */
void population::push_back(const decision_vector &x)
{
	if (!m_prob->verify_x(x)) {
		pagmo_throw(value_error,"decision vector is not compatible with problem");

	}
	// Store sizes temporarily.
	const fitness_vector::size_type f_size = m_prob->get_f_dimension();
	const constraint_vector::size_type c_size = m_prob->get_c_dimension();
	const decision_vector::size_type p_size = m_prob->get_dimension();
	// Push back an empty individual.
	m_container.push_back(individual_type());
	m_dom_list.push_back(std::vector<size_type>());
	// Resize individual's elements.
	m_container.back().cur_x.resize(p_size);
	m_container.back().cur_v.resize(p_size);
	m_container.back().cur_c.resize(c_size);
	m_container.back().cur_f.resize(f_size);
	m_container.back().best_x.resize(p_size);
	m_container.back().best_c.resize(c_size);
	m_container.back().best_f.resize(f_size);
	// Set the individual.
	set_x(m_container.size() - 1,x);
	// Initialise randomly the velocity vector.
	init_velocity(m_container.size() - 1);
	// Force update best values with current ones.
	m_container.back().best_x = m_container.back().cur_x;
	m_container.back().best_c = m_container.back().cur_c;
	m_container.back().best_f = m_container.back().cur_f;
}

/// Set the velocity vector of individual at position idx.
/**
 * Will fail if dimension of v differs from the problem dimension.
 *
 * @param[in] idx positional index of the individual to be set.
 * @param[in] v velocity vector to be set for the individual at position idx.
 */
void population::set_v(const size_type &idx, const decision_vector &v)
{
	if (idx >= size()) {
		pagmo_throw(index_error,"invalid individual position");
	}
	if (v.size() != this->problem().get_dimension()) {
		pagmo_throw(value_error,"decision vector is not compatible with problem");
	}
	// Set decision vector.
	m_container[idx].cur_v = v;
}

/// Get constant reference to internal problem::base object.
/**
 * @return const reference to internal problem::base object.
 */
const problem::base &population::problem() const
{
	return *m_prob;
}

/// Get constant reference to internal champion member.
/**
 * @return const reference to internal champion_type instance.
 */
const population::champion_type &population::champion() const
{
	if (!m_champion.x.size()) {
		pagmo_assert(!size());
		pagmo_throw(value_error,"champion has not been determined yet");
	}
	return m_champion;
}

/// Get the size of the population.
/**
 * @return number of individuals in the population.
 */
population::size_type population::size() const
{
	return m_container.size();
}

/// Iterator to the beginning of the population.
/**
 * @return iterator to the first individual.
 */
population::const_iterator population::begin() const
{
	return m_container.begin();
}

/// Iterator to the end of the population.
/**
 * @return iterator to the position one past the last individual.
 */
population::const_iterator population::end() const
{
	return m_container.end();
}

/// Number of dominated individuals.
/**
 * Get the number of individuals dominated by input individual ind.
 *
 * @param[in] ind input individual.
 *
 * @return number of individuals dominated by ind.
 */
population::size_type population::n_dominated(const individual_type &ind) const
{
	size_type retval = 0;
	for (size_type i = 0; i < m_container.size(); ++i) {
		if (m_prob->compare_fc(ind.cur_f,ind.cur_c,
			m_container[i].cur_f,m_container[i].cur_c))
		{
			++retval;
		}
	}
	return retval;
}

/// Overload stream operator for pagmo::population.
/**
 * Equivalent to printing pagmo::population::human_readable() to stream.
 *
 * @param[in] s stream to which the population will be sent.
 * @param[in] pop population to be sent to stream.
 *
 * @return reference to s.
 */
std::ostream &operator<<(std::ostream &s, const population &pop)
{
	s << pop.human_readable();
	return s;
}

/// Overload stream operator for pagmo::population::individual_type.
/**
 * Equivalent to printing pagmo::population::individual_type::human_readable() to stream.
 *
 * @param[in] s stream to which the individual will be sent.
 * @param[in] ind individual to be sent to stream.
 *
 * @return reference to s.
 */
std::ostream &operator<<(std::ostream &s, const population::individual_type &ind)
{
	s << ind.human_readable();
	return s;
}

/// Overload stream operator for pagmo::population::champion_type.
/**
 * Equivalent to printing pagmo::population::champion_type::human_readable() to stream.
 *
 * @param[in] s stream to which the champion will be sent.
 * @param[in] champ champion to be sent to stream.
 *
 * @return reference to s.
 */
std::ostream &operator<<(std::ostream &s, const population::champion_type &champ)
{
	s << champ.human_readable();
	return s;
}

}
