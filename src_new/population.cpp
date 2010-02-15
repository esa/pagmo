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

#include "problem/base.h"
#include "exceptions.h"
#include "population.h"
#include "rng.h"
#include "types.h"

namespace pagmo
{
// Function object to compare individuals according to their current fitness and constraints.
struct cur_fc_comp {
	cur_fc_comp(const problem::base &p):m_p(p) {}
	bool operator()(const population::individual_type &i1, const population::individual_type &i2) const
	{
		return m_p.compare_fc(i1.cur_f,i1.cur_c,i2.cur_f,i2.cur_c);
	}
	const problem::base &m_p;
};

/// Constructor from problem::base and number of individuals.
/**
 * Will store a copy of the problem and will initialise the population to n randomly-generated individuals.
 * Will fail if n is negative.
 *
 * @param[in] p problem::base that will be associated to the population.
 * @param[in] n integer number of individuals in the population.
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
	const decision_vector::size_type p_size = m_prob->get_dimension(), i_size = m_prob->get_i_dimension();
	pagmo_assert(p_size >= i_size);
	for (size_type i = 0; i < size; ++i) {
		// Push back an empty individual.
		m_container.push_back(individual_type());
		// Resize individual's elements.
		m_container.back().cur_x.resize(p_size);
		m_container.back().cur_v.resize(p_size);
		m_container.back().cur_c.resize(c_size);
		m_container.back().cur_f.resize(f_size);
		m_container.back().best_x.resize(p_size);
		m_container.back().best_c.resize(c_size);
		m_container.back().best_f.resize(f_size);
		// Initialise randomly the continuous part of the decision vector.
		for (decision_vector::size_type i = 0; i < p_size - i_size; ++i) {
			m_container.back().cur_x[i] = boost::uniform_real<double>(m_prob->get_lb()[i],m_prob->get_ub()[i])(m_drng);
		}
		// Initialise randomly the integer part of the decision vector.
		for (decision_vector::size_type i = p_size - i_size; i < p_size; ++i) {
			m_container.back().cur_x[i] = boost::uniform_int<int>(m_prob->get_lb()[i],m_prob->get_ub()[i])(m_urng);
		}
		// Initialise randomly the velocity vector.
		for (decision_vector::size_type i = 0; i < p_size; ++i) {
			// Initialise velocities so that in one tick the particles travel at most half the bounds distance.
			m_container.back().cur_v[i] = boost::uniform_real<double>(m_prob->get_lb()[i] / 2,m_prob->get_ub()[i] / 2)(m_drng);
			// Change randomly the sign of the velocity.
			m_container.back().cur_v[i] *= (m_drng() < .5) ? 1 : -1;
		}
		// Fill in the constraints part.
		m_prob->compute_constraints(m_container.back().cur_c,m_container.back().cur_x);
		// Compute the current fitness.
		m_prob->objfun(m_container.back().cur_f,m_container.back().cur_x);
		// Best decision vector is current decision vector, best fitness is current fitness, best constraints are current constraints.
		m_container.back().best_x = m_container.back().cur_x;
		m_container.back().best_f = m_container.back().cur_f;
		m_container.back().best_c = m_container.back().cur_c;
	}
	// Calculate the champion.
	container_type::iterator it = std::min_element(m_container.begin(),m_container.end(),cur_fc_comp(*m_prob));
	if (it != m_container.end()) {
		m_champion.x = it->cur_x;
		m_champion.f = it->cur_f;
		m_champion.c = it->cur_c;
	}
}

/// Copy constructor.
/**
 * Will perform a deep copy of all the elements.
 *
 * @param[in] p population used to initialise this.
 */
population::population(const population &p):m_prob(p.m_prob->clone()),m_container(p.m_container),m_champion(p.m_champion),m_drng(p.m_drng),m_urng(p.m_urng) {}

/// Default constructor.
/**
 * For use only by pagmo::island.
 */
population::population() {}

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
		// Perform the copies.
		m_prob = p.m_prob->clone();
		m_container = p.m_container;
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
 */
const population::individual_type &population::get_individual(const size_type &idx) const
{
	if (idx >= size()) {
		pagmo_throw(index_error,"invalid individual position");
	}
	return m_container[idx];
}

/// Get position of worst individual.
/**
 * problem::base::compare_fc() is used to rank the individuals according to their current fitnesses and constraints vectors.
 *
 * @return the positional index of the worst individual.
 */
population::size_type population::get_worst_idx() const
{
	if (!size()) {
		pagmo_throw(value_error,"empty population, cannot compute position of worst individual");
	}
	container_type::const_iterator it = std::max_element(m_container.begin(),m_container.end(),cur_fc_comp(*m_prob));
	return boost::numeric_cast<size_type>(std::distance(m_container.begin(),it));
}

/// Get position of best individual.
/**
 * problem::base::compare_fc() is used to rank the individuals according to their current fitnesses and constraints vectors.
 *
 * @return the positional index of the best individual.
 */
population::size_type population::get_best_idx() const
{
	if (!size()) {
		pagmo_throw(value_error,"empty population, cannot compute position of best individual");
	}
	container_type::const_iterator it = std::min_element(m_container.begin(),m_container.end(),cur_fc_comp(*m_prob));
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
		oss << "List of individuals:\n";
		for (size_type i = 0; i < size(); ++i) {
			oss << '#' << i << ":\n";
			oss << "\tDecision vector:\t\t" << m_container[i].cur_x << '\n';
			oss << "\tVelocity vector:\t\t" << m_container[i].cur_v << '\n';
			oss << "\tConstraint vector:\t\t" << m_container[i].cur_c << '\n';
			oss << "\tFitness vector:\t\t\t" << m_container[i].cur_f << '\n';
			oss << "\tBest decision vector:\t\t" << m_container[i].best_x << '\n';
			oss << "\tBest constraint vector:\t\t" << m_container[i].best_c << '\n';
			oss << "\tBest fitness vector:\t\t" << m_container[i].best_f << '\n';
		}
	}
	if (m_champion.x.size()) {
		oss << "Champion:\n";
		oss << "\tDecision vector:\t" << m_champion.x << '\n';
		oss << "\tConstraints vector:\t" << m_champion.c << '\n';
		oss << "\tFitness vector:\t\t" << m_champion.f << '\n';
	} else {
		pagmo_assert(!size());
		oss << "No champion yet.\n";
	}
	return oss.str();
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
	// If needed update the champion. Make sure with the assert that the champion exists. It
	// should be guaranteed at this point.
	pagmo_assert(m_champion.x.size());
	if (m_prob->compare_fc(m_container[idx].cur_f,m_container[idx].cur_c,m_champion.f,m_champion.c)) {
		m_champion.x = x;
		m_champion.f = m_container[idx].cur_f;
		m_champion.c = m_container[idx].cur_c;
	}
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

/// Overload stream operator for pagmo::population.
/**
 * Equivalent to printing population::human_readable() to stream.
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

}
