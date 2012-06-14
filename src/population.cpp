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
 * @throw value_error if n is negative.
 */
population::population(const problem::base &p, int n):m_prob(p.clone()), m_drng(rng_generator::get<rng_double>()),m_urng(rng_generator::get<rng_uint32>())
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
		m_dom_count.push_back(0);
		// Resize individual's elements.
		m_container.back().cur_x.resize(p_size);
		m_container.back().cur_v.resize(p_size);
		m_container.back().cur_c.resize(c_size);
		m_container.back().cur_f.resize(f_size);
		m_container.back().best_x.resize(p_size);
		m_container.back().best_c.resize(c_size);
		m_container.back().best_f.resize(f_size);
//std::cout << "constructor: " << m_dom_list.size() << " " << m_dom_count.size() << " " <<  size << std::endl;
		// Initialise randomly the individual.
		reinit(i);
	}
}

// Update the domination list and the domination count when the individual at position n has changed
void population::update_dom(const size_type &n)
{
	//The algorithm works as follow:
	// 1) For each element in m_dom_list[n] decrease the domination count by one. (m_dom_count[m_dom_list[n][j]] -= 1)
	// 2) We empty the dom_list and reinitialize m_dom_count[n] = 0-
	// 3) We loop over the population (j) and construct again m_dom_list[n] and m_dom_count, 
	//    taking care to also keep m_dom_list[j] correctly updated
	
	const size_type size = m_container.size();
//std::cout << "update_dom: " << m_dom_list.size() << " " << m_dom_count.size() << " " <<  size << std::endl;
	pagmo_assert(m_dom_list.size() == size && m_dom_count.size() == size && n < size);
	
	// Decrease the domination count for the individuals that were dominated
	for  (size_type i = 0; i < m_dom_list[n].size(); ++i) {
		m_dom_count[ m_dom_list[n][i] ]--;
	}
	
	// Empty the domination list of individual at position n.
	m_dom_list[n].clear();
	// Reset the dom_count of individual at position n.
	m_dom_count[n] = 0;
	
	for (size_type i = 0; i < size; ++i) {
		if (i != n) {
			// Check if individual in position i dominates individual in position n.
			if (m_prob->compare_fc(m_container[i].best_f,m_container[i].best_c,m_container[n].best_f,m_container[n].best_c)) {
				// Update the domination count in n. 
				m_dom_count[n]++;
				// Update the domination list in i. 
				//If n is already present, do nothing, otherwise push_back.
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
			if (m_prob->compare_fc(m_container[n].best_f,m_container[n].best_c,m_container[i].best_f,m_container[i].best_c)) {
				m_dom_list[n].push_back(i);
				m_dom_count[i]++;
			}
		}
	}
}

/// Computes and returns the population pareto fronts
/**
 * 
 * TBW
 */
std::vector<std::vector<population::size_type> > population::compute_pareto_fronts() const {
	std::vector<std::vector<population::size_type> > retval;
	std::vector<population::size_type> F,S;
	std::vector<population::size_type> dom_count_copy(m_dom_count);
	// We find the first Pareto Front
	for (population::size_type idx = 0; idx < m_container.size(); ++idx){
		if (m_dom_count[idx] == 0) {
			F.push_back(idx);
		}
	}
	// And if not empty, we push it back to retval
	if (F.size() != 0) retval.push_back(F);
	
	// We loop to find subsequent fronts
	while (F.size()!=0) {
		//For each individual F in the current front
		for (population::size_type i=0; i < F.size(); ++i) {
			//For each individual dominated by F
			for (population::size_type j=0; j<m_dom_list[F[i]].size(); ++j) {
				dom_count_copy[m_dom_list[F[i]][j]]--;
				if (dom_count_copy[m_dom_list[F[i]][j]] == 0) S.push_back(m_dom_list[F[i]][j]);
			}
		}
		F = S;
		S.clear();
		if (F.size() != 0) retval.push_back(F);
	}
	return retval;
	
}

// Init randomly the velocity of the individual in position idx.
void population::init_velocity(const size_type &idx)
{
	const decision_vector::size_type p_size = m_prob->get_dimension();
	double width = 0;
	for (decision_vector::size_type j = 0; j < p_size; ++j) {
		// Initialise velocities so that in one tick the particles travel
		// at most half the bounds distance.
		width = (m_prob->get_ub()[j] - m_prob->get_lb()[j]) / 2;
		m_container[idx].cur_v[j] = boost::uniform_real<double>(-width,width)(m_drng);
	}
}

double population::mean_velocity() const {
	const population::size_type pop_size(m_container.size());
	if (pop_size == 0) {
		pagmo_throw(zero_division_error,"Population has no individuals, no mean velocity can be computed.");
	}
	double ret=0, tmp;
	const decision_vector::size_type p_size = m_prob->get_dimension();

	for (population::size_type i = 0; i<pop_size; ++i) {
		tmp = 0;
		for (decision_vector::size_type j = 0; j < p_size; ++j) {
			tmp += m_container[i].cur_v[j]*m_container[i].cur_v[j];
		}
		ret += std::sqrt(tmp);
	}
	return (ret / pop_size);
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
 * will be evaluated, best_x and best_f are erasred and reset to the new values, the champion is updated.
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
	update_dom(idx);
}

/// Copy constructor.
/**
 * Will perform a deep copy of all the elements.
 *
 * @param[in] p population used to initialise this.
 */
population::population(const population &p):m_prob(p.m_prob->clone()),m_container(p.m_container),m_dom_list(p.m_dom_list),m_dom_count(p.m_dom_count),
	m_champion(p.m_champion),m_drng(p.m_drng),m_urng(p.m_urng)
{}

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
		m_prob = p.m_prob->clone();
		m_container = p.m_container;
		m_dom_list = p.m_dom_list;
		m_dom_count = p.m_dom_count;
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

/// Get domination count.
/**
 * Will return the domination count for the requested individual idx. That is the number of population individuals that dominate idx.
 *
 * @param[in] idx position of the individual whose domination count will be retrieved.
 *
 * @return the domination count
 *
 * @throws index_error if idx is not smaller than size().
 */
population::size_type population::get_domination_count(const size_type &idx) const
{
	if (idx >= size()) {
		pagmo_throw(index_error,"invalid index");
	}
	return m_dom_count[idx];
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

struct pair_order
{
	pair_order(const population &pop):m_pop(pop) {}
	bool operator()(const population::size_type &i1, const population::size_type &i2) const
	{
		return m_pop.get_domination_list(i1).size() > m_pop.get_domination_list(i2).size();
	}
	const population &m_pop;
};

/// Get positions of N best individuals.
/**
* The best individuals are the one dominating the highest number of other individuals
 * in the population.
 *
 * @return a std::vector of positional indexes of the best N individuals.
 * @throws value_error if N is larger than the population size or the population is empty
 */
std::vector<population::size_type> population::get_best_idx(const population::size_type& N) const
{
	if (!size()) {
		pagmo_throw(value_error,"empty population, cannot compute position of best individual");
	}
	if (N > size()) {
		pagmo_throw(value_error,"Best N individuals requested, but population has size smaller than N");
	}
	std::vector<population::size_type> retval;
	retval.reserve(size());
	for (population::size_type i=0; i<size(); ++i){
		retval.push_back(i);
	}
	pair_order po(*this);
	std::sort(retval.begin(),retval.end(),po);
	retval.resize(N);
	return retval;
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
			oss << "\tIs dominated by:\t\t" << m_dom_count[i] << "\tindividuals" << '\n';
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
	if (!m_champion.x.size() || m_prob->compare_fc(m_container[idx].best_f,m_container[idx].best_c,m_champion.f,m_champion.c)) {
		m_champion.x = m_container[idx].best_x;
		m_champion.f = m_container[idx].best_f;
		m_champion.c = m_container[idx].best_c;
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
	// NOTE: we update the bests in two cases:
	// - the bests are empty, meaning they are not defined and we are being called by push_back()
	// - the bests are defined, but they are worse than the currents.
	pagmo_assert((!m_container[idx].best_x.size() && !m_container[idx].best_f.size()) ||
		(m_container[idx].best_x.size() && m_container[idx].best_f.size()));
	if (!m_container[idx].best_x.size() ||
		m_prob->compare_fc(m_container[idx].cur_f,m_container[idx].cur_c,m_container[idx].best_f,m_container[idx].best_c))
	{
		m_container[idx].best_x = m_container[idx].cur_x;
		m_container[idx].best_f = m_container[idx].cur_f;
		m_container[idx].best_c = m_container[idx].cur_c;
	}
	// Update the champion.
	update_champion(idx);
	// Updated domination lists.
	update_dom(idx);
}

/// Erase individual idx
/**
 * The individual occupying position idx in the population will be erased from the population.
 * This method takes care to update accordingly the domination structures of the population
 * 
 * @param[in] idx index of the individual to be erased
 * 
 * @throws index_error if idx is out of range
 */

void population::erase(const population::size_type & idx) {

	pagmo_assert(m_dom_list.size() == size() && m_dom_count.size() == size() && idx < size());
	
	if (idx >= size()) {
		pagmo_throw(index_error,"invalid individual position");
	}
	for (population::size_type i = 0; i < m_dom_list[idx].size(); ++i) {
		m_dom_count[m_dom_list[idx][i]]--;
	}
//std::cout << "sizes: " << m_container.size() << " " << m_dom_count.size() << " " << m_dom_list.size() << std::endl;
	m_container.erase(m_container.begin() + idx);
	m_dom_count.erase(m_dom_count.begin() + idx);
	m_dom_list.erase(m_dom_list.begin() + idx);
	// Since an element is erased indexes in dom_list need an update
	for (population::size_type i=0; i<size(); ++i){
		for(population::size_type j=0; j<m_dom_list[i].size();++j) {
			if (m_dom_list[i][j] == idx) m_dom_list[i].erase(m_dom_list[i].begin()+j);
			else if (m_dom_list[i][j] > idx) m_dom_list[i][j]--;
  
		}
	}
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
	m_dom_count.push_back(0);
	// Resize individual's elements.
	m_container.back().cur_x.resize(p_size);
	m_container.back().cur_v.resize(p_size);
	m_container.back().cur_c.resize(c_size);
	m_container.back().cur_f.resize(f_size);
	// NOTE: do not allocate space for bests, as they are not defined yet. set_x will take
	// care of it.
	// Set the individual.
	set_x(m_container.size() - 1,x);
	// Initialise randomly the velocity vector.
	init_velocity(m_container.size() - 1);
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
		pagmo_throw(value_error,"velocity vector is not compatible with problem");
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

/// Clear population.
/**
 * Will clear the container of individuals, the domination lists and the champion. The problem and random number generators
 * are left untouched.
 */
void population::clear()
{
	m_container.clear();
	m_dom_list.clear();
	m_dom_count.clear();
	m_champion = champion_type();
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
 * Get the number of individuals in pop dominated by an input individual ind
 * If ind belongs to pop it is more efficient to use get_domination_list().size()
 *
 * @param[in] ind input individual.
 *
 * @return number of individuals dominated by ind.
 */
population::size_type population::n_dominated(const individual_type &ind) const
{
	size_type retval = 0;
	for (size_type i = 0; i < m_container.size(); ++i) {
		if (m_prob->compare_fc(ind.best_f,ind.best_c,
			m_container[i].best_f,m_container[i].best_c))
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
