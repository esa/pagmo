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

// 04/01/2009: Initial version by Francesco Biscani.

#include <algorithm>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include <limits>

#include "problem/base.h"
#include "problem/base_stochastic.h"
#include "exceptions.h"
#include "population.h"
#include "rng.h"
#include "types.h"
#include "util/racing.h"
#include "util/race_pop.h"

#include "algorithm/base.h"
#include "problem/con2uncon.h"

namespace pagmo
{

/// TODO: check if this is really needed!!!!
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
 * @param[in] seed rng seed (used to initialize the pop and in race)
 *
 * @throw value_error if n is negative.
 */
population::population(const problem::base &p, int n, const boost::uint32_t &seed):m_prob(p.clone()), m_pareto_rank(n), m_crowding_d(n), m_drng(seed),m_urng(seed)
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
		// Initialise randomly the individual.
		reinit(i);
	}
}

/// Copy constructor.
/**
 * Will perform a deep copy of all the elements.
 *
 * @param[in] p population used to initialise this.
 */
population::population(const population &p):m_prob(p.m_prob->clone()),m_container(p.m_container),m_dom_list(p.m_dom_list),m_dom_count(p.m_dom_count),
	m_champion(p.m_champion), m_pareto_rank(p.m_pareto_rank), m_crowding_d(p.m_crowding_d),m_drng(p.m_drng),m_urng(p.m_urng)
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
		m_pareto_rank = p.m_pareto_rank;
		m_crowding_d = p.m_crowding_d;
		m_drng = p.m_drng;
		m_urng = p.m_urng;
	}
	return *this;
}

// Update the domination list and the domination count when the individual at position n has changed
void population::update_dom(const size_type &n)
{
	// The algorithm works as follow:
	// 1) For each element in m_dom_list[n] decrease the domination count by one. (m_dom_count[m_dom_list[n][j]] -= 1)
	// 2) We empty the dom_list and reinitialize m_dom_count[n] = 0-
	// 3) We loop over the population (j) and construct again m_dom_list[n] and m_dom_count,
	//    taking care to also keep m_dom_list[j] correctly updated

	const size_type size = m_container.size();
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

/// Computes the mean curent velocity of all individuals in the population
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

/// Get Pareto rank
/**
 * Will return the Pareto rank for the requested individual idx (that is the Pareto front it belongs to, starting from 0,1,2....N).
 * A call to population::update_pareto_information() is needed if the population has
 * changed since the last time the Pareto rank was computed
 *
 * @param[in] idx position of the individual whose Pareto rank will be returned
 *
 * @return the Pareto rank of indiviual idx
 *
 * @throws index_error if idx is not smaller than m_pareto_rank.size().
 */
population::size_type population::get_pareto_rank(const size_type &idx) const
{
	if (idx >= m_pareto_rank.size()) {
		pagmo_throw(index_error,"invalid index");
	}
	return m_pareto_rank[idx];
}

/// Get Crowding Distance
/**
 * Will return the crowding distance for the requested individual idx. A call to population::update_pareto_information() is needed
 * if the population has changed since the last time the Crowding Distance was computed. The crowding distance is computed as
 * defined in Deb's work
 *
 * @see Deb, K. and Pratap, A. and Agarwal, S. and Meyarivan, T., "A fast and elitist multiobjective genetic algorithm: NSGA-II"
 *
 * @param[in] idx position of the individual whose Crowding Distance  will be returned
 *
 * @return the Crowding Distance of indiviual idx
 *
 * @throws index_error if idx is not smaller than m_crowding_distance.size().
 */
double population::get_crowding_d(const size_type &idx) const
{
	if (idx >= m_crowding_d.size()) {
		pagmo_throw(index_error,"invalid index");
	}
	return m_crowding_d[idx];
}

// This functor is used to sort (minimization assumed) along a particular fitness dimension.
// Needed for the computations of the crowding distance
struct one_dim_fit_comp {
	one_dim_fit_comp(const population &pop, fitness_vector::size_type dim):m_pop(pop), m_dim(dim) {};
	bool operator()(const population::size_type& idx1, const population::size_type& idx2) const
	{
		return m_pop.get_individual(idx1).cur_f[m_dim] < m_pop.get_individual(idx2).cur_f[m_dim];
	}
	const population& m_pop;
	fitness_vector::size_type m_dim;
};


/// Update Pareto Information
/**
 * Computes all pareto fronts, updates the pareto rank and the crowding distance of each individual.
 * Member variables for rank and crowding distance are set to zero and domination lists and
 * domination count are used to for the computation.
 */

void population::update_pareto_information() const {
	// Population size can change between calls and m_pareto_rank, m_crowding_d are updated if necessary
	m_pareto_rank.resize(size());
	m_crowding_d.resize(size());

	// We initialize the ranks and all distances to zero
	std::fill(m_pareto_rank.begin(), m_pareto_rank.end(), 0);
	std::fill(m_crowding_d.begin(), m_crowding_d.end(), 0);

	// We define some utility vectors .....
	std::vector<population::size_type> F,S;

	// And make a copy of the domination count (number of individuals that dominating one individual)
	std::vector<population::size_type> dom_count_copy(m_dom_count);

	// 1 - Find the first Pareto Front
	for (population::size_type idx = 0; idx < m_container.size(); ++idx){
		if (m_dom_count[idx] == 0) {
			F.push_back(idx);
		}
	}

	unsigned int irank = 1;

	// We loop to find subsequent fronts
	while (F.size()!=0) {
		// update crowding distance of the current pareto front
		population::update_crowding_d(F);
		//For each individual F in the current front
		for (population::size_type i=0; i < F.size(); ++i) {
			//For each individual dominated by F
			for (population::size_type j=0; j<m_dom_list[F[i]].size(); ++j) {
				dom_count_copy[m_dom_list[F[i]][j]]--;
				if (dom_count_copy[m_dom_list[F[i]][j]] == 0){
					S.push_back(m_dom_list[F[i]][j]);
					m_pareto_rank[m_dom_list[F[i]][j]] = irank;
				}
			}
		}
		F = S;
		S.clear();
		irank++;
	}
}


/// Update Crowding Distance
/**
 * This method computes the crowding distance of the entire population. This distance is saved
 * in the member variable m_crowding_d.
 *
 * @see Deb, K. and Pratap, A. and Agarwal, S. and Meyarivan, T., "A fast and elitist multiobjective genetic algorithm: NSGA-II"
 *
 * @param[in] I indices of the individuals of the pareto front
 *
 */
void population::update_crowding_d(std::vector<population::size_type> I) const {

	size_type lastidx = I.size() - 1;

	// we construct the comparison functor along the first fitness component
		one_dim_fit_comp funct(*this,0);

	// we loop along fitness components
		for (fitness_vector::size_type i = 0; i < problem().get_f_dimension(); ++i) {
			funct.m_dim = i;
				// we sort I along the fitness_dimension i
		std::sort(I.begin(),I.end(), funct );
			// assign Inf to the boundaries
			m_crowding_d[I[0]] = std::numeric_limits<double>::max();
			m_crowding_d[I[lastidx]] = std::numeric_limits<double>::max();
			//and compute the crowding distance
			double df = get_individual(I[lastidx]).cur_f[i] - get_individual(I[0]).cur_f[i];
			for (population::size_type j = 1; j < lastidx; ++j) {
				if (df == 0.0) { 						// handles the case in which the pareto front collapses to one single point
					m_crowding_d[I[j]] += 0.0;			// avoiding creation of nans that can't be serialized
				} else {
					m_crowding_d[I[j]] += (get_individual(I[j+1]).cur_f[i] - get_individual(I[j-1]).cur_f[i])/df;
				}
			}
		}
}

/// Computes and returns the population Pareto fronts
/**
 * This method computes all Pareto Fronts of the population, returning the positional indices
 * of the individuals belonging to each Pareto front.
 *
 * @return a vector containing, for each Pareto front, a vector of the individuals idx that belong to
 * each front
 */
std::vector<std::vector<population::size_type> > population::compute_pareto_fronts() const {
	std::vector<std::vector<population::size_type> > retval;

	// Be sure to have actual information about pareto rank
	population::update_pareto_information();

	for (population::size_type idx = 0; idx < size(); ++idx) {
		if (m_pareto_rank[idx] >= retval.size()) {
			retval.resize(m_pareto_rank[idx] + 1);
		}
		retval[m_pareto_rank[idx]].push_back(idx);
	}

	return retval;
}

/// Compute and return the ideal objective vector
/**
 * This method returns the ideal objective vector for the current optimal pareto set.
 * The components of the ideal objective vector are defined as:
 * \f[ z_i^{\mbox{ideal}} = \mbox{min}_{x \in X : x \mbox{ is Pareto Optimal}} f_i(x) \f]
 *
 * @return the ideal objective vector for the current optimal pareto set
 */
fitness_vector population::compute_ideal() const {
	update_pareto_information();

	fitness_vector ideal(problem().get_f_dimension(),std::numeric_limits<double>::max());
	for (population::size_type idx = 0; idx < size(); ++idx) {
		if (m_pareto_rank[idx] == 0) { //it is in the first pareto front
			for(fitness_vector::size_type i = 0; i < ideal.size(); ++i) {
				if (m_container[idx].cur_f[i] < ideal[i]) {
					ideal[i] = m_container[idx].cur_f[i];
				}
			}
		}
	}
	return ideal;
}

/// Compute and return the nadir objective vector
/**
 * This method returns the nadir objective vector for the current optimal pareto set.
 * The components of the nadir objective vector are defined as:
 * \f[ z_i^{\mbox{nadir}} = \mbox{max}_{x \in X : x \mbox{ is Pareto Optimal}} f_i(x) \f]
 *
 * @return the nadir objective vector for the current optimal pareto set
 */
fitness_vector population::compute_nadir() const {
	update_pareto_information();

	fitness_vector nadir(m_champion.f);
	for (population::size_type idx = 0; idx < size(); ++idx) {
		if (m_pareto_rank[idx] == 0) { //it is in the first pareto front
			for(fitness_vector::size_type i = 0; i < nadir.size(); ++i) {
				if (m_container[idx].cur_f[i] > nadir[i]) {
					nadir[i] = m_container[idx].cur_f[i];
				}
			}
		}
	}
	return nadir;
}

/// Crowded comparison functor.
/**
 * A binary functor that can be used to sort population individuals with respect
 * to the crowded comparison operator (assumes m_problem is multi-objective)
 *
 * An individual is better than a second one, if and only if:
 * (i) it has a lower Pareto rank; or
 * (ii) it has a larger crowding distance and the same Pareto rank
 *
 * @param[in] pop population over which the crowded comparison operator should operate
 *
 * @throws value_error if the number of individuals is less than 2, or if the underlying problem is not multi-objective
 */
population::crowded_comparison_operator::crowded_comparison_operator(const population &pop):m_pop(pop)
{
	if (m_pop.size() < 2) {
		pagmo_throw(value_error, "the population seems to contain one or zero individuals, comparisons between individuals do not make sense");
	}
	if (m_pop.problem().get_f_dimension() < 2) {
		pagmo_throw(value_error, "The crowded comparison operator can only operate on multi-objective problems.");
	}
}

/// Crowded comparison functor operator()(ind,ind)
/**
 *
 * NOTE: The user has to make sure that the Pareto information (which includes the crowding distance)
 * of m_pop is computed via update_pareto_information() before the use of this
 * operator.
 *
 * @param[in] i1 reference to the first individual to be compared
 * @param[in] i2 reference to the second individual to be compared
 *
 * @return true if the first individual is better than the second individual; false otherwise.
 *
 * @throw value_error if i1 or i2 are not within the population supplied
 */
bool population::crowded_comparison_operator::operator()(const individual_type &i1, const individual_type &i2) const
{
	if (!(&i1 >= &m_pop.m_container.front() && &i1 <= &m_pop.m_container.back())) {
		pagmo_throw(value_error, "operator called on individuals that do not belong to the population");
	}

	if (!(&i2 >= &m_pop.m_container.front() && &i2 <= &m_pop.m_container.back())) {
		pagmo_throw(value_error, "operator called on individuals that do not belong to the population");
	}
	const size_type idx1 = &i1 - &m_pop.m_container.front(), idx2 = &i2 - &m_pop.m_container.front();
	if (m_pop.m_pareto_rank[idx1] == m_pop.m_pareto_rank[idx2]) {
		return (m_pop.m_crowding_d[idx1] > m_pop.m_crowding_d[idx2]);
	}
	else {
		return (m_pop.m_pareto_rank[idx1] < m_pop.m_pareto_rank[idx2]);
	}
}

/// Crowded comparison functor operator()(idx,idx)
/**
 *
 * NOTE: The user has to make sure that the Pareto information (which includes the crowding distance)
 * of m_pop is computed via update_pareto_information() before the use of this
 * operator.
 *
 * @param[in] idx1 index of the first individual to be compared
 * @param[in] idx2 index of the second individual to be compared
 *
 * @return true if the first individual is better than the second individual; false otherwise.
 *
 * @throw index_error if idx1 or idx2 is out of bound
 */
bool population::crowded_comparison_operator::operator()(const size_type &idx1, const size_type &idx2) const
{
	if (idx1 >= m_pop.size() || idx2 >= m_pop.size()) {
		pagmo_throw(index_error, "operator called on out of bound indexes");
	}

	if (m_pop.m_pareto_rank[idx1] == m_pop.m_pareto_rank[idx2]) {
		return (m_pop.m_crowding_d[idx1] > m_pop.m_crowding_d[idx2]);
	}
	else {
		return (m_pop.m_pareto_rank[idx1] < m_pop.m_pareto_rank[idx2]);
	}
}

/// Trivial comparison functor
/**
 * A binary functor that can be used to sort population individuals in a
 * single objective case. (constraints are accounted for)
 *
 * NOTE: This functor uses the virtual method compare_fc assuming a weak strict
 * ordering implemented. If the user reimplements such a virtual method at the problem level,
 * he needs to make sure this condition is met (or pay the consequences :)
 *
 * @param[in] pop population over which the comparison operator should operate
 */
population::trivial_comparison_operator::trivial_comparison_operator(const population &pop):m_pop(pop) {}

/// Trivial comparison operator
/**
 * @param[in] i1 reference to the first individual to be compared
 * @param[in] i2 reference to the second individual to be compared
 *
 * @return true if the first individual is better than the second individual; false otherwise.
 *
 * @throw value_error if i1 or i2 is not within the population supplied
 */
bool population::trivial_comparison_operator::operator()(const individual_type &i1, const individual_type &i2) const
{
	if (!(&i1 >= &m_pop.m_container.front() && &i1 <= &m_pop.m_container.back())) {
		pagmo_throw(value_error, "operator called on individuals that do not belong to the population");
	}

	if (!(&i2 >= &m_pop.m_container.front() && &i2 <= &m_pop.m_container.back())) {
		pagmo_throw(value_error, "operator called on individuals that do not belong to population");
	}
	const size_type idx1 = &i1 - &m_pop.m_container.front(), idx2 = &i2 - &m_pop.m_container.front();
	return m_pop.problem().compare_fc(m_pop.get_individual(idx1).cur_f, m_pop.get_individual(idx1).cur_c, m_pop.get_individual(idx2).cur_f,m_pop.get_individual(idx2).cur_c);
}

/// Trivial comparison operator
/**
 * @param[in] idx1 index of the first individual to be compared
 * @param[in] idx2 index of the second individual to be compared
 *
 * @return true if the first individual is better than the second individual; false otherwise.
 *
 * @throw index_error if idx1 or idx2 is out of bound
 */
bool population::trivial_comparison_operator::operator()(const size_type &idx1, const size_type &idx2) const
{
	if (idx1 >= m_pop.size() || idx2 >= m_pop.size()) {
		pagmo_throw(index_error, "operator called on out of bound indexes");
	}
	return m_pop.problem().compare_fc(m_pop.get_individual(idx1).cur_f, m_pop.get_individual(idx1).cur_c, m_pop.get_individual(idx2).cur_f,m_pop.get_individual(idx2).cur_c);
}

/// Get position of worst individual.
/**
 * The definition of what makes an individual worst with respect to another differs in single objective
 * optimization from multiple objectives optimization.
 *
 * Single objective optimization: the comparison operator is the virtual method problem::compare_fc.
 *
 * Multi objective optimization: the crowded comparison operator is used. Note that with respect to
 * what originally defined by Deb in "A Fast and Elitist Multiobjective Genetic Algorithm: NSGA II",
 * we do not use the front rank,  but the m_dom_count (which is related but not identical).
 *
 *
 * NOTE: population.get_worst_idx assumes a weak strict ordering defined in problem::compare_fc. If the user
 * reimplements such a virtual method at the problem level, he needs to make sure this condition
 * is met (or pay the consequences :)
 *
 * @return the positional index of the worst individual.
 */
population::size_type population::get_worst_idx() const
{
	if (!size()) {
		pagmo_throw(value_error,"empty population, cannot compute position of worst individual");
	}
	container_type::const_iterator it;
	if (m_prob->get_f_dimension() == 1) {
		it = std::max_element(m_container.begin(),m_container.end(),trivial_comparison_operator(*this));
	}
	else {
		update_pareto_information();
		it = std::max_element(m_container.begin(),m_container.end(),crowded_comparison_operator(*this));
	}
	return boost::numeric_cast<size_type>(std::distance(m_container.begin(),it));
}

/// Get position of best individual.
/**
 * The definition of what makes an individual best with respect to another differs in single objective
 * optimization from multiple objectives optimization.
 *
 * Single objective optimization: the comparison operator is the virtual method problem::compare_fc.
 *
 * Multi objective optimization: the crowded comparison operator is used. Note that with respect to
 * what originally defined by Deb in "A Fast and Elitist Multiobjective Genetic Algorithm: NSGA II",
 * we do not use the front rank,  but the m_dom_count (which is related but not identical).
 *
 * NOTE: population.get_best_idx assumes a weak strict ordering defined in problem::compare_fc. If the user
 * reimplements such a virtual method at the problem level, he needs to make sure this condition
 * is met (or pay the consequences :)
 *
 * @return the positional index of the best individual.
 */
population::size_type population::get_best_idx() const
{
	if (!size()) {
		pagmo_throw(value_error,"empty population, cannot compute position of best individual");
	}
	container_type::const_iterator it;
	if (m_prob->get_f_dimension() == 1) {
		it = std::min_element(m_container.begin(),m_container.end(),trivial_comparison_operator(*this));
	}
	else {
		update_pareto_information();
		it = std::min_element(m_container.begin(),m_container.end(),crowded_comparison_operator(*this));
	}	return boost::numeric_cast<size_type>(std::distance(m_container.begin(),it));
}

/// Get positions of N best individuals.
/**
 * The definition of what makes an individual best with respect to another differs in single objective
 * optimization from multiple objectives optimization.
 *
 * Single objective optimization: the comparison operator is the virtual method problem::compare_fc.
 *
 * Multi objective optimization: the crowded comparison operator is used. Note that with respect to
 * what originally defined by Deb in "A Fast and Elitist Multiobjective Genetic Algorithm: NSGA II",
 * we do not use the front rank,  but the m_dom_count (which is related but not identical).
 *
 * NOTE: population.get_best_idx assumes a weak strict ordering defined in problem::compare_fc. If the user
 * reimplements such a virtual method at the problem level, he needs to make sure this condition
 * is met (or pay the consequences :)
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
	if (m_prob->get_f_dimension() == 1) {
		std::sort(retval.begin(),retval.end(),trivial_comparison_operator(*this));
	}
	else {
		update_pareto_information();
		std::sort(retval.begin(),retval.end(),crowded_comparison_operator(*this));
	}
	retval.resize(N);
	return retval;
}


/// Repairs the individual at the position idx.
/**
 * This methods repairs an infeasible individual to make it feasible. The method uses the repairing
 * algorithm provided by the user. It minimizes the constraints violation.
 *
 * @param[in] idx index of the individual to repair
 * @param[in] repair_algo algorithm to be used to repair the individual. Should be an algorithm
 * working with a population of size 1.
 *
 * @throws index_error if idx is larger than the population size
 */
void population::repair(const population::size_type &idx, const algorithm::base_ptr &repair_algo)
{
	if (idx >= size()) {
		pagmo_throw(index_error,"invalid individual position");
	}

	const decision_vector &current_x = m_container[idx].cur_x;
	const constraint_vector &current_c = m_container[idx].cur_c;

	// if feasible, nothing is done
	if(m_prob->feasibility_c(current_c)) {
		return;
	}

	problem::con2uncon feasibility_problem(*m_prob,problem::con2uncon::FEASIBILITY);

	population pop_repair(feasibility_problem);
	pop_repair.clear();
	pop_repair.push_back(current_x);

	repair_algo->evolve(pop_repair);

	this->set_x(idx,pop_repair.champion().x);

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
	m_container.erase(m_container.begin() + idx);
	m_dom_count.erase(m_dom_count.begin() + idx);
	m_dom_list.erase(m_dom_list.begin() + idx);
	// Since an element is erased indexes in dom_list need an update
	for (population::size_type i=0; i<m_dom_list.size(); ++i){
		for(population::size_type j=0; j<m_dom_list[i].size();++j) {
			if (m_dom_list[i][j] == idx) {
				m_dom_list[i].erase(m_dom_list[i].begin()+j);
				// If we did not erase the last individual
				if (m_dom_list[i].size() > j) {
					//check the next individual which would be skipped otherwise
					if (m_dom_list[i][j] > idx) m_dom_list[i][j]--;
				}
			}
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
 * Will clear the container of individuals, the domination lists, pareto ranks, crowding distance and the champion.
 * The problem and random number generators are left untouched.
 */
void population::clear()
{
	m_container.clear();
	m_dom_list.clear();
	m_dom_count.clear();
	m_crowding_d.clear();
	m_pareto_rank.clear();
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

/// Race the individuals in the population
/**
 * Perform racing on the individuals. Alternative to racing the whole population,
 * user can specify to race on only a subset of the individual. This is simply
 * a wrapper over the race_pop function in util::racing.
 *
 * @param[in] n_final Desired number of winners.
 * @param[in] min_trials Minimum number of trials to be executed before dropping individuals.
 * @param[in] max_count Maximum number of objective function evaluation before the race ends.
 * @param[in] delta Confidence level for statistical testing.
 * @param[in] active_set Indices of individuals that should participate in the race. If empty, race on the whole population.
 * @param[in] race_best If true winners are the best, otherwise winners are the worst
 * @param[in] screen_output If true some screen output is produced
 *
 * @return Indices of the individuals that remain in the race in the end, a.k.a the winners.
 *
 * @see pagmo::util::racing::race
 */
std::pair<std::vector<population::size_type>, unsigned int> population::race(const size_type n_final, const unsigned int min_trials, const unsigned int max_count, double delta, const std::vector<size_type>& active_set, const bool race_best, const bool screen_output) const
{
	unsigned int seed = m_urng();
	util::racing::race_pop m_race_pop(*this, seed);
	return m_race_pop.run(n_final, min_trials, max_count, delta, active_set, util::racing::race_pop::MAX_BUDGET, race_best, screen_output);
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
