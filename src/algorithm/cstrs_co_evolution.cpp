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

#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <string>
#include <vector>

#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../problem/cstrs_co_evolution.h"
#include "../types.h"
#include "base.h"
#include "cstrs_co_evolution.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Constructs a co-evolution adaptive algorithm
 *
 * @param[in] original_algo pagmo::algorithm to use as 'original' optimization method
 * @param[in] original_algo_penalties pagmo::algorithm to use as 'original' optimization method 
 * for population representing the penaltiesweights
 * @param[in] pop_penalties_size population size for the penalty encoding population.
 * @param[in] gen number of generations.
 * @param[in] method the method used for the population 2.
 * Three possibililties are available: SIMPLE, SPLIT_NEQ_EQ and SPLIT_CONSTRAINTS.
 * The simple one is the original version of the Coello/He implementation. The SPLIT_NEQ_EQ,
 * splits the equalities and inequalities constraints in two different sets for the
 * penalty weigths, containing respectively inequalities and equalities weigths. The
 * SPLIT_CONSTRAINTS splits the constraints in M set of weigths with M the number of
 * constraints.
 * @param[in] pen_lower_bound the lower boundary used for penalty.
 * @param[in] pen_upper_bound the upper boundary used for penalty.
 * @param[in] ftol stopping criteria on the f tolerance.
 * @param[in] xtol stopping criteria on the x tolerance.
 * @throws value_error if stop is negative
 */
cstrs_co_evolution::cstrs_co_evolution(const base &original_algo, 
									   const base &original_algo_penalties, int pop_penalties_size,
									   int gen,method_type method, double pen_lower_bound,
									   double pen_upper_bound,
									   double ftol, double xtol):
	base(),m_original_algo(original_algo.clone()), m_original_algo_penalties(original_algo_penalties.clone()),
	m_gen(gen),m_pop_penalties_size(pop_penalties_size),m_method(method),
	m_pen_lower_bound(pen_lower_bound),m_pen_upper_bound(pen_upper_bound),m_ftol(ftol),m_xtol(xtol)
{
	if(gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if(pop_penalties_size <= 0) {
		pagmo_throw(value_error,"the population size of the penalty weights must be greater than 0");
	}
	if(pen_lower_bound>=pen_upper_bound){
		pagmo_throw(value_error,"Lower Bound of penalty coefficients must be smaller than Upper Bound");
	}
}

/// Copy constructor.
cstrs_co_evolution::cstrs_co_evolution(const cstrs_co_evolution &algo):
	base(algo),m_original_algo(algo.m_original_algo->clone()),
	m_original_algo_penalties(algo.m_original_algo_penalties->clone()),m_gen(algo.m_gen),
	m_pop_penalties_size(algo.m_pop_penalties_size),m_method(algo.m_method),
	m_pen_lower_bound(algo.m_pen_lower_bound),m_pen_upper_bound(algo.m_pen_upper_bound),
	m_ftol(algo.m_ftol),m_xtol(algo.m_xtol)
{}

/// Clone method.
base_ptr cstrs_co_evolution::clone() const
{
	return base_ptr(new cstrs_co_evolution(*this));
}

/// Evolve implementation.
/**
 * Run the co-evolution algorithm
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */
void cstrs_co_evolution::evolve(population &pop) const
{	
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const population::size_type pop_size = pop.size();
	const problem::base::size_type prob_dimension = prob.get_dimension();

	// get the constraints dimension
	problem::base::c_size_type prob_c_dimension = prob.get_c_dimension();

	//We perform some checks to determine wether the problem/population are suitable for co-evolution
	if(prob_c_dimension < 1) {
		pagmo_throw(value_error,"The problem is not constrained and co-evolution is not suitable to solve it");
	}
	if(prob.get_f_dimension() != 1) {
		pagmo_throw(value_error,"The problem is multiobjective and co-evolution is not suitable to solve it");
	}

	// Get out if there is nothing to do.
	if(pop_size == 0) {
		return;
	}

	//get the dimension of the chromosome of P2
	unsigned int pop_2_dim = 0;
	switch(m_method)
	{
	case algorithm::cstrs_co_evolution::SIMPLE:
	{
		pop_2_dim = 2;
		break;
	}
	case algorithm::cstrs_co_evolution::SPLIT_NEQ_EQ:
	{
		pop_2_dim = 4;
		break;
	}
	case algorithm::cstrs_co_evolution::SPLIT_CONSTRAINTS:
		pop_2_dim = 2*prob.get_c_dimension();
		break;
	default: 
		pagmo_throw(value_error,"The constraints co-evolutionary method is not valid.");
		break;
	}

	// split the population into two populations
	// the population P1 associated to the modified problem with penalized fitness
	// and population P2 encoding the penalty weights
	
	// Populations size
	population::size_type pop_1_size = pop_size;
	population::size_type pop_2_size = m_pop_penalties_size;

	//Creates problem associated to P2
	problem::cstrs_co_evolution_penalty prob_2(prob,pop_2_dim,pop_2_size);
	prob_2.set_bounds(m_pen_lower_bound,m_pen_upper_bound);

	//random initialization of the P2 chromosome (needed for the fist generation)
	std::vector<decision_vector> pop_2_x(pop_2_size);
	std::vector<fitness_vector> pop_2_f(pop_2_size);

	for(population::size_type j=0; j<pop_2_size; j++) {
		pop_2_x[j] = decision_vector(pop_2_dim,0.);
		// choose random coefficients between lower bound and upper bound
		for(population::size_type i=0; i<pop_2_dim;i++) {
			pop_2_x[j][i] = boost::uniform_real<double>(m_pen_lower_bound,m_pen_upper_bound)(m_drng);
		}
	}

	//vector of the population P1. Initialized with clones of the original population
	std::vector<population> pop_1_vector;
	for(population::size_type i=0; i<pop_2_size; i++){
		pop_1_vector.push_back(population(pop));
	}

	// Main Co-Evolution loop
	for(int k=0; k<m_gen; k++) {
		// for each individuals of pop 2, evolve the current population,
		// and store the position of the feasible idx
		for(population::size_type j=0; j<pop_2_size; j++) {

			problem::cstrs_co_evolution prob_1(prob, pop_1_vector.at(j), m_method);

			// modify the problem by setting decision vector encoding penalty
			// coefficients w1 and w2 in prob 1
			prob_1.set_penalty_coeff(pop_2_x.at(j));

			// creating the POPULATION 1 instance based on the
			// updated prob 1

			// prob_1 is a BASE_META???? THE CLONE OF prob_1 IN POP_1 IS AT THE LEVEL OF
			// THE BASE CLASS AND NOT AT THE LEVEL OF THE BASE_META, NO?!?
			population pop_1(prob_1,0);

			// initialize P1 chromosomes. The fitnesses related to problem 1 are computed
			for(population::size_type i=0; i<pop_1_size; i++) {
				pop_1.push_back(pop_1_vector.at(j).get_individual(i).cur_x);
			}

			// evolve the P1 instance
			m_original_algo->evolve(pop_1);

			//updating the original problem population (computation of fitness and constraints)
			pop_1_vector.at(j).clear();
			for(population::size_type i=0; i<pop_1_size; i++){
				pop_1_vector.at(j).push_back(pop_1.get_individual(i).cur_x);
			}

			// set up penalization variables needs for the population 2
			// the constraints has not been evaluated yet.
			prob_2.update_penalty_coeff(j,pop_2_x.at(j),pop_1_vector.at(j));

		}
		// creating the POPULATION 2 instance based on the
		// updated prob 2
		population pop_2(prob_2,0);

		// compute the fitness values of the second population
		for(population::size_type i=0; i<pop_2_size; i++) {
			pop_2.push_back(pop_2_x[i]);
		}

		m_original_algo_penalties->evolve(pop_2);

		// store the new chromosomes
		for(population::size_type i=0; i<pop_2_size; i++) {
			pop_2_x[i] = pop_2.get_individual(i).cur_x;
			pop_2_f[i] = pop_2.get_individual(i).cur_f;
		}

		// Check the exit conditions (every 40 generations, just as DE)
		if(k % 40 == 0) {
			// finds the best population
			population::size_type best_idx = 0;
			for(population::size_type j=1; j<pop_2_size; j++) {
				if(pop_2_f[j][0] < pop_2_f[best_idx][0]) {
					best_idx = j;
				}
			}

			const population &current_population = pop_1_vector.at(best_idx);

			decision_vector tmp(prob_dimension);

			double dx = 0;
			for(decision_vector::size_type i=0; i<prob_dimension; i++) {
				tmp[i] = current_population.get_individual(current_population.get_worst_idx()).best_x[i] -
						current_population.get_individual(current_population.get_best_idx()).best_x[i];
				dx += std::fabs(tmp[i]);
			}

			if(dx < m_xtol ) {
				if (m_screen_output) {
					std::cout << "Exit condition -- xtol < " << m_xtol << std::endl;
				}
				break;
			}

			double mah = std::fabs(current_population.get_individual(current_population.get_worst_idx()).best_f[0] -
					current_population.get_individual(current_population.get_best_idx()).best_f[0]);

			if(mah < m_ftol) {
				if(m_screen_output) {
					std::cout << "Exit condition -- ftol < " << m_ftol << std::endl;
				}
				break;
			}

			// outputs current values
			if(m_screen_output) {
				std::cout << "Generation " << k << " ***" << std::endl;
				std::cout << "    Best global fitness: " << current_population.champion().f << std::endl;
				std::cout << "    xtol: " << dx << ", ftol: " << mah << std::endl;
			}
		}
	}

	// find the best fitness population in the final pop
	population::size_type best_idx = 0;
	for(population::size_type j=1; j<pop_2_size; j++) {
		if(pop_2_f[j][0] < pop_2_f[best_idx][0]) { 
			best_idx = j;
		}
	}

	// store the final population in the main population
	// can't avoid to recompute the vectors here, otherwise we
	// clone the problem stored in the population with the
	// pop = operator!
	pop.clear();
	for(population::size_type i=0; i<pop_1_size; i++) {
		pop.push_back(pop_1_vector.at(best_idx).get_individual(i).cur_x);
	}
}

/// Algorithm name
std::string cstrs_co_evolution::get_name() const
{
	return m_original_algo->get_name() + "[Co-Evolution]";
}

/// Get a copy of the internal local algorithm.
/**
 * @return algorithm::base_ptr to a copy of the internal local algorithm.
 */
base_ptr cstrs_co_evolution::get_algorithm() const
{
	return m_original_algo->clone();
}

/// Set algorithm.
/**
 * A copy of the input algorithm will be set as the internal local algorithm.
 *
 * @param[in] algo algorithm to be set as local algorithm.
 */
void cstrs_co_evolution::set_algorithm(const base &algo)
{
	m_original_algo = algo.clone();
}

/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string cstrs_co_evolution::human_readable_extra() const
{
	std::ostringstream s;
	s << "algorithms: " << m_original_algo->get_name() << " - " << m_original_algo_penalties->get_name() << " ";
	s << "\n\tConstraints handled with co-evolution algorithm";
	return s.str();
}
}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::cstrs_co_evolution)
