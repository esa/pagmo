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
#include "../problem/antibodies_problem.h"
#include "../problem/con2uncon.h"
#include "../types.h"
#include "base.h"
#include "cstrs_immune_system.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Constructs an immune system constraints handling algorithm
 *
 * @param[in] original_algo pagmo::algorithm to use as 'original' optimization method. Its number of
 * generations should be set to 1.
 * @param[in] original_algo_immune pagmo::algorithm to use as 'original' optimization method
 * for the immune system
 * @param[in] gen number of generations.
 * @param[in] select_method the method used for selecting the antibodies.
 * @param[in] inject_method the method used for reinjecting the antibodies.
 * @param[in] distance_method the method used to compute the antibodies distance 
 * @param[in] phi the feasible fraction selection to compute the mean value
 * @param[in] gamma number of antigens selected / number of total antigens
 * @param[in] sigma number of antibodies / number of antigens
 * @param[in] ftol stopping criteria on the f tolerance
 * @param[in] xtol stopping criteria on the x tolerance
 * @throws value_error if gen is negative
 */
cstrs_immune_system::cstrs_immune_system(const base &original_algo, const base &original_algo_immune, int gen,
										 select_method_type select_method,
										 inject_method_type inject_method,
										 distance_method_type distance_method,
										 double phi,
										 double gamma,
										 double sigma,
										 double ftol, double xtol):
	base(),m_gen(gen),m_select_method(select_method),m_inject_method(inject_method),m_distance_method(distance_method),
	m_phi(phi),m_gamma(gamma),m_sigma(sigma),m_ftol(ftol),m_xtol(xtol)
{
	m_original_algo = original_algo.clone();
	m_original_algo_immune = original_algo_immune.clone();

	if(gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
}

/// Copy constructor.
cstrs_immune_system::cstrs_immune_system(const cstrs_immune_system &algo):
	base(algo),m_original_algo(algo.m_original_algo->clone()),
	m_original_algo_immune(algo.m_original_algo_immune->clone()),m_gen(algo.m_gen),
	m_select_method(algo.m_select_method),m_inject_method(algo.m_inject_method),m_distance_method(algo.m_distance_method),
	m_phi(algo.m_phi),m_gamma(algo.m_gamma),m_sigma(algo.m_sigma),m_ftol(algo.m_ftol), m_xtol(algo.m_xtol)
{}

/// Clone method.
base_ptr cstrs_immune_system::clone() const
{
	return base_ptr(new cstrs_immune_system(*this));
}

/// Evolve implementation.
/**
 * Run the co-evolution algorithm
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */
void cstrs_immune_system::evolve(population &pop) const
{	
	// store useful variables
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

	// generates the unconstrained problem
	problem::con2uncon prob_unconstrained(prob);

	// associates the population to this problem
	population pop_mixed(prob_unconstrained);
	std::vector<decision_vector> pop_mixed_c(pop_size);

	// initializaton of antigens vector and antibodies population

	// antigens population
	std::vector<decision_vector> pop_antigens;

	// vector containing the best antigens position
	std::vector<population::size_type> pop_antigens_pool;
	std::vector<population::size_type> pop_antibodies_pool;

	pop_mixed.clear();
	// the initial popluation contains the initial random population
	for(population::size_type i=0; i<pop_size; i++) {
		pop_mixed.push_back(pop.get_individual(i).cur_x);
	}

	// Main Co-Evolution loop
	for(int k=0; k<m_gen; k++) {

		pop_antigens.clear();

		// clearing the pools
		pop_antigens_pool.clear();
		pop_antibodies_pool.clear();

		// first of all we compute the constraints
		for(population::size_type i=0; i<pop_size; i++) {
			const population::individual_type &current_individual = pop_mixed.get_individual(i);
			pop_mixed_c[i] = prob.compute_constraints(current_individual.cur_x);
		}

		// we find if there are feasible individuals
		bool has_feasible = false;

		for(population::size_type i=0; i<pop_size; i++) {
			if(prob.feasibility_c(pop_mixed_c[i])) {
				has_feasible = true;
				break;
			}
		}

		// if feasible solutions exist in the population, the
		// antigens are selected by being close to the average fitness
		// of a sub population

		if(has_feasible) {
			// the pop_antigens_pool is based on the feasible population
			for(population::size_type i=0; i<pop_size; i++) {
				if(prob.feasibility_c(pop_mixed_c[i])) {
					pop_antigens_pool.push_back(i);
				} else {
					pop_antibodies_pool.push_back(i);
				}
			}

			population::size_type pop_antigens_pool_size = pop_antigens_pool.size();

			// we sort the pop_antigens_pool according to the fitness
			for(population::size_type i=0; i<pop_antigens_pool_size-1; i++) {
				const population::size_type &current_idx = pop_antigens_pool.at(i);
				const population::individual_type &current_individual = pop_mixed.get_individual(current_idx);

				for(population::size_type j=i+1; j<pop_antigens_pool_size; j++) {
					const population::size_type &current_second_idx = pop_antigens_pool.at(j);
					const population::individual_type &current_second_individual = pop_mixed.get_individual(current_second_idx);

					if(prob.compare_fitness(current_second_individual.cur_f, current_individual.cur_f)) {
						std::swap(pop_antigens_pool[i],pop_antigens_pool[j]);
					}
				}
			}

			// a subset of the best antigens in the pool is selected to compute the fitness average
			population::size_type pop_antigens_subset_size = std::max((population::size_type)(m_phi * pop_antigens_pool_size),(population::size_type)1);

			// we compute the mean fitness value from the subset of the best individuals in the population
			double mean_fitness = 0.;
			for(population::size_type i=0; i<pop_antigens_subset_size; i++) {
				const population::size_type &current_idx = pop_antigens_pool.at(i);
				mean_fitness += pop_mixed.get_individual(current_idx).cur_f.at(0);
			}
			mean_fitness /= pop_antigens_subset_size;

			// the population around this mean fitness is selected to fill the antigens population
			// finds the position of the individual closest to the mean fitness

			// from the antigens pool we select a fraction from ones that are closest to the mean value
			// according to Hjale and Lee. The idea here is to get as much information about the
			// feasible domain as possible.

			int mean_position_idx=0;
			for(population::size_type i=0; i<pop_antigens_pool_size; i++) {
				const population::size_type &current_idx = pop_antigens_pool.at(i);
				if(mean_fitness <= pop_mixed.get_individual(current_idx).cur_f.at(0)) {
					mean_position_idx = i;
					break;
				}
			}

			// generates the antigens population
			int pop_antigens_size = std::max((int)(m_gamma * pop_antigens_pool_size), 1);

			population::size_type begin_antigen_idx = mean_position_idx - pop_antigens_size/2;
			population::size_type end_antigen_idx = mean_position_idx + pop_antigens_size/2;

			// move the selection range depending on the boundaries
			if(mean_position_idx - pop_antigens_size/2 < 0) {
				begin_antigen_idx = 0;
				end_antigen_idx = std::min((int)(end_antigen_idx + (pop_antigens_size/2 - mean_position_idx)),(int)pop_antigens_pool_size);
			}
			if(mean_position_idx + pop_antigens_size/2 >= pop_antigens_size) {
				begin_antigen_idx = std::max((int)(begin_antigen_idx - (mean_position_idx + pop_antigens_size/2 - pop_antigens_size)),0);
				end_antigen_idx = pop_antigens_size;
			}

			// we select individuals around the mean
			for(population::size_type i=begin_antigen_idx; i<end_antigen_idx; i++) {
				population::size_type current_individual_idx = pop_antigens_pool.at(i);
				pop_antigens.push_back(pop_mixed.get_individual(current_individual_idx).cur_x);
			}

			if(pop_antigens.size() == 0) {
				pop_antigens.push_back(pop_mixed.get_individual(pop_antigens_pool.at(mean_position_idx)).cur_x);
			}

		} else {
			// no feasible founds, the antigen population
			// is selected accordingly to selected method

			switch(m_select_method) {
			case(BEST_ANTIBODY): {
				// Coello method the antigen population contains only one antigen which is the
				// individual with lowest constraints violation
				population::size_type best_idx = 0;
				for(population::size_type i=1; i<pop_size; i++) {
					if(prob.compare_constraints(pop_mixed_c[i], pop_mixed_c[best_idx])) {
						best_idx = i;
					}
				}
				pop_antigens_pool.push_back(best_idx);
				pop_antigens.push_back(pop_mixed.get_individual(best_idx).cur_x);

				// antibodies
				for(population::size_type i=0; i<pop_size; i++) {
					if(i != best_idx){
						pop_antibodies_pool.push_back(i);
					}
				}
				break;
			}
			case(INFEASIBILITY): {
				// Personal method where the infeasibility is used
				// to select the antigen population
				for(population::size_type i=0; i<pop_size; i++) {
					pop_antigens_pool.push_back(i);
				}

				population::size_type pop_antigens_pool_size = pop_antigens_pool.size();

				// we sort the pop_antigens_pool according to the fitness
				for(population::size_type i=0; i<pop_antigens_pool_size-1; i++) {
					const population::size_type &current_idx = pop_antigens_pool.at(i);

					for(population::size_type j=i+1; j<pop_antigens_pool_size; j++) {
						const population::size_type &current_second_idx = pop_antigens_pool.at(j);

						if(prob.compare_constraints(pop_mixed_c[current_second_idx], pop_mixed_c[current_idx])) {
							std::swap(pop_antigens_pool[i],pop_antigens_pool[j]);
						}
					}
				}

				// fill the antigens with the best ones
				population::size_type pop_antigens_size = std::max((population::size_type)(m_gamma * pop_antigens_pool_size), (population::size_type)1);

				for(population::size_type i=0; i<pop_antigens_size; i++) {
					population::size_type current_individual_idx = pop_antigens_pool.at(i);
					pop_antigens.push_back(pop_mixed.get_individual(current_individual_idx).cur_x);
				}

				// antibodies
				// not the best implementation, trying to find a better one...
				for(population::size_type i=0; i<pop_size; i++) {
					if(std::find(pop_antigens_pool.begin(), pop_antigens_pool.begin()+pop_antigens_size, i) == pop_antigens_pool.begin()+pop_antigens_size) {
						pop_antibodies_pool.push_back(i);
					}
				}
				break;
			}
			default: {
				pagmo_throw(value_error,"The antibody selection method must be either BEST_ANTIBODY or INFEASIBILITY.");
				break;
			}
			}
		}

		population::size_type initial_pop_antibodies_pool_size = pop_antibodies_pool.size();

		if(initial_pop_antibodies_pool_size != 0) {

			// initial pool size
			// random shuffle of the antibodies as the antibodies population
			// are randomly selected without replacement from the antibodies pool
			if(initial_pop_antibodies_pool_size > 1) {
				for (population::size_type i=0; i<initial_pop_antibodies_pool_size; i++)
				{
					int j = boost::uniform_int<int>(0, initial_pop_antibodies_pool_size - 1)(m_urng);
					std::swap(pop_antibodies_pool[i], pop_antibodies_pool[j]);
				}
			}

			// select the antibodies population with the requested size:
			// antibodies population size = 1/3 antigens population size

			population::size_type pop_antigens_size = pop_antigens.size();

			population::size_type min_individual_for_algo = 8;

			population::size_type pop_antibodies_size = std::max( (int)(m_sigma * pop_antigens_size), (int)min_individual_for_algo);
			pop_antibodies_size = std::min(pop_antibodies_size, initial_pop_antibodies_pool_size);

			//population::size_type pop_antibodies_size = std::max((int)(0.5 * pop_antigens_size), 6);
			//population::size_type pop_antibodies_size = std::max((int)(pop_antibodies_pool.size()), 6);

			// the problem can be updated with antigenes, need to be done here to avoid a cast
			problem::antibodies_problem prob_antibodies(prob, m_distance_method);
			prob_antibodies.set_antigens(pop_antigens);

			// immune system initialization
			population pop_antibodies(prob_antibodies);
			pop_antibodies.clear();

			for(population::size_type i=0; i<pop_antibodies_size; i++) {
				pop_antibodies.push_back(pop_mixed.get_individual(pop_antibodies_pool.at(i)).cur_x);
			}

			// ensure that the antibodies population has at least 6 individuals for de, sga, 8 for jde...
			if(min_individual_for_algo>pop_antibodies_size){
				population::size_type extra_antibodies_size = min_individual_for_algo - pop_antibodies_size;

				// add extra needed by randomly selecting the individuals in the pool
				for(population::size_type i=0; i<extra_antibodies_size; i++) {
					if(initial_pop_antibodies_pool_size > 1) {
						int j = boost::uniform_int<int>(0, initial_pop_antibodies_pool_size - 1)(m_urng);
						pop_antibodies.push_back(pop_mixed.get_individual(pop_antibodies_pool.at(j)).cur_x);
					} else {
						pop_antibodies.push_back(pop_mixed.get_individual(pop_antibodies_pool.at(0)).cur_x);
					}
				}
			}

			pop_antigens_size = pop_antigens.size();
			pop_antibodies_size = pop_antibodies.size();
			
			// run the immune system
			m_original_algo_immune->evolve(pop_antibodies);

			// sets the mixed population with all current best designs
			// and the copies of constraint conditioned designs
			switch(m_inject_method) {
			case(CHAMPION): {
				for(population::size_type i=0; i<initial_pop_antibodies_pool_size; i++) {
					const population::size_type &current_idx = pop_antibodies_pool.at(i);
					// multiple copies of the best antibody
					pop_mixed.set_x(current_idx, pop_antibodies.champion().x);
				}
				break;
			}
			case(BEST25): {
				std::vector<population::size_type> pop_best_25(pop_antibodies_size);
				for(population::size_type i=0; i<pop_antibodies_size; i++) {
					pop_best_25[i] = i;
				}

				for(population::size_type i=0; i<pop_antibodies_size-1; i++) {
					const population::size_type &current_idx = pop_best_25.at(i);
					const population::individual_type &current_individual = pop_mixed.get_individual(current_idx);

					for(population::size_type j=i+1; j<pop_antibodies_size; j++) {
						const population::size_type &current_second_idx = pop_best_25.at(j);
						const population::individual_type &current_second_individual = pop_mixed.get_individual(current_second_idx);

						if(prob.compare_fitness(current_second_individual.cur_f, current_individual.cur_f)) {
							std::swap(pop_best_25[i],pop_best_25[j]);
						}
					}
				}

				int best25_size = pop_antibodies_size/4;

				// multiple copies of the 25% bests
				for(population::size_type i=0; i<initial_pop_antibodies_pool_size; i++) {
					const population::size_type &current_idx = pop_antibodies_pool.at(i);
					pop_mixed.set_x(current_idx, pop_antibodies.get_individual(pop_best_25.at(i % best25_size)).cur_x);
				}
				break;
			}
			default: {
				pagmo_throw(value_error,"The antibody injection method must be either CHAMPION or BEST25.");
				break;
			}
			}
		}

		// for individuals, evolve the mixed population
		// which is an unconstrained problem
		// only one iteration should be done...
		m_original_algo->evolve(pop_mixed);

		// Check the exit conditions (every 40 generations, just as DE)
		if(k % 40 == 0) {
			decision_vector tmp(prob_dimension);

			double dx = 0;
			for(decision_vector::size_type i=0; i<prob_dimension; i++) {
				tmp[i] = pop_mixed.get_individual(pop_mixed.get_worst_idx()).best_x[i] - pop_mixed.get_individual(pop_mixed.get_best_idx()).best_x[i];
				dx += std::fabs(tmp[i]);
			}

			if(dx < m_xtol ) {
				if (m_screen_output) {
					std::cout << "Exit condition -- xtol < " << m_xtol << std::endl;
				}
				break;
			}

			double mah = std::fabs(pop_mixed.get_individual(pop_mixed.get_worst_idx()).best_f[0] - pop_mixed.get_individual(pop_mixed.get_best_idx()).best_f[0]);

			if(mah < m_ftol) {
				if(m_screen_output) {
					std::cout << "Exit condition -- ftol < " << m_ftol << std::endl;
				}
				break;
			}

			// outputs current values
			if(m_screen_output) {
				std::cout << "Generation " << k << " ***" << std::endl;
				std::cout << "    Best global fitness: " << pop.champion().f << std::endl;
				std::cout << "    xtol: " << dx << ", ftol: " << mah << std::endl;
			}
		}
	}

	// store the final population in the main population
	pop.clear();
	for(population::size_type i=0; i<pop_size; i++) {
		pop.push_back(pop_mixed.get_individual(i).cur_x);
	}
}

/// Algorithm name
std::string cstrs_immune_system::get_name() const
{
	return m_original_algo->get_name() + "[Immune]";
}

/// Get a copy of the internal local algorithm.
/**
 * @return algorithm::base_ptr to a copy of the internal local algorithm.
 */
base_ptr cstrs_immune_system::get_algorithm() const
{
	return m_original_algo->clone();
}

/// Set algorithm.
/**
 * A copy of the input algorithm will be set as the internal local algorithm
 * for immune system.
 *
 * @param[in] algo algorithm to be set as local algorithm.
 */
void cstrs_immune_system::set_algorithm_immune(const base &algo)
{
	m_original_algo_immune = algo.clone();
}
/// Get a copy of the internal local algorithm.
/**
 * @return algorithm::base_ptr to a copy of the internal local algorithm
 * for immune system.
 */
base_ptr cstrs_immune_system::get_algorithm_immune() const
{
	return m_original_algo_immune->clone();
}

/// Set algorithm.
/**
 * A copy of the input algorithm will be set as the internal local algorithm.
 *
 * @param[in] algo algorithm to be set as local algorithm.
 */
void cstrs_immune_system::set_algorithm(const base &algo)
{
	m_original_algo = algo.clone();
}

/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string cstrs_immune_system::human_readable_extra() const
{
	std::ostringstream s;
	s << "algorithms: " << m_original_algo->get_name() << " - " << m_original_algo_immune->get_name() << " ";
	s << "\n\tConstraints handled with immune system algorithm";
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::cstrs_immune_system)
