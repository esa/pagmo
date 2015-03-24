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

// Test code for the noisy meta-problem

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cassert>
#include <numeric>

#include "../src/pagmo.h"
#include "../src/util/race_pop.h"

using namespace pagmo;
using namespace util::racing;

const double EPS = 10e-9;
rng_uint32 rng_seed_provider;

// This function takes a non-stochastic problem and creates a random
// population associated to it. Then it transforms the problem into
// a stochastic one (using the noisy meta-problem), and creates an
// associated population that contains the same individuals ordered
// with respect to the denoised problem objective function.
// It then tests that race returns, as it should [0,1,2,3,....,end_size-1].
// Such a check cannot be strict as for the stochastic nature of race.
// So some error is allowed. In particular
// it is only requested that a permutation of the correct indexes is returned
// and some small error in rankings is allowed:
//
// Examples: 
//			[0]     --> Pass (correct answer)
//			[1]     --> Pass (allowed error met)
//			[2]     --> Failed (allowed error not met)
//			[0,1]   --> Pass (correct answer)
//			[1,0]   --> Pass (correct permutation)
//			[0,2]   --> Pass (allowed error met)
//			[4,0]   --> Failed (allowed error not met)
//			[0,1,2] --> Pass (correct answer)
//			[0,2,1] --> Pass (correct permutation)
//			[0,3,2] --> Pass (allowed error met)
//			[1,3,5] --> Failed (allowed error not met)
//
// The error is defined as the difference between the true ranks sum and the
// returned ranks sum. For example if race returns [1,3,5], then the error is
// (1+3+5) - (0+1+2) = 6. The allowed error is twice the size of the returned list,
// (i.e. [2,3,4] is still valid, but [3,4,5] not)


int test_racing(const problem::base_ptr& prob, population::size_type pop_size, 
				population::size_type end_size, double noise_std_dev = 0.1)
{
	// sanity checks on test inputs
	if(pop_size < end_size){
		std::cout << "End size larger than start size! Test can never pass -- test case invalid." << std::endl;
		return 1;
	}

	unsigned int seed = 123;
	// we create a random population associated to the input problem and compute its ordering
	population pop_original_prob(*prob, pop_size,seed);
	std::vector<population::size_type> best_idx_order = pop_original_prob.get_best_idx(pop_size);

	// we create the noisy version of the problem
	problem::noisy prob_noisy(*prob, 1, 0, noise_std_dev, problem::noisy::NORMAL, seed);
	std::cout << prob_noisy << std::endl;

	// and a random empty population associated to it
	population pop(prob_noisy,0,seed);

	// We push_back the old population chromosomes in the new population according to their rank 
	// If we do this, by definition, race should return [0,1,2,3,4,5,...,end_size-1]
	for(population::size_type i = 0; i < pop_size; i++){
		pop.push_back(pop_original_prob.get_individual(best_idx_order[i]).cur_x);
	}

	std::pair<std::vector<population::size_type>, unsigned int> race_results = pop.race(end_size, 50, 5000, 0.05);
	std::vector<population::size_type> winners = race_results.first;
	double ground_truth = ((winners.size()-1)+1) * (winners.size()-1) / 2;
	double obtained = std::accumulate(winners.begin(),winners.end(),0.0,std::plus<population::size_type>());
    double acceptable = 2 * winners.size();

	std::cout << prob->get_name() << std::endl;
	std::cout << "\tRace winners: " << winners << std::endl;
	std::cout << "\tError: " << obtained-ground_truth << std::endl;
	std::cout << "\tAcceptable: " << acceptable << std::endl;

	if(winners.size() != end_size){
		std::cout << " Winner list size failed." << std::endl;
		return 1;
	}

	if (obtained-ground_truth > acceptable ) {
		std::cout << "\tFAILED" << std::endl;
		return 1;
	}

	std::cout << "\tPASSED" << std::endl<< std::endl;
	return 0;
}

int test_racing_worst(const problem::base_ptr& prob, population::size_type pop_size, 
		population::size_type end_size, double noise_std_dev = 0.1)
{
	std::cout << ":::test_racing_worst" << std::endl;
	// sanity checks on test inputs
	if(pop_size < end_size){
		std::cout << "End size larger than start size! Test can never pass -- test case invalid." << std::endl;
		return 1;
	}

	unsigned int seed = 1234;
	// we create a random population associated to the input problem and compute its ordering
	population pop_original_prob(*prob, pop_size,seed);
	std::vector<population::size_type> best_idx_order = pop_original_prob.get_best_idx(pop_size);

	// we create the noisy version of the problem
	problem::noisy prob_noisy(*prob, 1, 0, noise_std_dev, problem::noisy::NORMAL, seed);
	std::cout << prob_noisy << std::endl;

	// and a random empty population associated to it
	population pop(prob_noisy,0,seed);

	// We push_back the old population chromosomes in the new population according to their rank
	// If we do this, by definition, race should return [0,1,2,3,4,5,...,end_size-1]
	for(population::size_type i = 0; i < pop_size; i++){
		pop.push_back(pop_original_prob.get_individual(best_idx_order[pop_size-1-i]).cur_x);
	}

	// NOTE: now winners are actually losers
	std::pair<std::vector<population::size_type>, unsigned int> race_results = pop.race(end_size, 50, 5000, 0.05, std::vector<population::size_type>(), false);
	std::vector<population::size_type> winners = race_results.first;
	double ground_truth = ((winners.size()-1)+1) * (winners.size()-1) / 2;
	double obtained = std::accumulate(winners.begin(),winners.end(),0.0,std::plus<population::size_type>());
	double acceptable = winners.size();

	std::cout << prob->get_name() << std::endl;
	std::cout << "\tRace winners: " << winners << std::endl;
	std::cout << "\tError: " << obtained-ground_truth << std::endl;
	std::cout << "\tAcceptable: " << acceptable << std::endl;

	if(winners.size() != end_size){
		std::cout << " Winner list size failed." << std::endl;
		return 1;
	}

	if (obtained-ground_truth > acceptable ) {
		std::cout << "\tFAILED" << std::endl;
		return 1;
	}

	std::cout << "\tPASSED" << std::endl<< std::endl;
	return 0;
}

// This second test only checks the implementation of the active_set
// so that subset_size individuals are raced and all of them returned
int test_racing_subset(const problem::base_ptr& prob)
{
	std::cout << "Testing race() on a subset of population" << std::endl;

	unsigned int pop_size = 50;	
	unsigned int subset_size = 20;
	std::vector<pagmo::population::size_type> active_set(subset_size);
	for(unsigned int i = 0; i < subset_size; i++){
		active_set[i] = i*2;
	}

	problem::noisy prob_noisy(*prob, 1, 0, 0.3, problem::noisy::NORMAL);
	population pop(prob_noisy, pop_size);

	// Race until subset_size individuals remain, with subset_size individuals being active
	// at the beginning. This implies that the winners should be exactly the same as
	// active_set.
	std::pair<std::vector<pagmo::population::size_type>, unsigned int> race_results = pop.race(subset_size, 0, 500, 0.05, active_set);
	std::vector<pagmo::population::size_type> winners = race_results.first;
	if(winners.size() != subset_size){
		std::cout << " Winner list size failed." << std::endl;	
		return 1;
	}

	std::sort(active_set.begin(), active_set.end());		
	std::sort(winners.begin(), winners.end());

	if(winners != active_set){
		std::cout << " Racing subset result failed." << std::endl;
		return 1;
	}

	std::cout << prob->get_name() << " racing subset passed." << std::endl;

	return 0;
}

/// Check if the caching mechanism is working
int test_racing_cache(const problem::base_ptr& prob)
{
	std::cout << "Testing the caching mechanism of racing" << std::endl;
	
	unsigned int seed = 123;
	// we create the noisy version of the problem
	problem::noisy prob_noisy(*prob, 1, 0, 0.5, problem::noisy::NORMAL, seed);
	std::cout << prob_noisy << std::endl;

	population pop(prob_noisy, 5, seed);

	util::racing::race_pop race_pop_dev(pop, seed);

	// [0,1,2,3,4] then again
	std::vector<population::size_type> active_set;
	for(unsigned int i = 0; i <= 4; i++){
		active_set.push_back(i);
	}
	population::size_type n_final = 1;
	std::pair<std::vector<population::size_type>, unsigned int> res1 = race_pop_dev.run(n_final, 0, 500, 0.05, active_set, race_pop::MAX_BUDGET, true, false);
	std::pair<std::vector<population::size_type>, unsigned int> res2 = race_pop_dev.run(n_final, 0, 500, 0.05, active_set, race_pop::MAX_BUDGET, true, false);

	std::cout << "\tfevals: First race consumed " << res1.second << " and second consumed " << res2.second << std::endl;

	if(res1.first != res2.first){
		std::cout << "\tFAILED Caching: Winners are different!" << std::endl;
	}

	if(res2.second > 0){
		std::cout << "\tFAILED Caching: Second race consumed fevals!" << std::endl;
		return 1;
	}

	std::cout << "\tPASSED Caching." << std::endl;
	return 0;
}

int test_racing_cache_transfer(const problem::base_ptr &prob)
{
	std::cout << "Testing the caching mechanism (memory transfer aspect) of racing" << std::endl;
	
	unsigned int seed = 123;
	// we create the noisy version of the problem
	problem::noisy prob_noisy(*prob, 1, 0, 0.5, problem::noisy::NORMAL, seed);
	std::cout << prob_noisy << std::endl;

	population pop(prob_noisy, 5, seed);

	util::racing::race_pop race_pop_dev1(pop, seed);
	util::racing::race_pop race_pop_dev2(pop, seed);

	// Race [0,1,2,3,4] using race_pop_dev1, then transfer the memory to
	// race_pop_dev2, which should accept all the cache data as they contain
	// identical population, and same seed is being used.
	std::vector<population::size_type> active_set;
	for(unsigned int i = 0; i <= 4; i++){
		active_set.push_back(i);
	}
	population::size_type n_final = 1;

	// Race on device 1
	std::pair<std::vector<population::size_type>, unsigned int> res1 = race_pop_dev1.run(n_final, 0, 500, 0.05, active_set, race_pop::MAX_BUDGET, true, false);	

	// Transfer to device 2
	race_pop_dev2.inherit_memory(race_pop_dev1);

	// Race on device 2
	std::pair<std::vector<population::size_type>, unsigned int> res2 = race_pop_dev2.run(n_final, 0, 500, 0.05, active_set, race_pop::MAX_BUDGET, true, false);

	std::cout << "\tfevals: First race consumed " << res1.second << " and second consumed " << res2.second << std::endl;

	if(res1.first != res2.first){
		std::cout << "\tFAILED Caching: Winners are different!" << std::endl;
	}

	if(res2.second > 0){
		std::cout << "\tFAILED Caching: Second race consumed fevals!" << std::endl;
		return 1;
	}

	std::cout << "\tPASSED Caching (memory transfer)." << std::endl;
	return 0;
}

/// Check if the returned mean fitness by race_pop is sensible
int test_racing_get_mean_fitness(const problem::base_ptr &prob)
{
	std::cout << "Testing get_mean_fitness()" << std::endl;
	
	unsigned int seed = 123;
	// we create the noisy version of the problem
	problem::noisy prob_noisy(*prob, 1, 0, 0.5, problem::noisy::NORMAL, seed);
	std::cout << prob_noisy << std::endl;

	population pop(prob_noisy, 5, seed);

	util::racing::race_pop race_pop_dev(pop, seed);

	// Race everyone
	std::vector<population::size_type> active_set;
	for(unsigned int i = 0; i < pop.size(); i++){
		active_set.push_back(i);
	}
	population::size_type n_final = 1;

	std::pair<std::vector<population::size_type>, unsigned int> res = race_pop_dev.run(n_final, 0, 5000, 0.05, active_set, race_pop::MAX_BUDGET, true, true);

	// Check that the dimension is correct
	unsigned int winner_idx = res.first[0];
	active_set.clear();
	active_set.push_back(winner_idx);
	fitness_vector mean_fitness_race = race_pop_dev.get_mean_fitness(active_set)[0];
	if(mean_fitness_race.size() != prob->get_f_dimension()){
		std::cout << "\tFAILED get_mean_fitness: Wrong fitness dimension" << std::endl;
		return 1;
	}
	
	// Check that the returned mean fitness is about the same as the fitness
	// obtained by repeated evaluation.
	unsigned int n_evals = 100;
	fitness_vector mean_fitness_repeated_eval(prob->get_f_dimension(),0);
	for(unsigned int i = 0; i < n_evals; i++){
		fitness_vector cur_f = prob->objfun(pop.get_individual(winner_idx).cur_x);
		for(unsigned int j = 0; j < prob->get_f_dimension(); j++){
			mean_fitness_repeated_eval[j] += cur_f[j] / (double)n_evals;
		}
	}

	// Tolerance is a bit large here as race could not and should not perform
	// too much evaluation, hence affecting the accuracy.
	double eps = 0.2;
	for(unsigned int i = 0; i < mean_fitness_race.size(); i++){
		if(fabs(mean_fitness_repeated_eval[i]-mean_fitness_race[i]) > eps){	
			std::cout << "\tFAILED get_mean_fitness: Dimension #"  << i << " had too much deviation" << std::endl;
			std::cout << "From repeated evaluation: " << mean_fitness_repeated_eval << std::endl;
			std::cout << "From racing: " << mean_fitness_race << std::endl;
			return 1;
		}
	}

	std::cout << "\tPASSED get_mean_fitness" << std::endl;
	return 0;
}

/// Check if the second type of constructor (population deferred) is sane
int test_race_pop_constructor(const problem::base_ptr& prob)
{
	std::cout << "Testing the constructor without population" << std::endl;
	
	unsigned int seed = 123;
	// we create the noisy version of the problem
	problem::noisy prob_noisy(*prob, 1, 0, 0.5, problem::noisy::NORMAL, seed);
	std::cout << prob_noisy << std::endl;

	population pop(prob_noisy, 5, seed);

	// Two construction methods
	util::racing::race_pop race_pop_dev1(pop, seed);
	util::racing::race_pop race_pop_dev2(seed);

	// [0,1,2,3,4] then again
	std::vector<population::size_type> active_set;
	for(unsigned int i = 0; i <= 4; i++){
		active_set.push_back(i);
	}
	population::size_type n_final = 1;
	std::pair<std::vector<population::size_type>, unsigned int> res1 = race_pop_dev1.run(n_final, 0, 500, 0.05, active_set, race_pop::MAX_BUDGET, true, false);

	// Register the population
	race_pop_dev2.register_population(pop);
	// Try to inherit memory from the other race structure
	race_pop_dev2.inherit_memory(race_pop_dev1);

	std::pair<std::vector<population::size_type>, unsigned int> res2 = race_pop_dev2.run(n_final, 0, 500, 0.05, active_set, race_pop::MAX_BUDGET, true, false);

	std::cout << "\tfevals: First race consumed " << res1.second << " and second consumed " << res2.second << std::endl;

	if(res1.first != res2.first){
		std::cout << "\tFAILED Caching: Winners are different!" << std::endl;
	}

	if(res2.second > 0){
		std::cout << "\tFAILED Caching: Second race consumed fevals!" << std::endl;
		return 1;
	}

	std::cout << "\tPASSED race_pop constructor without population." << std::endl;
	return 0;
}


int main()
{
	int dimension = 10;
	problem::base_ptr prob_ackley(new problem::ackley(dimension));
	problem::base_ptr prob_cec2006(new problem::cec2006(5));
	problem::base_ptr prob_zdt1(new problem::zdt(1, dimension));
	
	return test_racing(prob_ackley, 10, 2) ||
		   test_racing(prob_ackley, 20, 2) ||
		   test_racing(prob_ackley, 100, 5) ||
		   test_racing(prob_ackley, 5, 1) ||
		   test_racing_subset(prob_ackley) ||
		   test_racing(prob_cec2006, 10, 2) ||
		   test_racing(prob_cec2006, 20, 2) ||
		   test_racing(prob_cec2006, 100, 5, 0.05) ||
		   test_racing(prob_cec2006, 5, 1) ||
		   test_racing_subset(prob_cec2006) ||
		   test_racing(prob_zdt1, 10, 5, 0.03) ||
		   test_racing(prob_zdt1, 20, 5, 0.03) ||
		   test_racing(prob_zdt1, 30, 5, 0.03) ||
		   test_racing_subset(prob_zdt1) ||

		   test_racing_worst(prob_ackley, 10, 2) || 
		   test_racing_worst(prob_ackley, 20, 2) ||
		   test_racing_worst(prob_ackley, 100, 5) ||
		   test_racing_worst(prob_ackley, 5, 1) ||
		   test_racing_worst(prob_cec2006, 10, 2) ||
		   test_racing_worst(prob_cec2006, 20, 2) ||
		   test_racing_worst(prob_cec2006, 100, 5) ||
		   test_racing_worst(prob_cec2006, 5, 1) ||
		   test_racing_worst(prob_zdt1, 10, 5, 0.03) ||
		   test_racing_worst(prob_zdt1, 20, 5, 0.03) ||
		   test_racing_worst(prob_zdt1, 30, 5, 0.03) ||

		   test_racing_cache(prob_ackley) ||
		   test_racing_cache_transfer(prob_ackley) ||

		   test_racing_get_mean_fitness(prob_ackley) ||

		   test_race_pop_constructor(prob_ackley);
}
