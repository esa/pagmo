/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
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

// Test code for the noisy meta-problem

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cassert>
#include "../src/pagmo.h"

using namespace pagmo;

const double EPS = 10e-9;
rng_uint32 rng_seed_provider;

bool is_eq(const fitness_vector & f1, const fitness_vector & f2, double eps){
	if(f1.size() != f2.size()) return false;
	for(unsigned int i = 0; i < f1.size(); i++){
		if(fabs(f1[i]-f2[i])>eps) return false;
	}
	return true;
}

// Receives a regular problem, and record the index of the best individual.
// Transform it to a stochastic one (currently using noisy), and then deliberately
// duplicate several copies of the best individual into the population. If race()
// is sane, all these duplicated copies should be in the final race winner list.
// Also check for basic sanity e.g. correct final race list size, etc.
int test_racing(const problem::base_ptr& prob, population::size_type start_size, population::size_type end_size, population::size_type num_best_copy)
{
	if(start_size < end_size){
		std::cout << "End size larger than start size! Test can never pass -- test case invalid." << std::endl;
		return 1;
	}
	if((1 + num_best_copy) > end_size){
		std::cout << "(Number of best copy + 1) larger than end size! Test can never pass -- test case invalid." << std::endl;
		return 1;
	}

	population pop_original_prob(*prob, start_size);

	population::size_type best_idx = pop_original_prob.get_best_idx();
	std::vector<population::size_type> best_idx_order = pop_original_prob.get_best_idx(start_size);

	std::cout << "original best order: " << best_idx_order << std::endl;

	problem::noisy prob_noisy(*prob, 1, 0, 0.1, problem::noisy::NORMAL);

	population pop(prob_noisy, start_size);

	for(population::size_type i = 0; i < start_size; i++){
		pop.set_x(i, pop_original_prob.get_individual(i).cur_x);
	}

	/*		
	for(population::size_type i = num_best_copy; i < start_size; i++){
		population::size_type some_bad_idx = best_idx_order[(start_size - 1)- (i % end_size)];
		pop.set_x(best_idx_order[i], pop_original_prob.get_individual(some_bad_idx).cur_x);
	}
	*/

	// If we do the following, by definition race should return [0, ... , num_best_copy-1]
	// PLUS best_idx, in no particular order.
	pop.set_x(best_idx, pop_original_prob.get_individual(best_idx).cur_x);
	for(population::size_type i = 0; i < num_best_copy; i++){
		pop.set_x(i, pop_original_prob.get_individual(best_idx).cur_x);
	}

	std::vector<population::size_type> winners = pop.race(end_size, 0, 500, 0.01);
	std::cout << "winners: " << winners << ", true best_idx = " << best_idx << std::endl;


	if(winners.size() != end_size){
		std::cout << " Winner list size failed." << std::endl;
		return 1;
	}

	std::vector<bool> found(start_size, false);

	for(population::size_type i = 0; i < winners.size(); i++){
		found[winners[i]] = true;
	}

	// Check if indeed every best copies are in place in the final winner list.
	for(population::size_type i = 0; i < pop.size(); i++){
		if( (i < num_best_copy || i == best_idx) && !found[i] ){
			std::cout << " Expected individual indexed " << i << " to be in the list but turned out not, failed! " << std::endl;
			std::cout << " Original best idx: " << best_idx << std::endl;
			std::cout << " Final winner list: " << winners << std::endl;
			std::cout << " Expected " << num_best_copy + 1 << " of the copies to be in the list." << std::endl;
			return 1;
		}
	}
	
	std::cout << prob->get_name() << " racing passed." << std::endl;
	return 0;
}

int test_racing_subset(const problem::base_ptr& prob)
{
	std::cout << "Testing race() on a subset of population" << std::endl;

	unsigned int N = 50;	
	unsigned int ActiveN = 20;
	std::vector<pagmo::population::size_type> active_set(20);
	for(unsigned int i = 0; i < ActiveN; i++){
		active_set[i] = i*2;
	}

	problem::noisy prob_noisy(*prob, 1, 0, 0.3, problem::noisy::NORMAL);
	population pop(prob_noisy, N);

	// Race until ActiveN individuals remain, with ActiveN individuals being active
	// at the beginning. This implies that the winners should be exactly the same as
	// active_set.
	std::vector<pagmo::population::size_type> winners = pop.race(ActiveN, 0, 500, 0.05, active_set);

	if(winners.size() != ActiveN){
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

// TODO: This test is fragile -- ``most" of the time it passes.
int main()
{
	int dimension = 10;
	problem::base_ptr prob_ackley(new problem::ackley(dimension));
	return test_racing(prob_ackley, 10, 2, 1) || 
		   test_racing(prob_ackley, 20, 2, 1) ||
		   test_racing(prob_ackley, 100, 50, 0) ||
		   test_racing(prob_ackley, 5, 1, 0) ||
		   test_racing_subset(prob_ackley);
}
