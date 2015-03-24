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

// Test code for util::race_pop

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cassert>
#include "../src/pagmo.h"
#include "../src/util/race_algo.h"

using namespace pagmo;

/* Test strategy:
 * Generally, higher gen_num gives better results for any specific
 * meta-heuristic. This test instantiates different instances of an algorithm,
 * each with different gen_num, and verifies if the winner(s) are indeed as
 * expected.
 */
int varied_n_gen(const problem::base& prob, unsigned int n_winner)
{
	std::vector<algorithm::base_ptr> algos;
	unsigned int gen_interval = 50;
	unsigned int num_instances = 5;

	std::cout << "Testing on problem " << prob.get_name() << ", n_winner = " << n_winner << std::endl;

	for(unsigned int i = 1; i <= num_instances - 1; i++){
		algos.push_back(algorithm::base_ptr(new algorithm::pso_generational(i * gen_interval, 0.7298, 2.05, 2.05, 0.5, 1, 2, 4)));
	}
	// Create a super-algo so that when n_winner == 1 test passes faster
	algos.push_back(algorithm::base_ptr(new algorithm::pso_generational(500, 0.7298, 2.05, 2.05, 0.5, 1, 2, 4)));

	unsigned int pop_size = 50;

	util::racing::race_algo race_dev(algos, prob, pop_size);
	
	std::pair<std::vector<unsigned int>, unsigned int> res = race_dev.run(n_winner, 1, 500, 0.05, std::vector<unsigned int>(), true, true);
	std::vector<unsigned int> winners = res.first;	

	if(winners.size() != n_winner){
		std::cout << "\tError in size of winner list!" << std::endl;
		return 1;
	}

	// Check if all the expected winners are in place

	std::vector<bool> found(num_instances, false);

	for(unsigned int i = 0; i < n_winner; i++){
		if(winners[i] >= num_instances){
			std::cout << "\tError in winner index!" << std::endl;
			return 1;
		}
		found[winners[i]] = true;
	}

	// Only those with large gen_num should win
	for(unsigned int i = 0; i < num_instances; i++){
		if(i < (num_instances - n_winner) && found[i]){
			std::cout << "\tAlgo[" << i << "] found to be winner, but it is not!" << std::endl;
			return 1;
		}
		if(i >= (num_instances - n_winner) && !found[i]){
			std::cout << "\tAlgo[" << i << "] found not to be winner, but it is!" << std::endl;
			return 1;
		}
	}

	std::cout << "Test passed [varied_n_gen]" << std::endl;
	return 0;
}

// Same as above but test on multiple problems
int varied_n_gen(const std::vector<problem::base_ptr> &probs, unsigned int n_winner)
{
	std::vector<algorithm::base_ptr> algos;
	unsigned int gen_interval = 50;
	unsigned int num_instances = 5;

	std::cout << "Testing on multiple problems  (" << probs.size() << " of them), n_winner = " << n_winner << std::endl;

	for(unsigned int i = 1; i <= num_instances - 1; i++){
		algos.push_back(algorithm::base_ptr(new algorithm::pso_generational(i * gen_interval, 0.7298, 2.05, 2.05, 0.5, 1, 2, 4)));
	}
	// Create a super-algo so that when n_winner == 1 test passes faster
	algos.push_back(algorithm::base_ptr(new algorithm::pso_generational(500, 0.7298, 2.05, 2.05, 0.5, 1, 2, 4)));

	unsigned int pop_size = 50;

	util::racing::race_algo race_dev(algos, probs, pop_size);
	
	std::pair<std::vector<unsigned int>, unsigned int> res = race_dev.run(n_winner, 1, 500, 0.05, std::vector<unsigned int>(), true, true);
	std::vector<unsigned int> winners = res.first;	

	if(winners.size() != n_winner){
		std::cout << "\tError in size of winner list!" << std::endl;
		return 1;
	}

	// Check if all the expected winners are in place

	std::vector<bool> found(num_instances, false);

	for(unsigned int i = 0; i < n_winner; i++){
		if(winners[i] >= num_instances){
			std::cout << "\tError in winner index!" << std::endl;
			return 1;
		}
		found[winners[i]] = true;
	}

	// Only those with large gen_num should win
	for(unsigned int i = 0; i < num_instances; i++){
		if(i < (num_instances - n_winner) && found[i]){
			std::cout << "\tAlgo[" << i << "] found to be winner, but it is not!" << std::endl;
			return 1;
		}
		if(i >= (num_instances - n_winner) && !found[i]){
			std::cout << "\tAlgo[" << i << "] found not to be winner, but it is!" << std::endl;
			return 1;
		}
	}

	std::cout << "Test passed [varied_n_gen on multiple problems]" << std::endl;
	return 0;
}

// Test that race_algo accept problems with different constraint vectors via
// the internal dummy constraint padding mechanism.
//
// NOTE: Still needs to command some constraint handling algorithms here to
// verify the results of race.
int test_heterogeneous_constraints()
{
	// Initialize a list of problems with different constraint dimension
	std::vector<problem::base_ptr> probs;
	probs.push_back(problem::cec2006(1).clone());
	probs.push_back(problem::cec2006(2).clone());
	probs.push_back(problem::cec2006(3).clone());
	
	// Set up some algorithms
	std::vector<algorithm::base_ptr> algos;
	unsigned int gen_interval = 25;
	unsigned int num_instances = 5;
	for(unsigned int i = 1; i <= num_instances - 1; i++){
		algos.push_back(algorithm::ihs(i * gen_interval).clone());
	}
	// Create a super-algo so that when n_winner == 1 test passes faster
	algos.push_back(algorithm::ihs(500).clone());

	unsigned int pop_size = 50;

	util::racing::race_algo race_dev(algos, probs, pop_size);

	unsigned int n_winner = 1;

	std::pair<std::vector<unsigned int>, unsigned int> res = race_dev.run(n_winner, 1, 500, 0.05, std::vector<unsigned int>(), true, true);
	std::vector<unsigned int> winners = res.first;	

	if(winners.size() != n_winner){
		std::cout << "\tError in size of winner list!" << std::endl;
		return 1;
	}

	// Check if all the expected winners are in place

	std::vector<bool> found(num_instances, false);

	for(unsigned int i = 0; i < n_winner; i++){
		if(winners[i] >= num_instances){
			std::cout << "\tError in winner index!" << std::endl;
			return 1;
		}
		found[winners[i]] = true;
	}

	// Only those with large gen_num should win
	for(unsigned int i = 0; i < num_instances; i++){
		if(i < (num_instances - n_winner) && found[i]){
			std::cout << "\tAlgo[" << i << "] found to be winner, but it is not!" << std::endl;
			return 1;
		}
		if(i >= (num_instances - n_winner) && !found[i]){
			std::cout << "\tAlgo[" << i << "] found not to be winner, but it is!" << std::endl;
			return 1;
		}
	}

	std::cout << "Test passed [heterogeneous constraints]" << std::endl;

	return 0;
}

/*
// TODO: Find out offline which variant works best and verify in this test?
int varied_pso_variant(const problem::base_ptr& prob)
{
	std::vector<algorithm::base_ptr> algos;
	for(int i = 1; i <= 6; i++){
		algos.push_back(algorithm::base_ptr(new algorithm::pso_generational(100, 0.7298, 2.05, 2.05, 0.5, i, 2, 4)));
	}
	algorithm::race race_obj(*prob, algos);
	std::vector<algorithm::race::size_type> winners = race_obj.run(1);
	return 0;
}
*/

int main()
{
	int dimension = 10;
	problem::ackley prob(dimension);
	problem::ackley prob2(dimension*2);
	problem::ackley prob3(dimension*3);
	std::vector<problem::base_ptr> prob_list;	
	prob_list.push_back(prob.clone());
	prob_list.push_back(prob2.clone());
	prob_list.push_back(prob3.clone());

	return
		varied_n_gen(prob, 1) || 
		varied_n_gen(prob, 2) ||
		varied_n_gen(prob_list, 1) ||
		varied_n_gen(prob_list, 2) ||
		test_heterogeneous_constraints();
}
