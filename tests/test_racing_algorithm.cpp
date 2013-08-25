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

// TODO: Test is too slow!
int main()
{
	int dimension = 10;
	problem::ackley prob(dimension);
	return
		varied_n_gen(prob, 1) || 
		varied_n_gen(prob, 2);
}
