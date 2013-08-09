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

#include <iostream>
#include <numeric>
#include <iomanip>
#include <limits>
#include <string>
#include "../src/pagmo.h"

using namespace pagmo;

int archi_best_idx(archipelago archi) {
	double min = archi.get_island(0)->get_population().champion().f[0];
	int idx=0;
	for (size_t i=1;i<archi.get_size();++i) {
		double cur = archi.get_island(i)->get_population().champion().f[0];
		if (cur < min) {
			min=cur;
			idx=i;
		}
	}
	return idx;
}


double run_experiment(const int n_isl, const int pop_size, const int n_gen, bool use_racing, int seed=123)
{
	// Amount of evaluation budget derived from conventional set-up
	const int max_fevals = n_gen * pop_size;

	// Number of splits to distribute the total fevals budget (reason is that
	// with racing number of fevals consumed in each generation will be
	// different. Ideally, in each split two variants of algorithms (with or
	// without racing) should consume the same number of fevals.
	const int n_split = 1;

	const int max_fevals_per_split = max_fevals / n_split;

	const int prob_n_eval_racing = 5;
	const int prob_n_eval = 5;

	algorithm::base_ptr algo;
	if(use_racing){
		std::cout << "Using racing: " << std::endl;
		//std::cout << "\t-> Splitting into " << n_split << " times of evolves, each consuming " << max_fevals_per_split << " fevals." << std::endl;
		//std::cout << "\t-> Total fevals = " << n_split * max_fevals_per_split << std::endl;
		algo = algorithm::pso_generational(std::numeric_limits<int>::max(), 0.7298, 2.05, 2.05, 0.05, 5, 2, 4, true, max_fevals_per_split).clone();
		//algo = algorithm::pso_generational(n_gen / n_split / 2, 0.7298, 2.05, 2.05, 0.05, 5, 2, 4, true).clone();
	}
	else{
		std::cout << "Not using racing: " << std::endl;
		//std::cout << "\t-> Splitting into " << n_split << " times of evolves, each consuming " << (n_gen / n_split) * pop_size << " fevals." << std::endl;
		//std::cout << "\t-> Total fevals = " << n_gen * pop_size << std::endl;
		algo = algorithm::pso_generational(n_gen / n_split / 2, 0.7298, 2.05, 2.05, 0.05).clone();
	}
	algo->reset_rngs(seed);

	std::cout << "Initializing ..." << std::endl;
	archipelago archi = archipelago(topology::fully_connected());

	for (int j = 0; j < n_isl; ++j) {

		problem::base_ptr p_prob;
		if(use_racing){
			p_prob = problem::spheres(prob_n_eval_racing, 10, 1e-6, seed, false).clone();
		}
		else{
			p_prob = problem::spheres(prob_n_eval, 10, 1e-6, seed, false).clone();
		}
		//problem::spheres prob(prob_n_eval, 10, 1e-6, rand(), false);
		// This instantiates a population within the original bounds (-1,1)
		population pop_temp(*p_prob,pop_size, seed);

		// We make the bounds larger to allow neurons weights to grow
		p_prob->set_bounds(-10,10);

		// We create an empty population on the new prolem (-10,10)
		population pop(*p_prob, 0, seed);

		// And we fill it up with (-1,1) individuals having zero velocities
		decision_vector v(p_prob->get_dimension(),0);
		for (int i =0; i<pop_size; ++i) {
			pop.push_back(pop_temp.get_individual(i).cur_x);
			pop.set_v(i,v);
		}

		algo->reset_rngs(seed + j);
		archi.push_back(island(*algo,pop));
	}

	std::vector<double> record_best_f;
	//Evolution is here started on the archipelago
	for (int i = 0; i< n_split; ++i){

		int idx = archi_best_idx(archi);
		//std::cout << "best so far after " << i * max_fevals_per_split << " evaluations on each island ......" << "\n" << archi.get_island(idx)->get_population().champion().x << std::endl;
		double best_f = archi.get_island(idx)->get_population().champion().f[0];
		record_best_f.push_back(best_f);
		std::cout << "fevals "<< std::setw(12) << i * max_fevals_per_split << std::setw(12) << best_f << std::endl;
		/*
		double mean = 0.0;
		mean = std::accumulate(buff.begin(),buff.end(),mean);
		mean /= (double)buff.size();
		std::cout << "gen: "<< std::setw(12) << i << std::setw(12) <<
		best_f << std::setw(12) <<
		archi.get_island(idx)->get_population().mean_velocity() << std::setw(12) <<
		mean <<	 std::endl;
		*/
		archi.evolve(1);
	 }

	int idx = archi_best_idx(archi);
	decision_vector winner_x = archi.get_island(idx)->get_population().champion().x;
	//std::cout << "and the winner is ......" << "\n" << archi.get_island(idx)->get_population().champion().x << std::endl;
	/*
	double final_f = archi.get_island(idx)->get_population().champion().f[0];
	std::cout << "Final f = " << final_f << std::endl;
	*/	
	int eval_count = 20000;
	problem::spheres prob_for_evaluation(eval_count, 10, 1e-6, rand(), false);
	double final_f = prob_for_evaluation.objfun(winner_x)[0];
	std::cout << "fevals "<< std::setw(12) << max_fevals_per_split << std::setw(12) << final_f << std::endl;

	return final_f;
}

int main(int argc, char* argv[])
{	
	/*
	// EXPERIMENT SET-UP //
	const int n_isl = 8;
	const int pop_size = 512;
	const int n_eval = 5;
	const int n_gen = 400;
	const bool use_racing = false;
	// END OF EXPERIMENT SET-UP //
	*/	
	const int n_isl = 1;
	const int pop_size = 20;
	const int n_gen = 400;
	const int n_repeat_experiments = 10;
	bool use_racing = false;

	if(argc > 1){
		if(argv[1][0] == '1'){
			use_racing = true;
		}
	}

	unsigned int seed = 5;

	double avg_final_f = 0;
	for(int i = 0; i < n_repeat_experiments; i++, seed+=123){
		double final_f = run_experiment(n_isl, pop_size, n_gen, use_racing, seed);
		avg_final_f += final_f / (double)n_repeat_experiments;
	}
	
	std::cout << "Final averaged fitness of the champions over " << n_repeat_experiments << " independent runs = " << avg_final_f << std::endl;

	return 0;
}

int main_orig()
{
	// EXPERIMENT SET-UP //
	const int n_isl = 8;
	const int pop_size = 512;
	const int n_eval = 5;
	const int n_gen = 400;
	// END OF EXPERIMENT SET-UP //

	// Buffer
	std::vector<double> buff;
	// We instantiate a PSO algorithm capable of coping with stochastic prolems
	algorithm::pso_generational algo(1,0.7298,2.05,2.05,0.05);

	// This instantiates the spheres problem
	std::cout << "Initializing ...." << std::endl;

	archipelago archi = archipelago(topology::fully_connected());

	for (int j=0;j<n_isl; ++j) {

		problem::spheres prob(n_eval,10,1e-6,rand(), false);
		// This instantiates a population within the original bounds (-1,1)
		population pop_temp(prob,pop_size);

		// We make the bounds larger to allow neurons weights to grow
		prob.set_bounds(-10,10);

		// We create an empty population on the new prolem (-10,10)
		population pop(prob);


		// And we fill it up with (-1,1) individuals having zero velocities
		decision_vector v(prob.get_dimension(),0);
		for (int i =0; i<pop_size; ++i) {
			pop.push_back(pop_temp.get_individual(i).cur_x);
			pop.set_v(i,v);
		}
		archi.push_back(island(algo,pop));
	}

	//Evolution is here started on the archipelago
	for (int i=0; i< n_gen; ++i){
		int idx = archi_best_idx(archi);
		if (!(i%100)) {
			std::cout << "best so far ......" << "\n" << archi.get_island(idx)->get_population().champion().x << std::endl;
		}
		double best_f = archi.get_island(idx)->get_population().champion().f[0];

		if (i<50) {
			 buff.push_back(best_f);
		}
		else {
			 (buff[i%50] = best_f);
		}
		double mean = 0.0;
		mean = std::accumulate(buff.begin(),buff.end(),mean);
		mean /= (double)buff.size();
		std::cout << "gen: "<< std::setw(12) << i << std::setw(12) <<
		best_f << std::setw(12) <<
		archi.get_island(idx)->get_population().mean_velocity() << std::setw(12) <<
		mean <<	 std::endl;
		archi.evolve(1);
	 }

	int idx = archi_best_idx(archi);
	std::cout << "and the winner is ......" << "\n" << archi.get_island(idx)->get_population().champion().x << std::endl;

	return 0;
}
