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
#include <fstream>

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

double post_evaluation(const problem::spheres& prob, int eval_count, const decision_vector& x, unsigned int seed=0)
{
	// Averaging method 1: Same seed, loop many times internally per objfun call
	//problem::spheres prob_for_(eval_count, 10, 1e-6, rand(), false);
	//double final_f = prob_for_evaluation.objfun(winner_x)[0];

	// Averaging method 2: Different seeds, evaluate one time per objfun call
	rng_uint32 urng(seed);
	double post_evaluation_f = 0;
	for(int i = 0; i < eval_count; i++){
		prob.set_seed(urng());
		post_evaluation_f += prob.objfun(x)[0] / (double)eval_count;
	}
	return post_evaluation_f;
}


// Evovle population up to for n_gen generations in one go. Evaluate the final
// champion based on its converged fitness value, obtained from a problem with
// a very large n_eval.
double evolve_in_one_go(int n_isl, int pop_size, int nr_eval_per_x, int n_gen, bool use_racing, int seed=123)
{

	algorithm::base_ptr algo;

	// Here as we can evolve in one go, we can set a maximum cap for the
	// allowed number of objective function calls. This number is based on the
	// behaviour of non-racing pso_gen.
	unsigned int max_fevals = n_gen * pop_size * 2 * nr_eval_per_x + pop_size * nr_eval_per_x;

	if(use_racing){
		std::cout << "---Using racing---" << std::endl;
		algo = algorithm::pso_generational_racing(std::numeric_limits<int>::max(), 0.7298, 2.05, 2.05, 0.05, 5, 2, 4, nr_eval_per_x, max_fevals).clone();
	}
	else{
		std::cout << "---Not using racing---" << std::endl;
		algo = algorithm::pso_generational(std::numeric_limits<int>::max(), 0.7298, 2.05, 2.05, 0.05, 5, 2, 4).clone();
	}
	algo->reset_rngs(seed);

	std::cout << "Initializing ..." << std::endl;
	archipelago archi = archipelago(topology::fully_connected());

	for (int j = 0; j < n_isl; ++j) {

		problem::spheres prob(1, 10, 1e-6, seed, false);

		// This instantiates a population within the original bounds (-1,1)
		population pop_temp(prob, pop_size, seed);

		// We make the bounds larger to allow neurons weights to grow
		prob.set_bounds(-10,10);

		// We create an empty population on the new prolem (-10,10)
		population pop(prob, 0, seed);

		// And we fill it up with (-1,1) individuals having zero velocities
		decision_vector v(prob.get_dimension(),0);
		for (int i =0; i<pop_size; ++i) {
			pop.push_back(pop_temp.get_individual(i).cur_x);
			pop.set_v(i,v);
		}

		algo->reset_rngs(seed + j);
		archi.push_back(island(*algo,pop));
	}

	std::vector<double> record_best_f;
	//Evolution is here started on the archipelago
	std::cout << "Start evolving for " << n_gen << " generations..." << std::endl;
	archi.evolve(1);
	int idx = archi_best_idx(archi);
	decision_vector winner_x = archi.get_island(idx)->get_population().champion().x;
	std::cout << "Evolution done." << std::endl;

	int eval_count = 50;

	problem::spheres prob_for_evaluation(1, 10, 1e-6, seed, false);
	double post_evaluation_f = post_evaluation(prob_for_evaluation, eval_count, winner_x, seed);
	return post_evaluation_f;
}

double evolve_original(int n_isl, int pop_size, int n_eval_per_x, int n_gen, bool use_racing, unsigned int seed = 123)
{
	// Buffer
	std::vector<double> buff;

	int gen_batch_size = 10;

	// We instantiate a PSO algorithm capable of coping with stochastic prolems
	algorithm::base_ptr algo_ptr;
	int pso_variant = 5;
	int pso_neighb_type = 2;
	int pso_neighb_param = 4;

	// The implicit maximum evaluation budget based on the parameters
	unsigned int max_fevals = gen_batch_size * pop_size * 2 * n_eval_per_x + pop_size * n_eval_per_x;

	if(use_racing){
		algo_ptr = algorithm::pso_generational_racing(gen_batch_size, 0.7298, 2.05, 2.05, 0.05, pso_variant, pso_neighb_type, pso_neighb_param, n_eval_per_x, max_fevals).clone();
	}
	else{
		algo_ptr = algorithm::pso_generational(gen_batch_size, 0.7298, 2.05, 2.05, 0.05, pso_variant, pso_neighb_type, pso_neighb_param).clone();
	}

	algorithm::base &algo = *algo_ptr;

	// This instantiates the spheres problem
	std::cout << "Initializing ...." << std::endl;

	archipelago archi = archipelago(topology::fully_connected());


	for (int j=0;j<n_isl; ++j) {

		problem::base_ptr prob_ptr;
		if(use_racing){
			prob_ptr = problem::spheres(1, 10, 1e-6, seed, false).clone();
		}
		else{
			prob_ptr = problem::spheres(n_eval_per_x, 10, 1e-6, seed, false).clone();
		}
		// This instantiates a population within the original bounds (-1,1)
		population pop_temp(*prob_ptr, pop_size, seed);

		// We make the bounds larger to allow neurons weights to grow
		prob_ptr->set_bounds(-10,10);

		// We create an empty population on the new prolem (-10,10)
		population pop(*prob_ptr);

		// And we fill it up with (-1,1) individuals having zero velocities
		decision_vector v(prob_ptr->get_dimension(),0);
		for (int i =0; i<pop_size; ++i) {
			pop.push_back(pop_temp.get_individual(i).cur_x);
			pop.set_v(i,v);
		}
		algo.reset_rngs(seed + 17 * j);
		archi.push_back(island(algo,pop));
	}

	int window_width = 10;
	int post_evaluation_n = 100;

	//Evolution is here started on the archipelago
	problem::spheres prob_eval(1, 10, 1e-6, seed, false);
	double mean = 0.0;
	for (int i=0; i< n_gen; ++i){
		int idx = archi_best_idx(archi);
		if (!(i%1)) {
			//std::cout << "best so far ......" << "\n" << archi.get_island(idx)->get_population().champion().x << std::endl;
			decision_vector best_x = archi.get_island(idx)->get_population().champion().x;
			double post_eval_f = post_evaluation(prob_eval, post_evaluation_n, best_x, seed);
			std::cout << "Post-evaluated current best f = " << post_eval_f << std::endl;
		}
		double best_f = archi.get_island(idx)->get_population().champion().f[0];

		if (i < window_width) {
			 buff.push_back(best_f);
		}
		else {
			 (buff[i%window_width] = best_f);
		}
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

	return mean;
}

// Evolve until max_fevals is hit, repeat for multiple times and average the
// final fitness value
int repeat_experiment(int n_trials, int n_isl, int pop_size, int nr_eval_per_x, int n_gen, bool use_racing, int seed=123)
{
	double averaged_outcome = 0;
	for(int i = 0; i < n_trials; i++, seed+=123){
		//double outcome = evolve_in_one_go(n_isl, pop_size, nr_eval_per_x, n_gen, use_racing, seed);
		double outcome = evolve_original(n_isl, pop_size, nr_eval_per_x, n_gen, use_racing, seed);
		averaged_outcome += outcome / (double)n_trials;
	}
	if(use_racing){
		std::cout << ":::PSO with racing:::" << std::endl;
	}
	else{
		std::cout << ":::PSO without racing:::" << std::endl;
	}
	std::cout << "Results after repeating the experiments for " << n_trials << " independent trials: " << averaged_outcome << std::endl;
	return 0;
}

// Usage:
// To evolve using pso_gen with racing, invoke the program as follows:
// $ examples/evolve_spheres 1
// Otherwise, pso_gen without racing will be used for the evolution.
int main(int argc, char* argv[])
{
	// EXPERIMENT SET-UP //
	const int n_isl = 1;
	const int pop_size = 80;
	const int nr_eval_per_x = 5;
	const int n_gen = 200;
	bool use_racing = false;
	// END OF EXPERIMENT SET-UP //

	if(argc > 1){
		if(argv[1][0] == '1'){
			use_racing = true;
		}
	}

	std::cout << "n_isl = " << n_isl << " ";
	std::cout << "pop_size = " << pop_size << " ";
	std::cout << "nr_eval_per_x = " << nr_eval_per_x << " ";
	std::cout << "n_gen = " << n_gen << " ";
	std::cout << "use_racing = " << use_racing << " ";
	std::cout << std::endl;

	unsigned int seed = 5;

	//run_experiment_original(n_isl, pop_size, nr_eval_per_x, n_gen, use_racing, seed);
	//run_experiment_alternative(n_isl, pop_size, nr_eval_per_x, n_gen, use_racing, seed);
	repeat_experiment(5, n_isl, pop_size, nr_eval_per_x, n_gen, use_racing, seed);

	return 0;
}

