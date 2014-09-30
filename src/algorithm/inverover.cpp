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

#include "inverover.h"

namespace pagmo { namespace algorithm {
    
/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] gen Number of generations to evolve.
 * @param[in] ri Probability of performing a random invert (mutation probability)
*/

inverover::inverover(int gen, double ri)
	:base(),m_gen(gen),m_ri(ri)
{
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if (ri > 1 || ri < 0) {
		pagmo_throw(value_error,"random invert probability must be in the [0,1] range");
	}

}

    
    /// Clone method.
    base_ptr inverover::clone() const
    {
	return base_ptr(new inverover(*this));
    }
    
    /// Evolve implementation.
    /**
     * Runs the ACO algorithm for the number of generations specified in the constructor.
     *
     * @param[in,out] pop input/output pagmo::population to be evolved.
     */
    void inverover::evolve(population &pop) const
    {

	

        // Let's store some useful variables.
        const problem::tsp &tsp_prob = dynamic_cast<const problem::tsp &>(pop.problem());
        const std::vector<std::vector<double> > &weights = tsp_prob.get_weights();
        const problem::base::size_type Nv = tsp_prob.get_n_vertices();
        const problem::base::size_type NP = pop.size();
	
        //create proability distributions
        std::default_random_engine generator(time(NULL));
	std::uniform_int_distribution<int> unif_NPless1(0, NP - 2);
	std::uniform_int_distribution<int> unif_Nv(0, Nv - 1);
	std::uniform_int_distribution<int> unif_Nvless1(0, Nv - 2);
	std::uniform_real_distribution<double> unif_01(0,1);


        //TODO: add check on the problem. The problem needs to be an integer problem TSP like, single objective.
        if (NP <= 1) pagmo_throw(value_error, "population size must be greater than one.");
          
	//create population
	std::vector<std::vector<int> > my_pop(NP, std::vector<int>(Nv));
	//vector that stores the length of each individual - ONLY necessary for the atomatic stop criterion
	std::vector<double>  fitness(NP, 0);
	double best_fitness = 0;

	//initialize with nearest neighbor algorithm
	/*
	int nxt_city;
	size_t min_idx;
	std::vector<int> not_visited(Nv);
	for (size_t i = 0; i < NP; i++) {
		for (size_t j = 0; j < Nv; j++) {
			not_visited[j] = j;
		}
		my_pop[i][0] = unif_Nv(generator);
		std::swap(not_visited[my_pop[i][0]],not_visited[Nv-1]);
		for (size_t j = 1; j < Nv-1; j++) {
			min_idx = 0;
			nxt_city = not_visited[0];
			for (size_t l = 1; l < Nv-j; l++) {
				if(weights[my_pop[i][j-1]][not_visited[l]] < weights[my_pop[i][j-1]][nxt_city]){
					min_idx = l;		
					nxt_city = not_visited[l];}
			}
			my_pop[i][j] = nxt_city;
			std::swap(not_visited[min_idx],not_visited[Nv-j-1]);
		}
		my_pop[i][Nv-1] = not_visited[0];
	}
	
	*/

	//random initialization (produces feasible individuals)
	for (size_t i = 0; i < NP; i++) {
		for (size_t j = 0; j < Nv; j++) {
		my_pop[i][j] = j;
		}
	}
	int tmp;
	size_t rnd_idx;
	for (size_t j = 1; j < Nv-1; j++) {
	        std::uniform_int_distribution<int> dist(j, Nv - 1);
		for (size_t i = 0; i < NP; i++) {
		rnd_idx = dist(generator);
		tmp = my_pop[i][j];
		my_pop[i][j] = my_pop[i][rnd_idx];
		my_pop[i][rnd_idx] = tmp;
		}	
	}
	
	//compute fitness
	for(size_t i=0; i < NP; i++){
    		fitness[i] = weights[my_pop[i][Nv-1]][my_pop[i][0]];
    		for(size_t k=0; k < Nv; k++){
        		fitness[i] = fitness[i] + weights[my_pop[i][k-1]][my_pop[i][k]];
    		}
		if(i==0 || fitness[i] < best_fitness){
			best_fitness = fitness[i];
		}
	}
	
	std::vector<int> tmp_tour(Nv);
	bool stop, stop_algo;
	size_t last_change, rnd_num, i2, pos1_c1, pos1_c2, pos2_c1, pos2_c2; //pos2_c1 denotes the position of city1 in parent2
	double fitness_change;
	int total_iter;

	stop_algo = false;
	total_iter = 0;
	last_change = 0;

	//InverOver main loop
	while(!stop_algo && total_iter < m_gen){
		for(size_t i1 = 0; i1 < NP; i1++){
			fitness_change = 0;
			tmp_tour = my_pop[i1];
			pos1_c1 = unif_Nv(generator);
			stop = false;
			while(!stop){
				if(unif_01(generator) < m_ri){
					rnd_num = unif_Nvless1(generator);
					pos1_c2 = (rnd_num == pos1_c1? Nv-1:rnd_num);
				}
				else{
					i2 = unif_NPless1(generator);
					i2 = (i2 == i1? NP-1:i2);
					pos2_c1 = std::find(my_pop[i2].begin(),my_pop[i2].end(),tmp_tour[pos1_c1])-my_pop[i2].begin();
					pos2_c2 = (pos2_c1 == Nv-1? 0:pos2_c1+1);
					pos1_c2 = std::find(tmp_tour.begin(),tmp_tour.end(),my_pop[i2][pos2_c2])-tmp_tour.begin();
				}
				stop = (abs(pos1_c1-pos1_c2)==1 || abs(pos1_c1-pos1_c2)==Nv-1);
				if(!stop){
					
					if(pos1_c1<pos1_c2){
						fitness_change -= weights[tmp_tour[pos1_c1]][tmp_tour[pos1_c1+1]] + weights[tmp_tour[pos1_c2]][tmp_tour[pos1_c2+1 - (pos1_c2+1 > Nv-1? Nv:0)]];
						for(size_t l=0; l < (double (pos1_c2-pos1_c1-1)/2); l++){
							std::swap(tmp_tour[pos1_c1+1+l],tmp_tour[pos1_c2-l]);}
						fitness_change += weights[tmp_tour[pos1_c1]][tmp_tour[pos1_c1+1]] + weights[tmp_tour[pos1_c2]][tmp_tour[pos1_c2+1 - (pos1_c2+1 > Nv-1? Nv:0)]];
					}
					else{
						/*
						//Methode1: invert section from c2 to c1
						fitness_change -= weights[tmp_tour[pos1_c1-1]][tmp_tour[pos1_c1]] + weights[tmp_tour[pos1_c2-1]][tmp_tour[pos1_c2]];
						std::reverse(&tmp_tour[pos1_c2],&tmp_tour[pos1_c1]);
						fitness_change += weights[tmp_tour[pos1_c1-1]][tmp_tour[pos1_c1]] + weights[tmp_tour[pos1_c2-1]][tmp_tour[pos1_c2]];
						*/

						//Methode2: invert section from c1 to c2 (performs better than Methode 1 - shorter tour in less time)
						fitness_change -= weights[tmp_tour[pos1_c1]][tmp_tour[pos1_c1+1 - (pos1_c1+1 > Nv-1? Nv:0)]] + weights[tmp_tour[pos1_c2]][tmp_tour[pos1_c2+1]];
						for(size_t l=0; l < (double (Nv-(pos1_c1-pos1_c2)-1)/2); l++){
							std::swap(tmp_tour[pos1_c1+1+l - (pos1_c1+1+l>Nv-1? Nv:0)],tmp_tour[pos1_c2-l + (pos1_c2<l? Nv:0)]);}
						fitness_change += weights[tmp_tour[pos1_c1]][tmp_tour[pos1_c1+1 - (pos1_c1+1 > Nv-1? Nv:0)]] + weights[tmp_tour[pos1_c2]][tmp_tour[pos1_c2+1]];
						
					}
					//pos1_c1 = (pos1_c1 == Nv-1? 0:pos1_c1+1); //original Inver-Over
					pos1_c1 = pos1_c2; //better performance than original Inver-Over (shorter tour in less time)
				}
			} //end of while loop (looping over a single indvidual)	
			if(fitness_change < 0){ //replace individual?
				my_pop[i1] = tmp_tour;
				fitness[i1] += fitness_change;
				if(fitness[i1] < best_fitness){
					best_fitness = fitness[i1];
					last_change = total_iter;
				}
			}

		} //end of loop over population
		if(total_iter > 10 && total_iter-last_change > last_change){ //atomatic stop criterion if no change is found after a certain amount of iterations
			stop_algo = true;
		}
		total_iter++;
	} //end of loop over generations


	//change representation of tour
    	for (size_t i = 0; i < NP; i++) {
		decision_vector individual(Nv*(Nv-1), 0);
		for (size_t j = 0; j < Nv; j++) {
			individual[(my_pop[i][j])*(Nv-1)+my_pop[i][(j+1>Nv-1? 0:j+1)] - (my_pop[i][(j+1>Nv-1? 0:j+1)]>my_pop[i][j]? 1:0)] = 1;
		}
		pop.set_x(i,individual);
	}

    } // end of evolve
    

	
    /// Algorithm name
    std::string inverover::get_name() const
    {
        return "InverOver Algorithm";
    }

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::inverover)
