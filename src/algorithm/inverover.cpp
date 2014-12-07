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

#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include "../config.h"
#include "../serialization.h"
#include "../population.h"
#include "../problem/base_tsp.h"
#include "../algorithm/nn_tsp.h"
#include "base.h"
#include "inverover.h"

namespace pagmo { namespace algorithm {
    
/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] gen Number of generations to evolve.
 * @param[in] ri Probability of performing a random invert (mutation probability)
*/

inverover::inverover(int gen, double ri, initialization_type ini_type)
	:base(),m_gen(gen),m_ri(ri),m_ini_type(ini_type)
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
     * Runs the Inverover algorithm for the number of generations specified in the constructor.
     *
     * @param[in,out] pop input/output pagmo::population to be evolved.
     */
    void inverover::evolve(population &pop) const
    {

	const problem::base_tsp* prob;
	//check if problem is of type pagmo::problem::base_tsp
	try
	{
		const problem::base_tsp& tsp_prob = dynamic_cast<const problem::base_tsp &>(pop.problem());
		prob = &tsp_prob;
	}
	catch (const std::bad_cast& e)
	{
		pagmo_throw(value_error,"Problem not of type pagmo::problem::base_tsp");
	}
	
	// Let's store some useful variables.

	const population::size_type NP = pop.size();
	const problem::base::size_type Nv = prob->get_n_cities();

	// Get out if there is nothing to do.
	if (m_gen == 0) {
		return;
	}
	
	// Initializing the random number generators
	boost::uniform_real<double> uniform(0.0,1.0);
	boost::variate_generator<boost::lagged_fibonacci607 &, boost::uniform_real<double> > unif_01(m_drng,uniform);
	boost::uniform_int<int> NPless1(0, NP - 2);
	boost::variate_generator<boost::mt19937 &, boost::uniform_int<int> > unif_NPless1(m_urng,NPless1);
	boost::uniform_int<int> Nv_(0, Nv - 1);
	boost::variate_generator<boost::mt19937 &, boost::uniform_int<int> > unif_Nv(m_urng,Nv_);
	boost::uniform_int<int> Nvless1(0, Nv - 2);
	boost::variate_generator<boost::mt19937 &, boost::uniform_int<int> > unif_Nvless1(m_urng,Nvless1);
	
	//create own local population
	std::vector<decision_vector> my_pop(NP, decision_vector(Nv));

	//check if some individuals in the population that is passed as a function input are feasible.
	bool feasible;
	std::vector<int> not_feasible;
	for (size_t i = 0; i < NP; i++) {
		feasible = prob->feasibility_x(pop.get_individual(i).cur_x);
		if(feasible){ //if feasible store it in my_pop
			switch( prob->get_encoding() ) {
			    case problem::base_tsp::FULL:
			        my_pop[i] = prob->full2cities(pop.get_individual(i).cur_x);
			        break;
			    case problem::base_tsp::RANDOMKEYS:
			        my_pop[i] = prob->randomkeys2cities(pop.get_individual(i).cur_x);
			        break;
			    case problem::base_tsp::CITIES:
			        my_pop[i] = pop.get_individual(i).cur_x;
			        break;
			}
		}
		else
		{
			not_feasible.push_back(i);
		}
	}

	//replace the not feasible individuals by feasible ones	
	int i;		
	switch (m_ini_type){
		case 0:
		{
		//random initialization (produces feasible individuals)
			for (size_t ii = 0; ii < not_feasible.size(); ii++) {
				i = not_feasible[ii];
				for (size_t j = 0; j < Nv; j++) {
					my_pop[i][j] = j;
				}
			}
			int tmp;
			size_t rnd_idx;
			for (size_t j = 1; j < Nv-1; j++) {
	        		boost::uniform_int<int> dist_(j, Nv - 1);
					boost::variate_generator<boost::mt19937 &, boost::uniform_int<int> > dist(m_urng,dist_);
					
				for (size_t ii = 0; ii < not_feasible.size(); ii++) {
					i = not_feasible[ii];
					rnd_idx = dist();
					tmp = my_pop[i][j];
					my_pop[i][j] = my_pop[i][rnd_idx];
					my_pop[i][rnd_idx] = tmp;
				}	

			}
			break;
		}
		case 1:
		{
		//initialize with nearest neighbor algorithm
		std::vector<int> starting_notes(std::max(Nv,not_feasible.size()));
			for (size_t j = 0; j < starting_notes.size(); j++) {
					starting_notes[j] = j;
			}
			//std::shuffle(starting_notes.begin(), starting_notes.end(), m_urng);
			for (size_t ii = 0; ii < not_feasible.size(); ii++) {
				i = not_feasible[ii];
				pagmo::population one_ind_pop(pop.problem(), 1);
				pagmo::algorithm::nn_tsp algo(starting_notes[i] % Nv);
				algo.evolve(one_ind_pop);
				my_pop[i] = one_ind_pop.get_individual(0).cur_x;
			}
			break;
		}
		default:
			pagmo_throw(value_error,"Invalid initialization type");
	}	
	
	std::vector<fitness_vector>  fitness(NP, fitness_vector(1));
	for(size_t i=0; i < NP; i++){
    		fitness[i] = prob->objfun(my_pop[i]);
	}
	
	decision_vector tmp_tour(Nv);
	bool stop, changed;
	size_t rnd_num, i2, pos1_c1, pos1_c2, pos2_c1, pos2_c2; //pos2_c1 denotes the position of city1 in parent2
	fitness_vector fitness_tmp;

	//InverOver main loop
	for(int iter = 0; iter < m_gen; iter++){
		for(size_t i1 = 0; i1 < NP; i1++){
			tmp_tour = my_pop[i1];
			pos1_c1 = unif_Nv();
			stop = false;
			changed = false;
			while(!stop){
				if(unif_01() < m_ri){
					rnd_num = unif_Nvless1();
					pos1_c2 = (rnd_num == pos1_c1? Nv-1:rnd_num);
				}
				else{
					i2 = unif_NPless1();
					i2 = (i2 == i1? NP-1:i2);
					pos2_c1 = std::find(my_pop[i2].begin(),my_pop[i2].end(),tmp_tour[pos1_c1])-my_pop[i2].begin();
					pos2_c2 = (pos2_c1 == Nv-1? 0:pos2_c1+1);
					pos1_c2 = std::find(tmp_tour.begin(),tmp_tour.end(),my_pop[i2][pos2_c2])-tmp_tour.begin();
				}
				stop = (_abs64(pos1_c1-pos1_c2)==1 || _abs64(pos1_c1-pos1_c2)==Nv-1);
				if(!stop){
					changed = true;
					if(pos1_c1<pos1_c2){
						for(size_t l=0; l < (double (pos1_c2-pos1_c1-1)/2); l++){
							std::swap(tmp_tour[pos1_c1+1+l],tmp_tour[pos1_c2-l]);}
					}
					else{
						//inverts the section from c1 to c2 (see documentation Note3)
						
						for(size_t l=0; l < (double (Nv-(pos1_c1-pos1_c2)-1)/2); l++){
							std::swap(tmp_tour[pos1_c1+1+l - (pos1_c1+1+l>Nv-1? Nv:0)],tmp_tour[pos1_c2-l + (pos1_c2<l? Nv:0)]);}
						
					}
					pos1_c1 = pos1_c2; //better performance than original Inver-Over (shorter tour in less time)
				}
			} //end of while loop (looping over a single indvidual)
			if(changed){
				fitness_tmp = prob->objfun(tmp_tour);
				if(prob->compare_fitness(fitness_tmp,fitness[i1])){ //replace individual?
					my_pop[i1] = tmp_tour;
					fitness[i1][0] = fitness_tmp[0];
				}
			}
		} //end of loop over population
	} //end of loop over generations


	//change representation of tour
    	for (size_t ii = 0; ii < NP; ii++) {
			switch( prob->get_encoding() ) {
			    case problem::base_tsp::FULL:
			        pop.set_x(ii,prob->cities2full(my_pop[ii]));
			        break;
			    case problem::base_tsp::RANDOMKEYS:
			        pop.set_x(ii,prob->cities2randomkeys(my_pop[ii],pop.get_individual(ii).cur_x));
			        break;
			    case problem::base_tsp::CITIES:
			        pop.set_x(ii,my_pop[ii]);
			        break;
			}
		}

    } // end of evolve
    

	
    /// Algorithm name
    std::string inverover::get_name() const
    {
        return "InverOver Algorithm";
    }

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::inverover)
