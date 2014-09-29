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

#include "nn.h"

namespace pagmo { namespace algorithm {
    
/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] start_city First City in the tour.
*/

nn::nn(int start_city)
	:base(),m_start_city(start_city)
{
}

    
    /// Clone method.
    base_ptr nn::clone() const
    {
	return base_ptr(new nn(*this));
    }
    
    /// Evolve implementation.
    /**
     * Runs the NN algorithm.
     *
     * @param[in,out] pop input/output pagmo::population to be evolved.
     */
    void nn::evolve(population &pop) const
    {

	

        // Let's store some useful variables.
        const problem::tsp &tsp_prob = dynamic_cast<const problem::tsp &>(pop.problem());
        const std::vector<std::vector<double> > &weights = tsp_prob.get_weights();
        const problem::base::size_type Nv = tsp_prob.get_n_vertices();
	 
	//create individuals
	std::vector<int> best_tour(Nv);
	std::vector<int> new_tour(Nv);
	
	size_t first_city, Nt;
	if(m_start_city == -1){
		first_city = 0;  
      		Nt = Nv;
	}
	else{
		first_city = m_start_city; 
		Nt = m_start_city+1;
	}

	int Lbest_tour, Lnew_tour;
	size_t nxt_city, min_idx;
	std::vector<int> not_visited(Nv);
	Lbest_tour = 0;
	
	//NN main loop
	for (size_t i = first_city; i < Nt; i++) {
		Lnew_tour = 0;
		for (size_t j = 0; j < Nv; j++) {
			not_visited[j] = j;
		}
		new_tour[0] = i;
		std::swap(not_visited[new_tour[0]],not_visited[Nv-1]);
		for (size_t j = 1; j < Nv-1; j++) {
			min_idx = 0;
			nxt_city = not_visited[0];
			for (size_t l = 1; l < Nv-j; l++) {
				if(weights[new_tour[j-1]][not_visited[l]] < weights[new_tour[j-1]][nxt_city]){
					min_idx = l;		
					nxt_city = not_visited[l];}
			}
			new_tour[j] = nxt_city;
			Lnew_tour += weights[new_tour[j-1]][nxt_city];
			std::swap(not_visited[min_idx],not_visited[Nv-j-1]);
		}
		new_tour[Nv-1] = not_visited[0];
		Lnew_tour += weights[new_tour[Nv-2]][new_tour[Nv-1]];
		Lnew_tour += weights[new_tour[Nv-1]][new_tour[0]];
		if(i == first_city || Lnew_tour < Lbest_tour){
			best_tour = new_tour;
			Lbest_tour = Lnew_tour;
		}
	}
		
	
	//change representation of tour
	decision_vector individual(Nv*(Nv-1), 0);
	for (size_t j = 0; j < Nv; j++) {
		individual[(best_tour[j])*(Nv-1)+best_tour[(j+1>Nv-1? 0:j+1)] - (best_tour[(j+1>Nv-1? 0:j+1)]>best_tour[j]? 1:0)] = 1;
	}
	population::size_type best_idx = pop.get_best_idx();
	pop.set_x(best_idx,individual);


    } // end of evolve
    

	
    /// Algorithm name
    std::string nn::get_name() const
    {
        return "Nearest neighbour algorithm";
    }

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::nn)
