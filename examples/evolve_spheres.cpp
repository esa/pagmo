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

#include <iostream>
#include<numeric>
#include<iomanip>
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


int main()
{
// EXPERIMENT SET-UP //
		const int n_isl = 8;
		const int pop_size = 512;
		const int n_eval = 5;
		const int n_gen = 1000;
		std::vector<double> sides(3,0.0); 
		sides[0] = 0.6; sides[1] = 0.7; sides[2] = 0.8;
// END OF EXPERIMENT SET-UP //

// Buffer
std::vector<double> buff;
// We instantiate a PSO algorithm capable of coping with stochastic prolems
algorithm::pso_generational algo(1,0.7298,2.05,2.05,0.05);

// This instantiates the spheres problem
std::cout << "Initializing ....";

archipelago archi = archipelago(topology::fully_connected());

for (int j=0;j<n_isl; ++j) {
	problem::spheres prob(n_eval,10,1e-9,rand(), true, 50.0, sides);
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
