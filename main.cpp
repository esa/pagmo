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
#include <iomanip>
#include "src/pagmo.h"

#include "src/algorithm/cstrs_co_evolution.h"

using namespace pagmo;

// Example in C++ of the use of PaGMO 1.1.4

int main()
{
	pagmo::problem::dtlz7 prob(20);
	pagmo::algorithm::spea2 alg(100);

	std::cout << alg << std::endl;
	std::cout << prob << std::endl;

	pagmo::island isl = island(alg, prob, 100);

	for (size_t i = 0; i<1; ++i){
		isl.evolve(100);
		isl.join();
		std::cout << "Distance from Pareto Front (p-distance): " << prob.p_distance(isl.get_population()) << std::endl;
		//std::cout << "Original fitness: " << isl.get_population() << std::endl;
		//std::cout << "Decomposed fitness: " << decomposed_problem.objfun(isl.get_population().champion().x) << std::endl;
		std::vector<population::size_type> front0_index = isl.get_population().compute_pareto_fronts()[0];
		for(int i=0; i<front0_index.size(); i++){
			std::cout<<isl.get_population().get_individual(front0_index[i]).cur_f<<std::endl;
		}
	}
	std::cout << "Finished" << std::endl;
	return 0;
}
