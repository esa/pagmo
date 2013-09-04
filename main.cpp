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

using namespace pagmo;

// Example in C++ of the use of PaGMO 1.1.4

int main()
{
	pagmo::problem::zdt1 orig_prob(10);

//	std::vector<double> weights(2,0.5);
//	pagmo::problem::decompose decomposed_problem(orig_prob, weights);

	// pagmo::algorithm::jde alg(50);

	pagmo::algorithm::vega alg(50);

	std::cout << alg << std::endl;
	std::cout << orig_prob << std::endl;
	//std::cout << decomposed_problem << std::endl;

	pagmo::island isl = island(alg, orig_prob, 50);
	//pagmo::population original_problem_pop = population(orig_prob, 50);

	std::cout << "Initial distance from Pareto Front (p-distance): " << orig_prob.p_distance(isl.get_population()) << std::endl;
	for (size_t i = 0; i< 10; ++i){
		isl.evolve();
		//original_problem_pop.set_x(0, isl.get_population().champion().x);
		std::cout << "Distance from Pareto Front (p-distance): " << orig_prob.p_distance(isl.get_population()) << std::endl;
		//std::cout << "Original fitness: " << orig_prob.objfun(isl.get_population().champion().x) << std::endl;
		//std::cout << "Decomposed fitness: " << decomposed_problem.objfun(isl.get_population().champion().x) << std::endl;

	}
	return 0;
}
