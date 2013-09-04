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
	pagmo::problem::cec2006 prob_constrained(5);
	
	//pagmo::algorithm::monte_carlo algo(1); //only one generation for the algo!
	//pagmo::algorithm::sga algo(1); //only one generation for the algo!
	pagmo::algorithm::sga_gray algo(1,0.9,0.003, 1,
									algorithm::sga_gray::mutation::UNIFORM,
									algorithm::sga_gray::selection::ROULETTE,
									algorithm::sga_gray::crossover::SINGLE_POINT); //only one generation for the algo!
	//pagmo::algorithm::cmaes algo(1); //only one generation for the algo!
	//pagmo::algorithm::de algo(1); //only one generation for the algo!
	//pagmo::algorithm::pso algo(1); //only one generation for the algo!
    pagmo::algorithm::cstrs_self_adaptive algo_constrained(algo, 5000);

	std::cout << algo_constrained;

	for (size_t i=0; i<20; ++i) {
		pagmo::population pop(prob_constrained,70);
		algo_constrained.evolve(pop);
		std::cout << pop.champion();
	}

	std::cout << algo_constrained << std::endl;
	std::cout << prob_constrained << std::endl;

//	// test sga_gray
//	pagmo::problem::branin prob;
//	pagmo::algorithm::sga_gray algo(1000,0.9,0.003, 1,
//									algorithm::sga_gray::mutation::UNIFORM,
//									algorithm::sga_gray::selection::ROULETTE,
//									algorithm::sga_gray::crossover::SINGLE_POINT);

//	pagmo::population pop(prob,70);
//	algo.evolve(pop);
//	std::cout << pop.champion();




//	pagmo::island isl = island(algo_constrained, prob_constrained, 70);

//	for (size_t i=0; i<20; ++i){
//		isl.evolve(1);
//		std::cout << isl.get_population().champion();
//	}
	return 0;
}
