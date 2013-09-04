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
	//pagmo::problem::welded_beam prob_constrained;
	//pagmo::problem::tens_comp_string prob_constrained;
	//pagmo::problem::pressure_vessel prob_constrained;

	pagmo::algorithm::de algo(1, 0.8, 0.9, 2, 1e-15, 1e-15);
//	pagmo::algorithm::cmaes algo(1);
//	pagmo::algorithm::cmaes algo_2(70);

#ifndef PAGMO_ENABLE_GSL
	pagmo_throw(assertion_error,"The GSL library must be installed to use the repair main example");
	return 1;
#endif

#ifdef PAGMO_ENABLE_GSL
	pagmo::algorithm::gsl_nm2 repair_algo(100, 1e-5, 0.02);
	pagmo::algorithm::cstrs_core algo_constrained(algo, repair_algo, 1000);
	algo_constrained.reset_rngs(100);

	std::cout << algo_constrained;

	for (size_t i=0; i<1; ++i) {
		pagmo::population pop(prob_constrained,90);
		pagmo::population pop_copy = pop;

//		std::cout << "---------------BEFORE------------" << std::endl;
//		std::cout << pop.get_individual(0);
//		std::cout << "---------------------------------" << std::endl;
//		pop.repair(0);
//		std::cout << "---------------AFTER-------------" << std::endl;
//		std::cout << pop.get_individual(0);
//		std::cout << "---------------------------------" << std::endl;

		algo_constrained.evolve(pop);
		std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<std::endl;
		std::cout<<"CHAMPION1"<<std::endl;
		std::cout << pop.champion();
		std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<std::endl;
//		algo_constrained.evolve(pop_copy);
//		std::cout<<"CHAMPION2"<<std::endl;
//		std::cout << pop_copy.champion();
	}

	std::cout << algo_constrained << std::endl;
	std::cout << prob_constrained << std::endl;

	std::cin.get();
//	pagmo::island isl = island(algo_constrained, prob_constrained, 70);

//	for (size_t i=0; i<20; ++i){
//		isl.evolve(1);
//		std::cout << isl.get_population().champion();
//	}

#endif
	return 0;
}
