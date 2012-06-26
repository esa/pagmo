/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
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
#include "src/pagmo.h"

using namespace pagmo;

int main()
{
// This instantiates a differential evolution algorithm that will run for 500 generations. Refer to the documentation to

// see what othert parameters do

pagmo::algorithm::cmaes algo(100);

algo.set_screen_output(true);

//This instantiate a 50 dimensional Rosenbrock problem
pagmo::problem::zdt1 prob;


std::cout << alg_mde << std::endl;

alg_mde.set_screen_output(true);
alg_jde.set_screen_output(false);
alg_de.set_screen_output(false);

//This instantiate a 10 dimensional Schwefel problem
pagmo::problem::rastrigin prob(10);
// prob.set_bounds(-100, 100);

// pagmo::population pop(prob, 100);

//This instantiate an island containing a population of 10 individuals initialized at random and having their fitness evaluated
//with respect to the Schwefel problem. The island will evolve its population using the instantiated algorithm
pagmo::island isla = island(alg_mde, prob, 100);
pagmo::island islb = island(alg_jde, prob, 100);
pagmo::island islc = island(alg_de, prob, 100);

//This prints on screen the instantiated Schwefel problem
std::cout << prob << std::endl;

// isla.evolve(1);
// isla.join();
// pagmo::population popa = isla.get_population();

// islb.evolve(1);
// islb.join();
// pagmo::population popb = islb.get_population();

isla.evolve(1);
isla.join();
islb.evolve(1);
islb.join();
islc.evolve(1);
islc.join();

// pagmo::population pop = islb.get_population();
// alg_jde.evolve(pop);

// islc.evolve(1);
// islc.join();
// pagmo::population popc = islc.get_population();

//Evolution is here started on the single island instantiated

for (int i=0; i< 30; ++i){
	algo.evolve(pop);
	std::cout << pop.champion().f[0] << " " << std::endl;
}


return 0;
}
