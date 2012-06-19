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

pagmo::algorithm::mde_pbx mde(1000);
pagmo::algorithm::jde jde(1000);

mde.set_screen_output(true);
jde.set_screen_output(true);

//This instantiate a 10 dimensional Schwefel problem
pagmo::problem::rosenbrock prob(10);

pagmo::population pop(prob, 100);

//This instantiate an island containing a population of 10 individuals initialized at random and having their fitness evaluated
//with respect to the Schwefel problem. The island will evolve its population using the instantiated algorithm
pagmo::island isla = island(mde, pop);
pagmo::island islb = island(jde, pop);

//This prints on screen the instantiated Schwefel problem
std::cout << prob << std::endl;

islb.evolve();
pagmo::population popb = islb.get_population();

isla.evolve();
pagmo::population popa = isla.get_population();


// double ds[] = {1.0, 2.0, 3.0};
// std::vector<double> v(ds, ds + sizeof(ds) / sizeof(double));
// 
// std::cout << "Powermean: " << mde.powermean(v) << std::endl;

// pagmo::population popa = isla.get_population();
// std::cout << popa.champion().f[0] << " evolving..." << std::endl;
// 
// mde.evolve(popa);


std::cout << "Result jDE: " << popb.champion().f[0] << " " << std::endl;
std::cout << "Result MDE_pBX: " << popa.champion().f[0] << " " << std::endl;

return 0;
}
