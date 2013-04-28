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
#include "src/pagmo.h"

using namespace pagmo;

// Example in C++ of the use of PaGMO 1.1.4

int main()
{
	pagmo::algorithm::mde_pbx alg(3000, 0.15, 1.5, 1e-30, 1e-30);
	//pagmo::algorithm::de alg(1500,0.8,0.5,8);
	alg.set_screen_output(true);
	std::cout << alg << std::endl;
	pagmo::problem::rosenbrock prob(30);
	prob.set_bounds(-100,100);
	std::cout << prob << std::endl;
	pagmo::island isl = island(alg, prob, 100);
	for (size_t i = 0; i< 1; ++i){
		isl.evolve(1);
		std::cout << isl.get_population().champion().f << std::endl;
		//std::cout << "Distance from Pareto Front (p-distance): " << prob.p_distance(isl.get_population()) << std::endl;
	}
	isl.join(); //This is important to make sure the program is not terminated before the threads
				//are destroyed, thus invoking illegal calls to pure virtual methods .....
	return 0;
}
