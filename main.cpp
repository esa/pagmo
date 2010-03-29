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
 *   the Free Software Foundation; either version 3 of the License, or       *
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

#include <algorithm>
#include <climits>
#include <iostream>
#include <string>
#include <vector>
#include <list>

#include "src/algorithms.h"
#include "src/archipelago.h"
#include "src/island.h"
#include "src/problems.h"
#include "src/topologies.h"
#include "src/keplerian_toolbox/planet.h"

using namespace pagmo;

int main()
{
	archipelago a(problem::cassini_1(),algorithm::gsl_bfgs2(40),1,20,topology::ring());
	a.push_back(island(problem::cassini_1(),algorithm::gsl_bfgs2(40),1));
	for (int i = 0; i < 10; ++i) {
		a.evolve();
		std::cout << a.get_island(0).get_population().champion().f[0] << '\n';
		std::cout << a.get_island(1).get_population().champion().f[0] << '\n';
// 		std::cout << a.get_island(2).get_population().champion().f[0] << '\n';
// 		std::cout << a.get_island(3).get_population().champion().f[0] << '\n';
	}
#if 0
	algorithm::ihs ihs_instance(10000,.999,.001);
	const std::string str = "hello world, cruel world oh noes hello world, cruel world oh noes hello world, cruel world oh noes hello world, cruel world oh noes hello world, cruel world oh noes";
	problem::string_match_mo prob(str);
	island isl = island(prob,ihs_instance,120 * 5);
	for (int i = 0; i < 100; ++i) {
		isl.evolve(); isl.join();
		std::string tmp(str);
		population pop = isl.get_population();
		std::copy(pop.get_individual(pop.get_best_idx()).cur_x.begin(),pop.get_individual(pop.get_best_idx()).cur_x.end(),tmp.begin());
		std::cout << tmp << '\n';
		std::cout << pop.get_individual(pop.get_best_idx()).cur_f << '\n';
	}
	std::cout << isl;
#endif
}
