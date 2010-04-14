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

#include <climits>
#include <iostream>
#include <vector>
#include <list>

#include "src/algorithms.h"
#include "src/archipelago.h"
#include "src/island.h"
#include "src/problems.h"
#include "src/topologies.h"
#include "src/topologies.h"
#include "src/problem/base.h"

using namespace pagmo;

int main()
{

double x[22] = {-779.45156677047862, 3.2540570452459043, 0.52689731549926921, 0.38285203950839292, 167.55049483170021, 424.23381445972808, 53.310310073616151, 589.77169016660821, 2199.9939278768511, 0.7755487837508489, 0.53666875224363364, 0.010546794232542851, 0.014653697948643706, 0.6594344155776396, 1.3505009019745324, 1.0500032115842122, 1.3066951686999049, 69.81346649423962, -1.5907330914593276, -1.9595633249627791, -1.5547363420059905, -1.5134315210520082};
	algorithm::nlopt_sbplx algo2(1000, 1e-10);
	algorithm::mbh algo(algo2,500,1e-6/5);
	algo.screen_output(true);
	problem::cassini_2 prob;

	island isl = island(prob,algo,20);
	std::cout << prob << std::endl;
	std::cout << algo << std::endl;
	std::vector<double> x0(x,x+22);
	isl.set_x(0,x0);

	for (int i=0; i< 20; ++i){
		//isl.set_algorithm(algorithm::sa_corana(10000,(isl.get_population().champion().f[0]-9),0.1));
		isl.evolve(); isl.join();
		std::cout << i << " - " << isl.get_population().champion().f << " " << problem::objfun_calls() << std::endl;
		std::cout << isl.get_population().champion().x <<std::endl;
	}

	return 0;
}
