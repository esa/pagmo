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
	double lb[26] = {1.900000000000000000e+03,3.000000000000000000e+00,
0.000000000000000000e+00,0.000000000000000000e+00,
3.500000000000000000e+02,2.200000000000000000e+02,
2.040000000000000000e+02,2.600000000000000000e+02,
3.460000000000000000e+02,5.190000000000000000e+02,
1.000000000000000021e-02,1.000000000000000021e-02,
1.000000000000000021e-02,1.000000000000000021e-02,
1.000000000000000021e-02,1.000000000000000021e-02,
1.100000000000000089e+00,1.100000000000000089e+00,
1.050000000000000044e+00,1.050000000000000044e+00,
1.050000000000000044e+00,-3.141592653589793116e+00,
-3.141592653589793116e+00,-3.141592653589793116e+00,
-3.141592653589793116e+00,-3.141592653589793116e+00};
	double ub[26] = {2.200000000000000000e+03,4.049999999999999822e+00,
1.000000000000000000e+00,1.000000000000000000e+00,
4.800000000000000000e+02,2.510000000000000000e+02,
3.200000000000000000e+02,2.670000000000000000e+02,
3.590000000000000000e+02,5.360000000000000000e+02,
9.899999999999999911e-01,9.899999999999999911e-01,
9.899999999999999911e-01,9.899999999999999911e-01,
9.899999999999999911e-01,9.899999999999999911e-01,
6.000000000000000000e+00,6.000000000000000000e+00,
6.000000000000000000e+00,6.000000000000000000e+00,
6.000000000000000000e+00,3.141592653589793116e+00,
3.141592653589793116e+00,3.141592653589793116e+00,
3.141592653589793116e+00,3.141592653589793116e+00};
	double x[26] = {2.052037404713874821e+03,3.954442941493409691e+00,
4.511253386461788195e-01,6.768791547431635136e-01,
4.325075205985705225e+02,2.274230910720812346e+02,
2.216682548627323399e+02,2.663337584131617746e+02,
3.583157984809776053e+02,5.343501740101467021e+02,
7.748270026863098847e-01,2.969998105989154480e-01,
6.750794167470981488e-01,6.760615415138029327e-01,
8.248926524434435636e-01,9.011521085766600603e-01,
1.357545130899178165e+00,1.139696652938946730e+00,
1.304198680067162774e+00,1.061357582516698761e+00,
1.050000000000000044e+00,2.865062314760247641e+00,
1.557300412042473159e+00,1.741399754610104100e+00,
1.578846812880740913e+00,1.569549219255385708e+00};

	algorithm::snopt algo2(10000,1e-14,1e-14);
	//algo2.screen_output(false);
	algorithm::mbh algo(algo2,200,0.01);
	algo.screen_output(true);
	problem::messenger_full prob;

	std::vector<double> lb0(lb,lb+26);
	std::vector<double> ub0(ub,ub+26);
	prob.set_bounds(lb0,ub0);
	island isl = island(prob,algo,1);
	std::cout << prob << std::endl;
	std::cout << algo << std::endl;

	std::vector<double> x0(x,x+26);
	isl.set_x(0,x0);

	for (int i=0; i< 20; ++i){
		//isl.set_algorithm(algorithm::sa_corana(10000,(isl.get_population().champion().f[0]-9),0.1));
		isl.evolve(); isl.join();
		std::cout << isl.get_population().champion().f << " " << problem::objfun_calls() << std::endl;
	}
	return 0;
}
