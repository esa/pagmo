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

#include "src_new/algorithm/ihs.h"
#include "src_new/algorithm/monte_carlo.h"
#include "src_new/algorithm/null.h"
#include "src_new/archipelago.h"
#include "src_new/island.h"
#include "src_new/problem/golomb_ruler.h"
#include "src_new/problem/himmelblau.h"
#include "src_new/problem/knapsack.h"
#include "src_new/problem/paraboloid.h"
#include "src_new/problem/rastrigin.h"
#include "src_new/topology/ring.h"
#include "src_new/topology/one_way_ring.h"
#include "src_new/topology/unconnected.h"
#include "src_new/topology/fully_connected.h"
#include "src_new/topology/erdos_renyi.h"
#include "src_new/topology/barabasi_albert.h"

using namespace pagmo;

int main()
{
	topology::barabasi_albert t;
	for (int i = 0; i < 100; ++i) {
		t.push_back(i);
	}
	//t.push_back(4);
	//t.push_back(5);
	//t.push_back(6);
	//std::cout << t << '\n';
	//return 0;

	double lb1[] = {-1,-1};
	double ub1[] = {-1,-1};
	std::cout << problem::paraboloid(lb1,ub1) << '\n';

	std::vector<double> lb2(4,0);
	std::vector<double> ub2(4,10);
	problem::paraboloid p2(lb2,ub2);
	p2.set_bounds(lb2,ub2);
	p2.set_bounds(lb2.begin(),lb2.end(),ub2.begin(),ub2.end());
	p2.set_bounds(lb2.begin(),lb2.end(),ub2.begin(),ub2.end());
	double lb2a[] = {-1,-1,0,1};
	double ub2a[] = {-1,-1,1,1};
	p2.set_bounds(lb2a,ub2a);
	p2.set_lb(0,-4.);
	std::list<double> ub2b(ub2.begin(),ub2.end());
	p2.set_bounds(lb2a,lb2a + 4,ub2b.begin(),ub2b.end());
	p2.set_bounds(lb2,ub2);
	fitness_vector f(1);
	p2.objfun(f,ub2);
	ub2[0] -= 1E-3;
	lb2[0] += 1E-3;
	p2.objfun(f,ub2);
	ub2[0] -= 1E-3;
	lb2[0] += 1E-3;
	p2.objfun(f,lb2);
	ub2[0] -= 1E-3;
	lb2[0] += 1E-3;
	p2.objfun(f,ub2);
	ub2[0] -= 1E-3;
	lb2[0] += 1E-3;
	p2.objfun(f,lb2);
	ub2[0] -= 1E-3;
	lb2[0] += 1E-3;
	p2.objfun(f,ub2);
	ub2[0] -= 1E-3;
	lb2[0] += 1E-3;
	p2.objfun(f,lb2);
	ub2[0] -= 1E-3;
	lb2[0] += 1E-3;
	p2.objfun(f,ub2);
	ub2[0] -= 1E-3;
	lb2[0] += 1E-3;
	p2.objfun(f,lb2);
	ub2[0] -= 1E-3;
	lb2[0] += 1E-3;
	p2.objfun(f,ub2);
	ub2[0] -= 1E-3;
	lb2[0] += 1E-3;
	p2.objfun(f,lb2);
	ub2[0] -= 1E-3;
	lb2[0] += 1E-3;
	p2.objfun(f,ub2);
	ub2[0] -= 1E-3;
	lb2[0] += 1E-3;
	p2.objfun(f,lb2);
	p2.set_lb(0,-1.234);
	std::cout << f << '\n';

	double values[] = {1,2,3,4,5,6,7,8,9,11,15};
	double weights[] = {4,7,9,1,2,3,5,8,6,12,17};

//  	population pop(problem::knapsack(values,weights,30),10);
// 	population pop(problem::paraboloid(lb2,ub2),10);
//	population pop(problem::golomb_ruler(13,169),10);


	archipelago archi = archipelago(topology::barabasi_albert());
	archi.push_back(island(problem::knapsack(values,weights,30),algorithm::ihs(1000),10));
	archi.push_back(island(problem::knapsack(values,weights,30),algorithm::ihs(1000),10));

	std::cout << archi;

	//algorithm::ihs algo(1000000);

	//isl.evolve_t(1000 * 10);
	
	//std::cout << isl;

	//std::cout << pop << '\n';

	//std::cout << algorithm::ihs(5) << '\n';

	//algorithm::null().evolve(isl);

	//std::cout << isl << '\n';


	//std::cout << problem::paraboloid(lb2,ub2).objfun(lb2) << '\n';
}
