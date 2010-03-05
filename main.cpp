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

#include "src/algorithm/de.h"
#include "src/algorithm/pso.h"
#include "src/algorithm/sa_corana.h"
#include "src/algorithm/ihs.h"
#include "src/algorithm/monte_carlo.h"
#include "src/algorithm/null.h"
#include "src/archipelago.h"
#include "src/island.h"
#include "src/problem/golomb_ruler.h"
#include "src/problem/himmelblau.h"
#include "src/problem/knapsack.h"
#include "src/problem/paraboloid.h"
#include "src/problem/rastrigin.h"
#include "src/problem/rosenbrock.h"
#include "src/topology/ring.h"
#include "src/topology/one_way_ring.h"
#include "src/topology/unconnected.h"
#include "src/topology/fully_connected.h"
#include "src/topology/custom.h"
#include "src/topology/erdos_renyi.h"
#include "src/topology/barabasi_albert.h"
#include "src/topology/watts_strogatz.h"

using namespace pagmo;

int main()
{
	island isl = island(problem::rosenbrock(25),algorithm::sa_corana(10000,2,0.001),20);
	std::cout << isl.get_population().champion().f << std::endl;
	for (int i=0; i< 50; ++i){
		isl.evolve();
		std::cout << isl.get_population().champion().f << std::endl;
	}
}
