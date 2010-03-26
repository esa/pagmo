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
#include "src/keplerian_toolbox/planet.h"

using namespace pagmo;

int main()
{
	algorithm::snopt ipopt_instance(10000,1e-4,1e-4);
	ipopt_instance.screen_output(true);
// 	snopt_instance.file_output(false);
	problem::earth_planet prob(1,kep_toolbox::planet::MARS);

	island isl = island(prob,ipopt_instance,1);
	std::cout << prob << std::endl;
	std::cout << isl << std::endl;

	for (int i=0; i< 1; ++i){
		isl.evolve(); isl.join();
		std::cout << isl.get_population().champion().f[0] << " " << problem::objfun_calls() << std::endl;
	}
	//std::cout << isl << std::endl;
}
