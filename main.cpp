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

using namespace pagmo;

// Example in C++ of the use of PaGMO 1.1.5

int main() {
	/*std::vector<kep_toolbox::planet_ptr> seq;
	seq.push_back(kep_toolbox::planet_js("callisto").clone());
	seq.push_back(kep_toolbox::planet_js("ganymede").clone());
	seq.push_back(kep_toolbox::planet_js("ganymede").clone());
	seq.push_back(kep_toolbox::planet_js("ganymede").clone());
	
	std::vector<std::vector<double> > tofs;
	std::vector<double> dumb(2);
			dumb[0] = 180;dumb[1] = 200;
			tofs.push_back(dumb);
			dumb[0] = 0.1;dumb[1] = 5;
			tofs.push_back(dumb);
			dumb[0] = 10;dumb[1] = 150;
			tofs.push_back(dumb);
			dumb[0] = 10;dumb[1] = 40;
			tofs.push_back(dumb);
	problem::mga_incipit_cstrs prob(seq, kep_toolbox::epoch(10460.0), kep_toolbox::epoch(104803.0),tofs);
	
	algorithm::jde algo(50), algo2(1);
	algorithm::cstrs_co_evolution algo_coevo(algo,algo2,10,10,algorithm::cstrs_co_evolution::SPLIT_CONSTRAINTS,0,1);

	for(int j=0;j<100;j++){
		population pop(prob,20);
		for(int i=0;i<100;i++){
			algo_coevo.evolve(pop);
		}
		std::cout<<pop.champion().f<<std::endl;
	}*/

	problem::zdt prob(1,30);
	population pop(prob,100);
	algorithm::moead algo;
	algo.evolve(pop);
	std::cout << pop << std::endl;
	return 0;
}
