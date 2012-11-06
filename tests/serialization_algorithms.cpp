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
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

// include headers that implement a archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "../src/pagmo.h"

using namespace pagmo;
int main()
{
	unsigned int pop_size = 20;
	unsigned int gen = 10;

	// create two containers of pagmo::algorithms
	std::vector<algorithm::base_ptr> algos;
	std::vector<algorithm::base_ptr> algos_new;
	
	// fill it up with algorithm
	algos.push_back(algorithm::jde(gen,7,2).clone());
	algos_new.push_back(algorithm::jde().clone());
	algos.push_back(algorithm::bee_colony(gen,19).clone());
	algos_new.push_back(algorithm::bee_colony().clone());
	algos.push_back(algorithm::cmaes(gen,0.5, 0.5, 0.5, 0.5, 0.7, 1e-5, 1e-5, false).clone());
	algos_new.push_back(algorithm::cmaes().clone());
	
	// pick a box-constrained, single objective, continuous problem
	problem::ackley prob(10);
	
	// make a population out of it
	population pop_original(prob,pop_size);

	//serialize algos and deserialize into algos_new checking they are then identical
	for (size_t i=0; i< algos.size(); ++i) {
		{
		// create and open a character archive for output
		std::ofstream ofs("test.ar");
		// save data to archive
		boost::archive::text_oarchive oa(ofs);
		// write class instance to archive
		oa & algos[i];
		// archive and stream closed when destructors are called
		}
	
		{
		// create and open an archive for input
		std::ifstream ifs("test.ar");
		boost::archive::text_iarchive ia(ifs);
		// read class state from archive
		ia & algos_new[i];
		// archive and stream closed when destructors are called
		}
		
		{
		//copy the original population
		population pop1(pop_original), pop2(pop_original);
		algos[i]->evolve(pop1);
		algos_new[i]->evolve(pop2);
		std::cout << std::endl << std::setw(40) << algos[i]->get_name();
		decision_vector x1(pop1.champion().x), x2(pop2.champion().x);
		if (std::equal(x1.begin(),x1.end(),x2.begin())) {
			std::cout << ": pass";
		} else {
			std::cout << ": Champion is different:";
			return 1;
		}
		}
	}
	std::cout << std::endl;
	return 0;
}
