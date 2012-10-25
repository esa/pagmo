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
#include <fstream>
#include <iostream>
#include <algorithm>

// include headers that implement a archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "../src/pagmo.h"

using namespace pagmo;

int main()
{
	unsigned int dimension = 10;

	// create a container of pagmo::problems
	std::vector<problem::base_ptr> probs;
	std::vector<problem::base_ptr> probs_new;
	
	// fill it up with problems
	probs.push_back(problem::ackley(dimension).clone());
	probs_new.push_back(problem::ackley().clone());
	probs.push_back(problem::rosenbrock(dimension).clone());
	probs_new.push_back(problem::rosenbrock().clone());

	for (size_t i=0; i< probs.size(); ++i) {
		// create and open a character archive for output
		std::ofstream ofs("test.ar");
		// save data to archive
		{
		boost::archive::text_oarchive oa(ofs);
		// write class instance to archive
		oa << (*(probs[i]));
		// archive and stream closed when destructors are called
		}
	
		{
		// create and open an archive for input
		std::ifstream ifs("test.ar");
		boost::archive::text_iarchive ia(ifs);
		// read class state from archive
		ia >> *(probs_new[i]);
		// archive and stream closed when destructors are called
		}
		
		{
		decision_vector x(probs[i]->get_dimension(),0);
		fitness_vector f1(probs[i]->get_f_dimension(),0), f2(probs[i]->get_f_dimension(),0);
		population pop(*(probs[i]),1);
		x = pop.champion().x;
		probs[i]->objfun(f1,x);
		probs_new[i]->objfun(f2,x);
		std::cout << f1 << " " << f2 << std::endl;
		if (std::equal(f1.begin(),f1.end(),f2.begin())) {
			continue;
		} else { return 1;}
		}
	
	}
	return 0;
}
