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
	// create and open a character archive for output
	std::ofstream ofs("test.ar");

	unsigned int dimension = 10;
	// create class instance
	const pagmo::problem::ackley prob(dimension);

	// save data to archive
	{
		boost::archive::text_oarchive oa(ofs);
		// write class instance to archive
		oa << prob;
		// archive and stream closed when destructors are called
	}
	
	
	pagmo::problem::ackley prob_new;
	{
		// create and open an archive for input
		std::ifstream ifs("test.ar");
		boost::archive::text_iarchive ia(ifs);
		// read class state from archive
		ia >> prob_new;
		// archive and stream closed when destructors are called
	}
	decision_vector x(prob.get_dimension(),0);
	fitness_vector f1(prob.get_f_dimension(),0), f2(prob.get_f_dimension(),0);
	population(prob,1);
	prob.objfun(f1,x);
	prob_new.objfun(f2,x);
	return (!std::equal(f1.begin(),f1.end(),f2.begin()));
}
