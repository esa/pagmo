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

#include <algorithm>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/utility/result_of.hpp>
#include <cmath>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <numeric>
#include <string>
#include <vector>

#include "../src/pagmo.h"
#include "test_functions.h"

int main(int argc, char* argv[])
{
	boost::mpi::environment env(argc, argv);
	boost::mpi::communicator world;
	pagmo::algorithm::de de = pagmo::algorithm::de(500,.8,.8,2);
	const std::vector<pagmo::problem::base_ptr> probs(pagmo::get_test_problems());
	std::cout << "Testing algorithm: " << de << '\n';
	
	for (std::vector<pagmo::problem::base_ptr>::const_iterator it = probs.begin(); it != probs.end(); ++it)
	{
		std::cout << "\tTesting problem: " << (**it) << '\n';
		std::vector<double> champs;
		pagmo::archipelago a = pagmo::archipelago(pagmo::topology::rim());
			
		a.push_back(pagmo::mpi_island(**it,de,20,1));
		a.push_back(pagmo::mpi_island(**it,de,20,2));		
		a.push_back(pagmo::mpi_island(**it,de,20,3));
				
		a.evolve(10);
		a.join();
	
	
		if (world.rank() == 0)
		{
			std::cout << a.get_island(0)->get_population().champion().f << " " << a.get_island(0)->get_population().champion().x << std::endl;
		}

	}
	return 0;
}
