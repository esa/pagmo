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

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <iostream>
#include <fstream>

#include "src/mpi_environment.h"
#include "src/mpi_island.h"
#include "src/archipelago.h"
#include "src/topology/ring.h"
#include "src/problem/schwefel.h"
#include "src/algorithm/de.h"

using namespace pagmo;

int main()
{
// 	std::ofstream ofs("schwefel.txt");
// 	boost::archive::binary_oarchive oa(ofs);
// 	population pop(problem::schwefel(10),1);
// 	oa << pop;
// 	ofs.close();
// 	std::ifstream ifs("schwefel.txt");
// 	boost::archive::binary_iarchive ia(ifs);
// 	ia >> pop;
// 	std::cout << pop << '\n';
// 	return 0;

	mpi_environment env;
	archipelago a = archipelago(topology::ring());
	for (int i = 1; i < boost::mpi::communicator().size(); ++i) {
		a.push_back(mpi_island(problem::schwefel(300),algorithm::de(500),10));
	}
	a.evolve(2);
	a.join();
std::cout << "finished evolving\n";
}
