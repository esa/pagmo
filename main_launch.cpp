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

#include <iostream>
#include <fstream>

#include "src/pagmo.h"

using namespace pagmo;

int main()
{
	mpi_environment env;

	int n_segments = 10;

	algorithm::snopt algo(1000);

	const double pert_epoch = 1000;
	const double pert_nondim = 1E-1;
	const double pert_mass = 200;
	const double pert_vinf = 1000;
	std::vector<double> perturb(n_segments * 3 + 6,pert_nondim);
	perturb[0] = pert_epoch;
	perturb[1] = pert_epoch;
	perturb[2] = pert_mass;
	perturb[3] = pert_vinf;
	perturb[4] = pert_vinf;
	perturb[5] = pert_vinf;

	algorithm::mbh algo2(algo,20,perturb);
	algorithm::ms algo3(algo2,5);
	//algo2.screen_output(true);

	// Outer loop.
	int cur_id = 1;
	for (int i = 1; i <= 25; ++i) {
		archipelago a;
		for (int j = 0; j < 283; ++j) {
			a.push_back(mpi_island(problem::gtoc5_launch(n_segments,cur_id + j,problem::gtoc5_launch::TIME),algo3,5));
		}
		a.evolve(1);
		a.join();
		std::ofstream ofs((std::string("time_") + boost::lexical_cast<std::string>(cur_id)).c_str());
		boost::archive::text_oarchive oa(ofs);
		oa << static_cast<const archipelago &>(a);
		cur_id += 283;
	}

	cur_id = 1;
	for (int i = 1; i <= 25; ++i) {
		archipelago a;
		for (int j = 0; j < 283; ++j) {
			a.push_back(mpi_island(problem::gtoc5_launch(n_segments,cur_id + j,problem::gtoc5_launch::MASS),algo3,5));
		}
		a.evolve(1);
		a.join();
		std::ofstream ofs((std::string("mass_") + boost::lexical_cast<std::string>(cur_id)).c_str());
		boost::archive::text_oarchive oa(ofs);
		oa << static_cast<const archipelago &>(a);
		cur_id += 283;
	}
}
