/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
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

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <string>

#include <keplerian_toolbox/keplerian_toolbox.h>
#ifdef PAGMO_ENABLE_NLOPT
	#include"../src/algorithm/nlopt_sbplx.h"
#else
	#include"../src/algorithm/cs.h"
#endif
#include"../src/algorithm/sa_corana.h"
#include"../src/algorithm/de.h"
#include"../src/algorithm/cs.h"
#include"../src/archipelago.h"
#include"../src/island.h"
#include"../src/problem/sample_return.h"
#include"../src/topologies.h"
#include"../src/migration/random_r_policy.h"

using namespace pagmo;
using namespace kep_toolbox;

/*
WARNING: This example requires the file MPCORB.DAT to be placed in the same directory as
the executable. Such a file is freely downloadable from the internet and contains
all the information on the minor plantes needed to define their physical properties.

PROBLEM DESCRITPION: A Human Mission to Asteroids. In this example we solve a large number of different instances of the
pagmo::problem::sample_return. In particular we solve one instance for each minor planet in the database.
A filter is put on the semi-major axis as the entire database contains 500.000 objects that could make the CPU busy
for weeks with not much use. During the optimization we output on screen the current asteroid and
the objective function found in a particular multistart. All results are stored in out.pagmo file
The best asteroids found are suitable for a human mission in the 2020-2050 farme.

OPTIMIZATION STRATEGY: we use an archipelago with seven islands having populations evolved with differential
evolution, simulated annealing and either compass search or nlopt_sbplx, according to the option specified in ccmake.
The topology is arim with the local optimizer at the center.

CPU TIME: the example opens 7 threads and needs to perform a large number of optimization tasks, so it does
occupy the machine CPUs for a long time (days) according to the pruning level adopted
*/

int main()
{
	//Number of multistart to be done per asteroid
	int n_multistart = 1;

	//Open the output file
	std::ofstream myfile;
	myfile.open((std::string("pagmo_") +
		boost::lexical_cast<std::string>(rng_generator::get<rng_uint32>()()) + ".out").c_str());

	algorithm::sa_corana algo1(10000,1,0.01);
	algorithm::de algo2(500,0.8,0.8,3);
#ifdef PAGMO_ENABLE_NLOPT
	algorithm::nlopt_sbplx algo3(500,1e-4);
#else
	algorithm::cs algo3(500,0.0001,0.1);
#endif

	double Tmax = 600;
	
	//Opening the MPCORB.DAT file
	std::ifstream mpcorbfile("MPCORB.DAT");
	std::string line;
	if (!mpcorbfile) {
		throw_value_error("Could not find MPCORB.DAT is it in the current directory?");
	}
	
	//Skipping the first lines
	do {
		std::getline(mpcorbfile,line);
	} while (!boost::algorithm::find_first(line,"-----------------"));
	
	while(!mpcorbfile.eof()) {
		try {
			std::getline(mpcorbfile,line);
			planet::mpcorb target(line);

			//Pruning based on the asteroid elements
			array6D elem = target.get_elements();
			if (elem[0] / ASTRO_AU < 1.8) {

				//build the problem
				problem::sample_return prob(target,Tmax);

				for (int k=0;k<n_multistart;++k){
					std::cout << "\tTarget is: " << target.get_name() << ", Trial: " << k << std::endl;
					archipelago a = pagmo::archipelago(topology::rim());
					a.push_back(pagmo::island(algo3,prob,1,migration::best_s_policy(),migration::random_r_policy()));
					a.push_back(pagmo::island(algo1,prob,1));
					a.push_back(pagmo::island(algo2,prob,20));
					a.push_back(pagmo::island(algo1,prob,1));
					a.push_back(pagmo::island(algo2,prob,20));
					a.push_back(pagmo::island(algo1,prob,1));
					a.push_back(pagmo::island(algo2,prob,20));
					a.evolve_t(10000);
					a.join();
					std::cout << "\tBest:" << a.get_island(0)->get_population().champion().f << std::endl;

					//log
					std::vector<double> x = a.get_island(0)->get_population().champion().x;
					double time = x[4] * Tmax + x[6] + x[10] * ((1 - x[4]) * Tmax - x[6]);
					myfile << "[" << target.get_name() << "] %" << "[" << time << "] " << a.get_island(0)->get_population().champion().f << " " << a.get_island(0)->get_population().champion().x << std::endl;
				}
			}
		} catch (value_error) {
			std::cout << "The End!!!" << std::endl;
			return 0;
		} catch (boost::bad_lexical_cast) {
			// An empty line might be present in MPCORB.DAT
		}
	}
	//close file
	mpcorbfile.close();
	return 0;
}
