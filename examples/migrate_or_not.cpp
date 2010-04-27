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
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include"../src/pagmo.h"

using namespace pagmo;
using namespace kep_toolbox;

double mean(archipelago a) {
	double retval = 0;
	for (archipelago::size_type i = 0; i< a.get_size(); ++i) {
		retval += a.get_island(i).get_population().champion().f[0];
	}
	return retval / a.get_size();
}

double std_dev(archipelago a, double mean) {
	double retval = 0;
	for (archipelago::size_type i = 0; i< a.get_size(); ++i) {
		retval += pow((a.get_island(i).get_population().champion().f[0] - mean),2);
	}
	return sqrt(retval / a.get_size());
}

int main()
{
	//1 - We instantiate the problems
	problem::levy5 prob1(10);
	problem::griewank prob2(10);
	problem::ackley prob3(10);
	problem::rastrigin prob4(10);

	//2 - We instantiate the algorithms
	algorithm::de algo1(100);
	algorithm::sga algo2(100,0.8,0.05,1);
	algorithm::sa_corana algo3(2000,1,0.001);
	algorithm::pso algo4(100);

	//3 - We build a container of algorithms
	std::vector<algorithm::base_ptr> algo;
	algo.push_back(algo1.clone());
	algo.push_back(algo2.clone());
	algo.push_back(algo3.clone());
	algo.push_back(algo4.clone());

	//4 - And a container of problems
	std::vector<problem::base_ptr> prob;
	prob.push_back(prob1.clone());
	prob.push_back(prob2.clone());
	prob.push_back(prob3.clone());
	prob.push_back(prob4.clone());

	for (int pr=0; pr<4;++pr) {
		std::cout << std::endl << "Problem: " << prob[pr]->get_name() << std::endl;

		for (int al =0; al<4; ++al) {
			std::cout << *algo[al] << '\n' << '\n';
			std::cout << "\t\tMean" << "\t\tStd Deviation" << std::endl;

			//We start with an unconnected topology (i.e. no migration)
			{
				archipelago a = pagmo::archipelago(topology::unconnected());
				for (int i=0; i<20; ++i) {
					a.push_back(island(*prob[pr],*algo[al],20));
				}
				a.evolve(10);
				a.join();
				std::cout << "No topology:\t" << mean(a) << "\t" << std_dev(a,mean(a)) << std::endl;
			}
			//We repeat with ring topology
			{
				archipelago a = pagmo::archipelago(topology::ring());
				for (int i=0; i<20; ++i) {
					a.push_back(island(*prob[pr],*algo[al],20));
				}
				a.evolve(10);
				a.join();
				std::cout << "Ring topology:\t" << mean(a) << "\t" << std_dev(a,mean(a)) << std::endl << std::endl;
			}
		}
	}
	return 0;
}
