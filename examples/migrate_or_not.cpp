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

/**
DESCRITPION: This example tests several algorithms, topologies and problems and prints out a latex
table summarizing all results. The purpose is to show the effect of the generalized migration operator
and of the topology in the archipelago.

CPU TIME: each optimization open 20 threads. The example completes in minutes.
*/

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

void print_row(std::ostream &f, const std::string &topo_name, const double &mean, const double &dev)
{
	f << topo_name << " & " << mean << " & " << dev << "\\\\";
}

int main()
{
	//Open the output file
	ofstream myfile;
	myfile.open((std::string("pagmo_") +
		boost::lexical_cast<std::string>(rng_generator::get<rng_uint32>()()) + ".tex").c_str());

	myfile << "\\documentclass{article}\n";
	myfile << "\\usepackage{xtab}\n";
	myfile << "\\begin{document}";
	myfile << "\\begin{xtabular}{lll}\n";

	//0 - Experiment parameters
	int number_of_islands = 20;
	int number_of_individuals = 20;
	int evolution_time = 1000;

	//1 - We instantiate the problems
	problem::cassini_1 prob1;
	problem::griewank prob2(50);
	problem::ackley prob3(50);
	problem::rastrigin prob4(50);

	//2 - We instantiate the algorithms
	algorithm::de algo1(100);
	algorithm::sga algo2(100,0.8,0.05,1);
	algorithm::sa_corana algo3(2000,1,0.001);
	algorithm::pso algo4(100);
	algorithm::ihs algo5(2000);

	//b - We instantiate the topologies
	topology::unconnected topo1;
	topology::ring topo2;
	topology::fully_connected topo3;
	topology::watts_strogatz topo4;

	//3 - We build a container of algorithms
	std::vector<algorithm::base_ptr> algo;
	algo.push_back(algo1.clone());
	algo.push_back(algo2.clone());
	algo.push_back(algo3.clone());
	algo.push_back(algo4.clone());
	algo.push_back(algo5.clone());

	//4 - And a container of problems
	std::vector<problem::base_ptr> prob;
	prob.push_back(prob1.clone());
	prob.push_back(prob2.clone());
	prob.push_back(prob3.clone());
	prob.push_back(prob4.clone());

	//5 - And a container of topologies
	std::vector<topology::base_ptr> topo;
	topo.push_back(topo1.clone());
	topo.push_back(topo2.clone());
	topo.push_back(topo3.clone());
	topo.push_back(topo4.clone());

	for (unsigned int pr=0; pr<prob.size();++pr) {
		std::cout << std::endl << "Problem: " << prob[pr]->get_name() << std::endl;

		for (unsigned int al =0; al<algo.size()+1; ++al) {
			const std::string algo_name = ((al==algo.size()) ? std::string("Coop") : algo[al]->get_name());
			std::cout << algo_name << '\n' << '\n';
			std::cout << "\t\tMean" << "\t\tStd Deviation" << std::endl;
			myfile << "\\hline\n" << "\\multicolumn{3}{c}{" << prob[pr]->get_name() << ", " << algo_name << "}" << "\\\\ \n \\hline\n";

			for (unsigned int to=0; to<topo.size(); ++to) {

				archipelago a = pagmo::archipelago(*topo[to]);
				for (int i=0; i<number_of_islands; ++i) {
					if (al == algo.size())
						a.push_back(island(*prob[pr],*algo[i%al],number_of_individuals));
					else
						a.push_back(island(*prob[pr],*algo[al],number_of_individuals));
				}
				a.evolve_t(evolution_time);
				a.join();
				std::cout << topo[to]->get_name() << ":\t " << mean(a) << "\t" << std_dev(a,mean(a)) << std::endl;
				print_row(myfile,topo[to]->get_name(),mean(a),std_dev(a,mean(a)));
			}
		}
	}

	myfile << "\\end{xtabular}\n";
	myfile << "\\end{document}";
	return 0;
}
