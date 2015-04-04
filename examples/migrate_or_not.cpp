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
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "../src/pagmo.h"


/**
DESCRITPION: This example tests several algorithms, topologies and problems and prints out a latex
table summarizing all results. The purpose is to show the effect of the generalized migration operator
and of the topology in the archipelago.

CPU TIME: each optimization open 20 threads. The example completes in minutes.
*/

using namespace pagmo;

double mean(archipelago a) {
	double retval = 0;
	for (archipelago::size_type i = 0; i< a.get_size(); ++i) {
		retval += a.get_island(i)->get_population().champion().f[0];
	}
	return retval / a.get_size();
}

double std_dev(archipelago a, double mean) {
	double retval = 0;
	for (archipelago::size_type i = 0; i< a.get_size(); ++i) {
		retval += pow((a.get_island(i)->get_population().champion().f[0] - mean),2);
	}
	return sqrt(retval / a.get_size());
}

std::string getSolutions(archipelago a) {
	std::ostringstream sol;
	int solSize = a.get_island(0)->get_population().champion().x.size();
	for (archipelago::size_type i = 0; i< a.get_size(); ++i) {
		sol << "island " << i << ": (";
		for(int j = 0; j < solSize; ++j) {
			sol << a.get_island(i)->get_population().champion().x[j] << ",";
		}
		sol << ")" << std::endl;
	}	
	return sol.str();
}


void print_row(std::ostream &f, const std::string &topo_name, const double &mean, const double &dev)
{
	f << topo_name << " & " << mean << " & " << dev << "\\\\";
}

int main()
{
	//Open the output file
	std::ofstream myfile;
	myfile.open((std::string("pagmo_") +
		boost::lexical_cast<std::string>(rng_generator::get<rng_uint32>()()) + ".tex").c_str());
	myfile << std::setprecision(5);
	myfile << "\\documentclass{article}\n";
	myfile << "\\usepackage{xtab}\n";
	myfile << "\\usepackage{rotating}\n";
	myfile << "\\begin{document}\n";
	myfile << "\\begin{sidewaystable}\n";
	myfile << "\\begin{xtabular}{l|ll|ll|ll|ll|ll|ll}\n";

	//0 - Experiment parameters
	int number_of_islands = 8;
	int number_of_individuals = 20;
	int function_evaluations = 100;
	int number_of_migrations = 1;


	//1 - We instantiate the problems
	problem::griewank prob1(10);
	problem::ackley prob2(10);
	problem::rastrigin prob3(10);
	problem::rosenbrock prob4(10);
	problem::schwefel prob5(10);
	problem::levy5 prob6(10);
	problem::cassini_1 prob7;
	problem::cassini_2 prob8;
	problem::messenger_full prob9;
	problem::rosetta prob10;


	//2 - We instantiate the algorithms
	int gen = function_evaluations/number_of_individuals/number_of_migrations;
	algorithm::de algo1(gen);
	algorithm::sga algo2(gen,0.8,0.05,1);
	algorithm::sa_corana algo3(function_evaluations/number_of_migrations,1,0.001);
	algorithm::ihs algo4(function_evaluations/number_of_migrations);
	algorithm::bee_colony algo5(gen/2);
	algorithm::pso algo6(gen);

	//b - We instantiate the topologies
	topology::unconnected topo1;
	topology::fully_connected topo2;
	topology::ring topo3;

//	topology::watts_strogatz topo4;

	//3 - We build a container of algorithms
	std::vector<algorithm::base_ptr> algo;
	algo.push_back(algo1.clone());
	algo.push_back(algo2.clone());
	algo.push_back(algo3.clone());
	algo.push_back(algo4.clone());
	algo.push_back(algo5.clone());
	algo.push_back(algo6.clone());

	//4 - And a container of problems
	std::vector<problem::base_ptr> prob;
	prob.push_back(prob1.clone());
	prob.push_back(prob2.clone());
	prob.push_back(prob3.clone());
	prob.push_back(prob4.clone());
	prob.push_back(prob5.clone());
	prob.push_back(prob6.clone());
	prob.push_back(prob7.clone());
	prob.push_back(prob8.clone());
	prob.push_back(prob9.clone());
	prob.push_back(prob10.clone());

	//5 - And a container of topologies
	std::vector<topology::base_ptr> topo;
	topo.push_back(topo1.clone());
	topo.push_back(topo2.clone());
	topo.push_back(topo3.clone());


	for (unsigned int al =0; al<algo.size(); ++al) {
		myfile  << " & \\multicolumn{2}{c}{" << algo[al]->get_name() << "} ";
	}
	myfile << "\\\\";
	std::cout << std::setprecision(5);

	for (unsigned int pr=0; pr<prob.size();++pr) {
		std::cout << std::endl << "Problem: " << prob[pr]->get_name()<<std::endl;
		myfile << "\\hline\n" << "\\multicolumn{13}{c}{" << prob[pr]->get_name() << "}" <<  std::endl;
		myfile << std::endl << " \\\\ \\hline" << std::endl;
		for (unsigned int to=0; to<topo.size(); ++to) {
			myfile << topo[to]->get_name();
			std::cout << topo[to]->get_name() << std::endl;
			for (unsigned int al =0; al<algo.size(); ++al) {
			pagmo::archipelago a = pagmo::archipelago(*topo[to]);
				for (int i=0; i<number_of_islands; ++i) {
					a.push_back(island(*algo[al],*prob[pr],number_of_individuals));
				}
				a.evolve(number_of_migrations);
				a.join();
				myfile << " & " << mean(a) << " & " << std_dev(a,mean(a));
				std::cout << algo[al]->get_name() << ":\t " << mean(a) << "\t" << std_dev(a,mean(a)) << std::endl;
			}
			myfile << "\\\\" << std::endl;
		}
	}

	myfile << "\\hline\n";
	myfile << "\\end{xtabular}\n";
	myfile << "\\end{sidewaystable}\n";
	myfile << "\\end{document}";
	return 0;
}
