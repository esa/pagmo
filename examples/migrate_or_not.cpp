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
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "../src/keplerian_toolbox/keplerian_toolbox.h"
#include "../src/pagmo.h"


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

	myfile << "\\documentclass{article}\n";
	myfile << "\\usepackage{xtab}\n";
	myfile << "\\begin{document}";
	myfile << "\\begin{xtabular}{lll}\n";

	//0 - Experiment parameters
	int number_of_islands = 3;
	int number_of_individuals = 20;
	//int evolution_time = 1000;
	int number_of_migrations = 2;

	/*
	std::vector<double> r1;
	r1.push_back(10);
	r1.push_back(5.3);
	std::vector<double> r2;
	r2.push_back(4);
	r2.push_back(1.5);
	std::vector<std::vector<double> > weights;
	weights.push_back(r1);
	weights.push_back(r2);
	*/
	
	std::vector<double> tmp(4,0);
	std::vector<std::vector<double> > w(4,tmp);
	w[0][0] = 0;
	w[0][1] = 1;
	w[0][2] = 100;
	w[0][3] = 1;
	w[1][0] = 1;
	w[1][1] = 0;
	w[1][2] = 1;
	w[1][3] = 100;
	w[2][0] = 100;
	w[2][1] = 1;
	w[2][2] = 0;
	w[2][3] = 1;
	w[3][0] = 1;
	w[3][1] = 100;
	w[3][2] = 1;
	w[3][3] = 0;
	
	std::vector<double> values(5,0);
	values[0] = 3;
	values[1] = 7;
	values[2] = 5;
	values[3] = 3;
	values[4] = 5;

	std::vector<double> weights(5,0);
	weights[0] = 10;
	weights[1] = 10;
	weights[2] = 1;
	weights[3] = 2;
	weights[4] = 3;

	double max_weight = 11;


	//1 - We instantiate the problems
	problem::cassini_1 prob1;
	problem::griewank prob2(5);
	problem::ackley prob3(5);
	problem::rastrigin prob4(5);
	problem::michalewicz prob5(5);
	problem::tsp prob6(w);
	problem::knapsack prob7(values, weights, max_weight);
	problem::sch prob8;
	problem::fon prob9;
	problem::pol prob10;
	problem::kur prob11;
	problem::zdt1 prob12;
	problem::zdt2 prob13;
	problem::zdt3 prob14;
	problem::zdt4 prob15;
	problem::zdt6 prob16;

	//2 - We instantiate the algorithms
	algorithm::de algo1(100);
	algorithm::sga algo2(100,0.8,0.05,1);
	algorithm::sa_corana algo3(2000,1,0.001);
	algorithm::pso algo4(100);
	algorithm::ihs algo5(2000);
	algorithm::bee_colony algo6(50);
	algorithm::cross_entropy algo7(100,0.3);
	algorithm::aco algo8(10);
	algorithm::nsga2 algo9(250,0.9,0.3);
	algorithm::firefly algo10(50);

	//b - We instantiate the topologies
	topology::unconnected topo1;
	topology::ring topo2;
	topology::fully_connected topo3;
	topology::watts_strogatz topo4;

	//3 - We build a container of algorithms
	std::vector<algorithm::base_ptr> algo;
	//algo.push_back(algo1.clone());
	//algo.push_back(algo2.clone());
	//algo.push_back(algo3.clone());
	//algo.push_back(algo4.clone());
	//algo.push_back(algo5.clone());
	//algo.push_back(algo6.clone());
	//algo.push_back(algo7.clone());
	//algo.push_back(algo8.clone());
	algo.push_back(algo9.clone());
	//algo.push_back(algo10.clone());

	//4 - And a container of problems
	std::vector<problem::base_ptr> prob;
	//prob.push_back(prob1.clone());
	//prob.push_back(prob2.clone());
	//prob.push_back(prob3.clone());
	//prob.push_back(prob4.clone());
	//prob.push_back(prob5.clone());
	//prob.push_back(prob6.clone());
	//prob.push_back(prob7.clone());
	prob.push_back(prob8.clone());
	prob.push_back(prob9.clone());
	//prob.push_back(prob10.clone());
	//prob.push_back(prob11.clone());
	//prob.push_back(prob12.clone());
	//prob.push_back(prob13.clone());
	//prob.push_back(prob14.clone());
	//prob.push_back(prob15.clone());
	//prob.push_back(prob16.clone());

	//5 - And a container of topologies
	std::vector<topology::base_ptr> topo;
	topo.push_back(topo1.clone());
	//topo.push_back(topo2.clone());
	//topo.push_back(topo3.clone());
	//topo.push_back(topo4.clone());

	for (unsigned int pr=0; pr<prob.size();++pr) {
		std::cout << std::endl << "Problem: " << prob[pr]->get_name() << std::endl;

		for (unsigned int al =0; al<algo.size(); ++al) {
			const std::string algo_name = ((al==algo.size()) ? std::string("Coop") : algo[al]->get_name());
			std::cout << algo_name << '\n' << '\n';
			std::cout << "\t\tMean" << "\t\tStd Deviation" << std::endl;
			myfile << "\\hline\n" << "\\multicolumn{3}{c}{" << prob[pr]->get_name() << ", " << algo_name << "}" << "\\\\ \n \\hline\n";

			for (unsigned int to=0; to<topo.size(); ++to) {

				archipelago a = pagmo::archipelago(*topo[to]);
				for (int i=0; i<number_of_islands; ++i) {
					if (al == algo.size())
						a.push_back(island(*algo[i%al],*prob[pr],number_of_individuals));
					else
						a.push_back(island(*algo[al],*prob[pr],number_of_individuals));
				}
				a.evolve(number_of_migrations);
				a.join();
				std::cout << topo[to]->get_name() << ":\t " << mean(a) << "\t" << std_dev(a,mean(a)) << std::endl;
				for(int i = 0; i < number_of_islands; ++i) {
					std::cout << "Island " << i << std::endl << a.get_island(i)->get_population().champion().human_readable() << std::endl;					
				}
				print_row(myfile,topo[to]->get_name(),mean(a),std_dev(a,mean(a)));
			}
		}
	}

	myfile << "\\end{xtabular}\n";
	myfile << "\\end{document}";
	return 0;
}
