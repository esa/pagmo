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
table summarizing all results. The purpose is to show the effect of the constraints handling techniques.
*/

using namespace pagmo;

double best(archipelago a, problem::base_ptr original_problem) {
	double retval = boost::numeric::bounds<double>::highest();
	for (archipelago::size_type i = 0; i< a.get_size(); ++i) {
		// test feasibility
		if(original_problem->feasibility_x(a.get_island(i)->get_population().champion().x))
			retval = std::min(retval, a.get_island(i)->get_population().champion().f[0]);
	}
	return retval;
}
double worst(archipelago a, problem::base_ptr original_problem) {
	double retval = - boost::numeric::bounds<double>::highest();
	for (archipelago::size_type i = 0; i< a.get_size(); ++i) {
		// test feasibility
		if(original_problem->feasibility_x(a.get_island(i)->get_population().champion().x))
			retval = std::max(retval, a.get_island(i)->get_population().champion().f[0]);
	}
	return retval;
}
double mean(archipelago a, problem::base_ptr original_problem) {
	double retval = 0;
	for (archipelago::size_type i = 0; i< a.get_size(); ++i) {
		// test feasibility
		if(original_problem->feasibility_x(a.get_island(i)->get_population().champion().x))
			retval += a.get_island(i)->get_population().champion().f[0];
	}
	return retval / a.get_size();
}
double std_dev(archipelago a, double mean, problem::base_ptr original_problem) {
	double retval = 0;
	for (archipelago::size_type i = 0; i< a.get_size(); ++i) {
		// test feasibility
		if(original_problem->feasibility_x(a.get_island(i)->get_population().champion().x))
			retval += pow((a.get_island(i)->get_population().champion().f[0] - mean),2);
	}
	return sqrt(retval / a.get_size());
}

std::string get_solutions(archipelago a) {
	std::ostringstream sol;
	int sol_size = a.get_island(0)->get_population().champion().x.size();
	for (archipelago::size_type i = 0; i< a.get_size(); ++i) {
		sol << "island " << i << ": (";
		for(int j = 0; j < sol_size; ++j) {
			sol << a.get_island(i)->get_population().champion().x[j] << ",";
		}
		sol << ")" << std::endl;
	}
	return sol.str();
}

problem::base_ptr get_constrained_prob(problem::base_ptr prob, int i) {
	switch(i) {
	case(0):
		return problem::death_penalty(*prob,problem::death_penalty::SIMPLE).clone();
		break;
	case(1):
		return problem::death_penalty(*prob,problem::death_penalty::KURI).clone();
		break;
	case(2):
		return problem::con2mo(*prob,problem::con2mo::OBJ_CSTRS).clone();
		break;
	case(3):
		return problem::con2mo(*prob,problem::con2mo::OBJ_CSTRSVIO).clone();
		break;
	case(4):
		return problem::con2mo(*prob,problem::con2mo::OBJ_EQVIO_INEQVIO).clone();
		break;
	default:
		return prob;
		break;
	}
}

algorithm::base_ptr get_constrained_algo(algorithm::base_ptr algo, int i) {
	switch(i) {
//	case(3):
//		return algorithm::self_adaptive(*algo).clone();
//		break;
	default:
		return algo;
		break;
	}
}

std::string get_constrained_name(int i) {
	switch(i) {
	case(0):
		return "Death_penalty";
		break;
	case(1):
		return "Death_penalty_Kuri";
		break;
	case(2):
		return "con2mo_obj_cstrs";
		break;
	case(3):
		return "con2mo_obj_cstrsvio";
		break;
	case(4):
		return "con2mo_obj_eqvio_ineqvio";
		break;
	case(5):
		return "Self adaptive";
		break;
	default:
		return "No_constraint_technique";
		break;
	}
}

void print_row(std::ostream &f, const std::string &topo_name, const double &mean, const double &dev)
{
	f << topo_name << " & " << mean << " & " << dev << "\\\\";
}

std::string underscore_to_space(std::string text) {
	for(std::string::iterator it = text.begin(); it != text.end(); ++it) {
		if(*it == '_') {
			*it = ' ';
		}
	}
	return text;
}

int main()
{
	std::cout << std::setprecision(5);

	//0 - Experiment parameters

	// in CEC2006, the max number of function evaluation = 5,000; 50,000; 500,000
	// for each run
	size_t number_of_islands = 25;
	size_t number_of_individuals = 60; // was 60
	size_t function_evaluations = 5000;
	size_t number_of_migrations = 1;

	//1 - We instantiate the problems and store the problems
	std::vector<problem::base_ptr> probs;
	for(int i=0; i<24; i++) {
		probs.push_back(problem::cec2006(i+1).clone());
	}

	//2 - We instantiate the algorithms
	std::vector<algorithm::base_ptr> algos;
	int gen = function_evaluations/number_of_individuals/number_of_migrations;
	algos.push_back(algorithm::vega(gen).clone());
	algos.push_back(algorithm::nsga2(gen).clone());
	algos.push_back(algorithm::bee_colony(gen/2.).clone());
	algos.push_back(algorithm::cmaes(gen).clone());
	algos.push_back(algorithm::cs(gen*10,0.02,0.3,0.3).clone());
	algos.push_back(algorithm::de(gen).clone());
	algos.push_back(algorithm::de_1220(gen).clone());
	algos.push_back(algorithm::ihs(gen*number_of_individuals).clone());
	algos.push_back(algorithm::jde(gen).clone());
	//    algos.push_back(algorithm::mbh(algorithm::de(gen),2,0.03).clone());
	algos.push_back(algorithm::mde_pbx(gen).clone());
	algos.push_back(algorithm::monte_carlo(gen).clone());
	//    algos.push_back(algorithm::ms(algorithm::monte_carlo(gen),5).clone());
	//    algos.push_back(algorithm::null().clone());
	algos.push_back(algorithm::pso(gen).clone());
	//    algos.push_back(algorithm::pso_generational(gen,0.5,0.5,0.5,0.5,3,3,3).clone());
	algos.push_back(algorithm::sa_corana(gen*number_of_individuals).clone());
	algos.push_back(algorithm::sga(gen).clone());

	//b - We instantiate the topologies
	std::vector<topology::base_ptr> topo;
	topo.push_back(topology::unconnected().clone());
	// topo.push_back(topology::fully_connected().clone());
	// topo.push_back(topology::ring().clone());

	std::vector<double> stored_best_optimum(probs.size(), boost::numeric::bounds<double>::highest());
	std::vector<std::string> stored_best_algorithm(probs.size(), "");
	std::vector<std::string> stored_best_constraint_technique(probs.size(), "");

	// for each constraints handling technique:
	for(unsigned int ch=2; ch < 5; ch++) {

		//Open the output file
		std::ofstream myfile;
		myfile.open((std::string("pagmo_constraints_") + get_constrained_name(ch) + ".tex").c_str());
		myfile << std::setprecision(5);

		myfile << "\\documentclass{article}\n";
		myfile << "\\usepackage{array}\n";
		myfile << "\\usepackage{rotating}\n";
		myfile << "\\begin{document}\n";
		for (unsigned int pr=0; pr<probs.size();++pr) {
			myfile << "\\begin{sidewaystable}\n";
			myfile << "\\centering\n";
			myfile << "\\begin{tabular}{l|l|l|l|l|l|l}" << "\n";
			myfile << probs[pr]->get_name();
			myfile << "& Best & Worst & Mean & Std & Feasible rate & Success Rate \\\\";
			myfile << "\\hline\n";

			// we generate the constrained problem
			problem::base_ptr constrained_problem = get_constrained_prob(probs[pr], ch);

			std::cout << std::endl << "Problem: " << constrained_problem->get_name()<<std::endl;

			for (unsigned int to=0; to<topo.size(); ++to) {

				std::cout << topo[to]->get_name() << std::endl;

				for (unsigned int al=0; al<algos.size(); ++al) {
					// we generate the constrained algorithm
					algorithm::base_ptr constrained_algorithm = get_constrained_algo(algos[al], ch);

					pagmo::archipelago a = pagmo::archipelago(*topo[to]);

					for (size_t i=0; i<number_of_islands; ++i) {
						a.push_back(island(*constrained_algorithm, *constrained_problem, number_of_individuals));
					}
					a.evolve(number_of_migrations);
					a.join();

					myfile << "$" << constrained_algorithm->get_name() << "$";
					myfile << " & " << best(a, probs[pr]) << " & " << worst(a, probs[pr]) << " & " << mean(a, probs[pr]) << " & " << std_dev(a,mean(a, probs[pr]), probs[pr]) << "& &";
					myfile << "\\\\";
					std::cout << ":\t " << mean(a, probs[pr]) << "\t" << std_dev(a,mean(a, probs[pr]), probs[pr]) << std::endl;

					// store the best solution for the recap
					if(best(a, probs[pr]) < stored_best_optimum[pr]) {
						stored_best_optimum[pr] = best(a, probs[pr]);
						stored_best_algorithm[pr] = constrained_algorithm->get_name();
						stored_best_constraint_technique[pr] = get_constrained_name(ch);
					}
				}
			}
			myfile << "\\hline\n";
			myfile << "\\end{tabular}\n";
			myfile << "\\end{sidewaystable}\n";
			myfile << "\\clearpage\n";
		}

		myfile << "\\end{document}";
		myfile.close();
	}

	// recap
	std::ofstream recap_file;
	recap_file.open(std::string("pagmo_constraints_recap.tex").c_str());
	recap_file << std::setprecision(5);

	recap_file << "\\documentclass{article}\n";
	recap_file << "\\usepackage{array}\n";
	recap_file << "\\usepackage{rotating}\n";
	recap_file << "\\begin{document}\n";
	recap_file << "\\begin{sidewaystable}\n";
	recap_file << "\\centering\n";
	recap_file << "\\begin{tabular}{l|l|l|l|l|l}" << "\n";
	recap_file << "Prob & Nbr inequality cstr & Nbr equality cstr & Best known optimum & Optimum found & Best algorithm \\\\";
	recap_file << "\\hline\n";
	for (unsigned int pr=0; pr<probs.size();++pr) {
		recap_file << probs[pr]->get_name() << " & " << " & " << " & " << probs[pr]->get_best_f()[0][0] << " & " << stored_best_optimum[pr] <<  " & " << underscore_to_space(stored_best_algorithm[pr]) << "\\\\";

		recap_file << "" << " & " << "" << " & " << "" << " & " << "" << " & " << "" << " & " << underscore_to_space(stored_best_constraint_technique[pr]) << "\\\\";
		recap_file << "\\hline\n";
	}
	recap_file << "\\end{tabular}\n";
	recap_file << "\\end{sidewaystable}\n";
	recap_file << "\\end{document}";

	return 0;
}


