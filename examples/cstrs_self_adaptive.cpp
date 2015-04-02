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

template<typename T>
std::string to_string( const T & value )
{
	// utiliser un flux de sortie pour créer la chaîne
	std::ostringstream oss;
	// écrire la valeur dans le flux
	oss << value;
	// renvoyer une string
	return oss.str();
}

double mean(archipelago a, problem::base_ptr original_problem) {
	double retval = 0;
	int count_feasible_arch = 0;
	for (archipelago::size_type i = 0; i< a.get_size(); ++i) {
		// test feasibility
		if(original_problem->feasibility_x(a.get_island(i)->get_population().champion().x)) {
			retval += a.get_island(i)->get_population().champion().f[0];
			count_feasible_arch++;
		}
	}
	if(count_feasible_arch != 0)
		return retval / count_feasible_arch;
	else
		return 0.;
}
double std_dev(archipelago a, double mean, problem::base_ptr original_problem) {
	double retval = 0;
	int count_feasible_arch = 0;
	for (archipelago::size_type i = 0; i< a.get_size(); ++i) {
		// test feasibility
		if(original_problem->feasibility_x(a.get_island(i)->get_population().champion().x)) {
			retval += pow((a.get_island(i)->get_population().champion().f[0] - mean),2);
			count_feasible_arch++;
		}
	}

	if(count_feasible_arch != 0)
		return sqrt(retval / count_feasible_arch);
	else
		return 0.;
}

int number_violated_constraints(const constraint_vector &c, problem::base_ptr original_problem)
{
	int viol = 0;
	constraint_vector::size_type c_dim = original_problem->get_c_dimension();

	for(constraint_vector::size_type i=0; i<c_dim; i++) {
		if (!original_problem->test_constraint(c,i)) {
			viol++;
		}
	}
	return viol;
}

double mean_violated_constraints(const constraint_vector &c, problem::base_ptr original_problem) {
	double viol = 0;

	constraint_vector::size_type c_dim = original_problem->get_c_dimension();
	problem::base::c_size_type number_of_eq_constraints =
			original_problem->get_c_dimension() -
			original_problem->get_ic_dimension();

	const std::vector<double> &c_tol = original_problem->get_c_tol();

	for(constraint_vector::size_type j=0; j<number_of_eq_constraints; j++) {
		viol += std::max(0.,(std::abs(c.at(j)) - c_tol.at(j)));
	}
	for(constraint_vector::size_type j=number_of_eq_constraints; j<c_dim; j++) {
		viol += std::max(0.,c.at(j));
	}

	viol /= c_dim;

	return viol;
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

problem::base_ptr get_constrained_prob(problem::base_ptr prob) {
	return prob;
}

algorithm::base_ptr get_constrained_algo(algorithm::base_ptr algo, int i) {
	switch(i) {
	case(0):
        return algorithm::cstrs_self_adaptive(*algo, 100).clone();
		break;
	case(1):
        return algorithm::cstrs_self_adaptive(*algo, 5000).clone();
		break;
	case(2):
        return algorithm::cstrs_self_adaptive(*algo, 20000).clone();
		break;
	default:
		return algo;
		break;
	}
}

std::string get_constrained_name(int ch) {
	switch(ch) {
	case(0):
        return "cstrs_self_adaptive_valid_gen_100";
		break;
	case(1):
        return "cstrs_self_adaptive_valid_gen_5000";
		break;
	case(2):
        return "cstrs_self_adaptive_valid_gen_20000";
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
	size_t number_of_islands = 20;
	size_t number_of_individuals = 70;
	// size_t function_evaluations = 100;
	size_t number_of_migrations = 1;
	// int gen = function_evaluations/number_of_individuals/number_of_migrations;

	//1 - We instantiate the problems and store the problems
	std::vector<problem::base_ptr> probs;
	for(int i=0; i<11; i++) {
		probs.push_back(problem::cec2006(i+1).clone());
	}

	//2 - We instantiate the algorithms
	std::vector<algorithm::base_ptr> algos;

	// avoiding elitism here
//	algos.push_back(pagmo::algorithm::sga(1 ,0.9, 0.03, 1000,
//										  algorithm::sga::mutation::RANDOM, 0.1,
//										  algorithm::sga::selection::ROULETTE,
//										  algorithm::sga::crossover::BINOMIAL).clone());

	//	algos.push_back(pagmo::algorithm::sga(1 ,0.9, 0.05, 1000,
	//										  algorithm::sga::mutation::RANDOM, 0.1,
	//										  algorithm::sga::selection::ROULETTE,
	//										  algorithm::sga::crossover::BINOMIAL).clone());

	//	algos.push_back(pagmo::algorithm::sga(1 ,0.9, 0.05, 1000,
	//										  algorithm::sga::mutation::RANDOM, 0.1,
	//										  algorithm::sga::selection::ROULETTE,
	//										  algorithm::sga::crossover::EXPONENTIAL).clone());

	//	algos.push_back(pagmo::algorithm::sga(1 ,0.9, 0.05, 1000,
	//										  algorithm::sga::mutation::GAUSSIAN, 0.1,
	//										  algorithm::sga::selection::ROULETTE,
	//										  algorithm::sga::crossover::BINOMIAL).clone());

	//algos.push_back(pagmo::algorithm::de(1).clone());

	algos.push_back(pagmo::algorithm::sga_gray(1,0.9,0.003,10000,
											   algorithm::sga_gray::mutation::UNIFORM,
											   algorithm::sga_gray::selection::ROULETTE,
											   algorithm::sga_gray::crossover::SINGLE_POINT).clone());

	//b - We instantiate the topologies
	std::vector<topology::base_ptr> topo;
	topo.push_back(topology::unconnected().clone());

	// for each constraints handling version:
	for(unsigned int ch=0; ch < 3; ch++) {

		//Open the output file
		std::ofstream myfile;
		myfile.open((std::string("pagmo_constraints_") + get_constrained_name(ch) + ".tex").c_str());
		myfile << std::setprecision(5);

		myfile << "\\documentclass{article}\n";
		myfile << "\\usepackage{array}\n";
		myfile << "\\usepackage{rotating}\n";
		myfile << "\\begin{document}\n";

		myfile << "\\clearpage\n";

		myfile << "\\begin{sidewaystable}\n";
		myfile << "\\centering\n";
		myfile << "\\begin{tabular}{| c | c | c | c | c | c | c | c | c | c |}" << "\n";
		myfile << "FES & Problem & Best know & Best & Median & Worst & $\\bar{v}$ & Mean champ & Std champ\\\\";

		for (unsigned int to=0; to<topo.size(); ++to) {
			std::cout << topo[to]->get_name() << std::endl;

			for (unsigned int pr=0; pr<probs.size();++pr) {
				myfile << "\\hline\n";

				//			if(pr == 0) {
				//				myfile << gen;
				//			}

				myfile << " & " << probs[pr]->get_name() << " & " << probs[pr]->get_best_f()[0][0];

				// we generate the constrained problem
				problem::base_ptr constrained_problem = get_constrained_prob(probs[pr]);

				std::cout << std::endl << "Problem: " << constrained_problem->get_name()<<std::endl;

				for (unsigned int al=0; al<algos.size(); ++al) {
					// we generate the constrained algorithm
					algorithm::base_ptr constrained_algorithm = get_constrained_algo(algos[al], ch);

					pagmo::archipelago a = pagmo::archipelago(*topo[to]);

					for (size_t i=0; i<number_of_islands; ++i) {
						a.push_back(island(*constrained_algorithm, *constrained_problem, number_of_individuals));
					}
					a.evolve(number_of_migrations);
					a.join();

					// get the populations and fill a vector with it:
					std::vector<population::individual_type> final_population;

					// store all the individuals in a vector that will be used to sort the population and so on
					for (archipelago::size_type i=0; i<a.get_size(); ++i) {
						const population &current_population = a.get_island(i)->get_population();

						for (population::size_type j=0; j<(current_population.size()); j++) {
							final_population.push_back(current_population.get_individual(j));
						}
					}

					// sorting the population
					population::size_type final_population_size = final_population.size();

					// in practice, sorting the final population should be:
					// 1. feasible solutions in front of infeasible solutions
					// 2. feasible solutions according to their function error f(x) - f(x*)
					// 3. infeasible solutions according to their mean value of the violation of all constraints
					// instead, we prefer use the compare_constraints_impl of the problem
					for(population::size_type i=0; i<(final_population_size-1); i++) {
						for(population::size_type j=i+1; j<final_population_size; ++j) {
							if(probs[pr]->compare_fc(final_population.at(j).cur_f, final_population.at(j).cur_c,
													 final_population.at(i).cur_f, final_population.at(i).cur_c) ) {
								//swap individuals
								population::individual_type temp_ind = final_population[i];
								final_population[i] = final_population[j];
								final_population[j] = temp_ind;
							}
						}
					}

					// show the best solution:
					// std::cout << final_population.at(0);

					population::size_type median = 0;
					// finds the median individual in the current population
					// (gets the first one in the second middle if even numbered)
					if (final_population_size % 2 == 0) {
						median = (final_population_size + 1) / 2;
					}
					else {
						median = final_population_size / 2;
					}

					// number of violated constraints at the median solution:
					int viol_median = number_violated_constraints(final_population.at(median).cur_c, probs[pr]);
					// number of violated constraints at the best solution:
					int viol_best = number_violated_constraints(final_population.at(0).cur_c, probs[pr]);
					// number of violated constraints at the worst solution:
					int viol_worst = number_violated_constraints(final_population.at(final_population_size-1).cur_c, probs[pr]);

					double mean_viol_median = mean_violated_constraints(final_population.at(median).cur_c, probs[pr]);

					double mean_among_feas_champions = mean(a, probs[pr]);
					double std_dev_among_feas_champions = std_dev(a,mean_among_feas_champions, probs[pr]);

					myfile << " & " << final_population.at(0).cur_f[0] << "(" << viol_best << ")" << " & "
						   << final_population.at(median).cur_f[0] << "(" << viol_median << ")"  <<  " & "
						   << final_population.at(final_population_size-1).cur_f[0] << "(" << viol_worst << ")" << " & "
						   << mean_viol_median << " & "
						   << mean_among_feas_champions << " & "
						   << std_dev_among_feas_champions;
					myfile << "\\\\";
					std::cout << ":\t " << mean(a, probs[pr]) << "\t" << std_dev(a,mean(a, probs[pr]), probs[pr]) << std::endl;
				}
			}
		}
		myfile << "\\hline\n";
		myfile << "\\end{tabular}\n";
		myfile << "\\end{sidewaystable}\n";
		myfile << "\\clearpage\n";

		myfile << "\\end{document}";
		myfile.close();
	}
	return 0;
}


