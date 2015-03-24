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


#include <string>
#include <vector>
#include <algorithm>

#include <boost/math/special_functions/binomial.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../archipelago.h"
#include "../island.h"
#include "../population.h"
#include "../topology/fully_connected.h"
#include "../topology/custom.h"
#include "../topology/watts_strogatz.h"
#include "../problem/decompose.h"
#include "../util/discrepancy.h"
#include "../util/neighbourhood.h"
#include "../migration/worst_r_policy.h"
#include "../migration/best_s_policy.h"
#include "../types.h"
#include "base.h"
#include "pade.h"

namespace pagmo { namespace algorithm {
/// Constructor
 /**
 * Constructs a PaDe algorithm
 *
 * @param[in] gen Number of generations to evolve.
 * @param[in] threads the amounts of threads that will be used
 * @param[in] method the decomposition method to use (Weighted, Tchebycheff or BI)
 * @param[in] solver the algorithm to solve the single objective problems.
 * @param[in] T the size of the population on each subproblem (must be an even number)
 * @param[in] weight_generation the method to generate the weight vectors (RANDOM, GRID or LOW-DISCREPANCY)
 * @param[in] z the reference point used for decomposition (with Tchebycheff and BI)
 *
 * @throws value_error if gen is negative, weight_generation is not sane
 * @see pagmo::problem::decompose::method_type
 */
pade::pade(int gen, unsigned int threads, pagmo::problem::decompose::method_type method,
		   const pagmo::algorithm::base & solver, population::size_type T, weight_generation_type weight_generation,
		   const fitness_vector &z)
	  :base(),
	  m_gen(gen),
	  m_threads(threads),
	  m_method(method),
	  m_solver(solver.clone()),
	  m_T(T),
	  m_weight_generation(weight_generation),
	  m_z(z)
{
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}

	//0 - Check whether method is implemented
	if(m_weight_generation != RANDOM && m_weight_generation != GRID && m_weight_generation != LOW_DISCREPANCY) {
		pagmo_throw(value_error,"non existing weight generation method");
	}
}

/// Copy constructor. Performs a deep copy. Necessary as a pointer to a base algorithm is here contained
pade::pade(const pade &algo):
	  base(algo),
	  m_gen(algo.m_gen),
	  m_threads(algo.m_threads),
	  m_method(algo.m_method),
	  m_solver(algo.m_solver->clone()),
	  m_T(algo.m_T),
	  m_weight_generation(algo.m_weight_generation),
	  m_z(algo.m_z)
{}

/// Clone method.
base_ptr pade::clone() const
{
	return base_ptr(new pade(*this));
}

//Recursive function building all m-ple of elements of X summing to s
void pade::reksum(std::vector<std::vector<double> > &retval,
		   const std::vector<unsigned int>& X,
		   unsigned int m,
		   unsigned int s,
		   std::vector<double> eggs) const {

	if (m==1) {
		if (std::find(X.begin(),X.end(),s) == X.end()) { //not found
			return;
		} else {
			eggs.push_back(s);
			retval.push_back(eggs);
		}
	} else {
		for (unsigned int i=0; i<X.size(); ++i) {
			eggs.push_back(X[i]);
			reksum(retval,X,m-1,s-X[i],eggs);
			eggs.pop_back();
		}
	}
}

/// Generates the weights used in the problem decomposition
/**
 * Run the PaDe algorithm for the number of generations specified in the constructors.
 *
 * @param[in] n_f diemension of the fitness space
 * @param[in] n_w number of weights to be produced
 */
 
std::vector<fitness_vector> pade::generate_weights(const unsigned int n_f, const unsigned int n_w) const {

	// Sanity check
	if (n_f > n_w) {
		pagmo_throw(value_error,"To allow weight be generated correctly the number of weights must be strictly larger than the number of objectives");
	}

	// Definition of useful probability distributions
	boost::uniform_real<double> uniform(0.0,1.0);
	boost::variate_generator<boost::lagged_fibonacci607 &, boost::uniform_real<double> > r_dist(m_drng,uniform);

	std::vector<fitness_vector> retval;
	if(m_weight_generation == GRID) {
			//find the largest H resulting in a population smaller or equal to NP
			unsigned int H;
			if (n_f == 2) {
				H = n_w-1;
			} else if (n_f == 3) {
				H = floor(0.5 * (sqrt(8*n_w + 1) - 3));
			} else {
				H = 1;
				while(boost::math::binomial_coefficient<double>(H+n_f-1, n_f-1) <= n_w) {
					++H;
				}
				H--;
			}

			// We check that NP equals the population size resulting from H
			if (fabs(n_w-(boost::math::binomial_coefficient<double>(H+n_f-1, n_f-1))) > 1E-8) {
				std::ostringstream error_message;
				error_message << "Invalid population size. Select " << boost::math::binomial_coefficient<double>(H+n_f-1, n_f-1)
						<< " or " << boost::math::binomial_coefficient<double>(H+1+n_f-1, n_f-1)
						<< ".";
				pagmo_throw(value_error,error_message.str());
			}
	
			// We generate the weights
			std::vector<unsigned int> range;
			for (unsigned int i=0; i<H+1;++i) {
				range.push_back(i);
			}
			reksum(retval, range, n_f, H);
			for(unsigned int i=0; i< retval.size(); ++i) {
				for(unsigned int j=0; j< retval[i].size(); ++j) {
					retval[i][j] /= H;
				}
			}
	
		} else if(m_weight_generation == LOW_DISCREPANCY) {
			for(unsigned int i = 0; i< n_f; ++i) {
				retval.push_back(fitness_vector(n_f,0.0));
				retval[i][i] = 1.0;
			}
			pagmo::util::discrepancy::simplex generator(n_f,1);
			for(unsigned int i = n_f; i <n_w; ++i) {
				retval.push_back(generator());
			}
	
		} else if(m_weight_generation == RANDOM) {
			pagmo::util::discrepancy::project_2_simplex projection(n_f);
			for (unsigned int i = 0; i<n_w; ++i) {
				fitness_vector dummy(n_f-1,0.0);
				for(unsigned int j = 0; j <n_f-1; ++j) {
					dummy[j] = r_dist();
				}
				retval.push_back(projection(dummy));
			}
		}
		return retval;
 }

/// Evolve implementation.
/**
 * Run the PaDe algorithm for the number of generations specified in the constructors.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */
void pade::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const population::size_type NP = pop.size();

	// And make some sanity checks
	if ( prob.get_f_dimension() < 2 ) {
		pagmo_throw(value_error, "The problem is not multiobjective, try some other algorithm than PaDE");
	}
	
	if ( m_T > NP-1 ) {
		pagmo_throw(value_error, "Too many neighbours specified");
	}

	// Get out if there is nothing to do.
	if (m_gen == 0) {
		return;
	}

	decision_vector dummy(pop.get_individual(0).cur_x.size(),0); //used for initialisation purposes
	std::vector<decision_vector> X(NP,dummy); //set of population chromosomes

	// Copy the population chromosomes into X
	for ( population::size_type i = 0; i<NP; i++ ) {
		X[i]	=	pop.get_individual(i).cur_x;
	}

	// Generate the weights for the NP decomposed problems
	std::vector<fitness_vector> weights = generate_weights(prob.get_f_dimension(), NP);
	
	// We compute, for each weight vector, the neighbouring ones (this will form the topology later on)
	std::vector<std::vector<population::size_type> > indices;
	pagmo::util::neighbourhood::euclidian::compute_neighbours(indices, weights);

	// Create the archipelago of NP islands:
	// each island in the archipelago solves a different single-objective problem.
	// We use here the broadcast migration model. This will force, at each migration,
	// to have individuals from all connected island to be inserted.
	pagmo::archipelago arch(pagmo::archipelago::broadcast);

	// Sets random number generators of the archipelago using the algorithm urng to obtain
	// a deterministic behaviour upon copy.
	arch.set_seeds(m_urng());
	
	// Best individual will be selected for migration
	const pagmo::migration::best_s_policy  selection_policy;
	
	// As m_T neighbours are connected, we replace m_T individuals on the island
	const pagmo::migration::worst_r_policy replacement_policy(m_T);

	//We create all the decomposed problems (one for each individual)
	std::vector<pagmo::problem::base_ptr> problems_vector;
	for(pagmo::population::size_type i=0; i<NP;++i) {
		problems_vector.push_back(pagmo::problem::decompose(prob, m_method,weights[i],m_z).clone());
	}

	//We create a pseudo-random permutation of the problem indexes
	std::vector<population::size_type> shuffle(NP);
	for(pagmo::population::size_type i=0; i < NP; ++i) {
			shuffle[i] = i;
	}
	boost::uniform_int<int> pop_idx(0,NP-1);
	boost::variate_generator<boost::mt19937 &, boost::uniform_int<int> > p_idx(m_urng,pop_idx);
	std::random_shuffle(shuffle.begin(), shuffle.end(), p_idx);

	//We assign each problem to the individual which has minimum fitness on that problem
	//This allows greater performance .... check without on dtlz2 for example.
	std::vector<int> assignation_list(NP); //problem i is assigned to the individual assignation_list[i]
	std::vector<bool> selected_list(NP,false);	//keep track of the individuals already assigned to a problem
	fitness_vector dec_fit(1);	//temporary stores the decomposed fitness
	for(pagmo::population::size_type i=0; i<NP;++i) { //for each problem i, select an individual j
		unsigned int j = 0;
		while(selected_list[j]) j++; //get to the first not already selected individual

		dynamic_cast<const pagmo::problem::decompose &>(*problems_vector[shuffle[i]]).compute_decomposed_fitness(dec_fit, pop.get_individual(j).cur_f);
		double minFit = dec_fit[0];
		int minFitPos = j;

		for(;j < NP; ++j) { //find the minimum fitness individual for problem i
			if(!selected_list[j]) { //just consider individuals which have not been selected already
				dynamic_cast<const pagmo::problem::decompose &>(*problems_vector[shuffle[i]]).compute_decomposed_fitness(dec_fit, pop.get_individual(j).cur_f);
				if(dec_fit[0] < minFit) {
					minFit = dec_fit[0];
					minFitPos = j;
				}
			}
		}
		assignation_list[shuffle[i]] = minFitPos;
		selected_list[minFitPos] = true;
	}

	for(pagmo::population::size_type i=0; i<NP;++i) { //for each island/problem i
		pagmo::population decomposed_pop(*problems_vector[i], 0, m_urng()); //Create a population for each decomposed problem

		//Set the individuals of the new population as one individual of the original population
		// (according to assignation_list) plus m_T neighbours individuals
		if(m_T < NP-1) {
			decomposed_pop.push_back(X[assignation_list[i]]); //assign to the island the correct individual according to the assignation list
			for(pagmo::population::size_type  j = 1; j <= m_T; ++j) { //add the neighbours
				decomposed_pop.push_back(X[assignation_list[indices[i][j]]]); //add the individual assigned to the island indices[i][j]
			}
		} else { //complete topology
			for(pagmo::population::size_type  j = 0 ; j < NP; ++j) {
				decomposed_pop.push_back(X[j]);
			}
		}
		arch.push_back(pagmo::island(*m_solver,decomposed_pop, selection_policy, replacement_policy));
	}

	topology::custom topo;
	if(m_T >= NP-1) {
		topo = topology::fully_connected();
	} else {
		for(unsigned int i = 0; i < NP; ++i) {
			topo.push_back();
		}
		for(unsigned int i = 0; i < NP; ++i) { //connect each island with the T closest neighbours
			for(unsigned int j = 1; j <= m_T; ++j) { //start from 1 to avoid to connect with itself
				topo.add_edge(i,indices[i][j]);
			}
		}
	}
	arch.set_topology(topo);


	//Evolve the archipelago for m_gen generations
	if(m_threads >= NP) { //asynchronous island evolution
		arch.evolve(m_gen);
		arch.join();
	} else {
		for(int g = 0; g < m_gen; ++g) { //batched island evolution
			arch.evolve_batch(1, m_threads);
		}
	}

	// Finally, we assemble the evolved population selecting from the original one + the evolved one
	// the best NP (crowding distance)
	population popnew(pop);
	for(pagmo::population::size_type i=0; i<arch.get_size() ;++i) {
		popnew.push_back(arch.get_island(i)->get_population().champion().x);
	}
	std::vector<population::size_type> selected_idx = popnew.get_best_idx(NP);
	// We completely clear the population (NOTE: memory of all individuals and the notion of
	// champion is thus destroyed)
	pop.clear();
	// And we recreate it with the best NP among the evolved and the new population
	for (population::size_type i=0; i < NP; ++i) {
		pop.push_back(popnew.get_individual(selected_idx[i]).cur_x);
	}
}

/// Algorithm name
std::string pade::get_name() const
{
	return "Parallel Decomposition (PaDe)";
}

/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string pade::human_readable_extra() const
{
	std::ostringstream s;
	s << "gen:" << m_gen << ' ';
	s << "threads:" << m_threads << ' ';
	s << "solver:" << m_solver->get_name() << ' ';
	s << "neighbours:" << m_T << ' ';
	s << "decomposition:";
	switch (m_method)
	{
		case pagmo::problem::decompose::BI : s << "BI" << ' ';
			break;
		case pagmo::problem::decompose::WEIGHTED : s << "WEIGHTED" << ' ';
			break;
		case pagmo::problem::decompose::TCHEBYCHEFF : s << "TCHEBYCHEFF" << ' ';
			break;
	}
	s << "weights:";
	switch (m_weight_generation)
	{
		case RANDOM : s << "RANDOM" << ' ';
			break;
		case LOW_DISCREPANCY : s << "LOW_DISCREPANCY" << ' ';
			break;
		case GRID : s << "GRID" << ' ';
			break;
	}
	s << "ref. point" << m_z;
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::pade)
