/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
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


#include <string>
#include <vector>
#include <algorithm>

#include <boost/math/special_functions/binomial.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../archipelago.h"
#include "../island.h"
#include "../population.h"
#include "../topology/fully_connected.h"
#include "../topology/custom.h"
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
 * @param[in] max_parallelism limits the amounts of threads spawned
 * @param[in] method the decomposition method to use (Weighted, Tchebycheff or BI)
 * @param[in] solver the algorithm to solve the single objective problems.
 * @param[in] T the size of the population on each subproblem (must be an even number)
 * @param[in] weight_generation the method to generate the weight vectors (RANDOM, GRID or LOW-DISCREPANCY)
 * @param[in] z the reference point used for decomposition (with Tchebycheff and BI)
 *
 * @throws value_error if gen is negative, weight_generation is not sane
 * @see pagmo::problem::decompose::method_type
 */
pade::pade(int gen, unsigned int max_parallelism, pagmo::problem::decompose::method_type method,
		   const pagmo::algorithm::base & solver, population::size_type T, weight_generation_type weight_generation,
		   const fitness_vector &z)
	:base(),m_gen(gen),m_max_parallelism(max_parallelism),
	  m_method(method),m_solver(solver.clone()),m_T(T),m_weight_generation(weight_generation),m_z(z)
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
pade::pade(const pade &algo):base(algo), m_gen(algo.m_gen),m_max_parallelism(algo.m_max_parallelism),
	m_method(algo.m_method),m_solver(algo.m_solver->clone()),m_T(algo.m_T),
	m_weight_generation(algo.m_weight_generation),m_z(algo.m_z)
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

	// Definition of useful probability distributions
	boost::uniform_real<double> uniform(0.0,1.0);
	boost::variate_generator<boost::lagged_fibonacci607 &, boost::uniform_real<double> > r_dist(m_drng,uniform);

	decision_vector dummy(pop.get_individual(0).cur_x.size(),0); //used for initialisation purposes
	std::vector<decision_vector> X(NP,dummy); //set of population chromosomes

	// Copy the population chromosomes into X
	for ( population::size_type i = 0; i<NP; i++ ) {
		X[i]	=	pop.get_individual(i).cur_x;
	}
	// Clear the current population (TODO: is this necessary?)
	pop.clear();

	// Generate the weights for the decomposed problems
	std::vector<fitness_vector> weights;

	if(m_weight_generation == GRID) {
		//find the largest H resulting in a population smaller or equal to NP
		unsigned int H;
		if (prob.get_f_dimension() == 2) {
			H = NP-1;
		} else if (prob.get_f_dimension() == 3) {
			H = floor(0.5 * (sqrt(8*NP + 1) - 3));
		} else {
			std::cout << "Fitness dimension is " << prob.get_f_dimension() << std::endl;
			H = 1;
			while(boost::math::binomial_coefficient<double>(H+prob.get_f_dimension()-1, prob.get_f_dimension()-1) <= NP) {
				++H;
			}
			H--;
		}

		// We check that NP equals the population size rsulting from H
		if (fabs(NP-(boost::math::binomial_coefficient<double>(H+prob.get_f_dimension()-1, prob.get_f_dimension()-1))) > 1E-8) {
			std::ostringstream error_message;
			error_message << "Invalid population size. Select " << boost::math::binomial_coefficient<double>(H+prob.get_f_dimension()-1, prob.get_f_dimension()-1)
					<< " or " << boost::math::binomial_coefficient<double>(H+1+prob.get_f_dimension()-1, prob.get_f_dimension()-1)
					<< ".";
			pagmo_throw(value_error,error_message.str());
		}

		// We generate the weights
		std::vector<unsigned int> range;
		for (unsigned int i=0; i<H+1;++i) {
			range.push_back(i);
		}
		double epsilon = 1E-6;
		reksum(weights, range, prob.get_f_dimension(), H);
		for(unsigned int i=0; i< weights.size(); ++i) {
			for(unsigned int j=0; j< weights[i].size(); ++j) {
				weights[i][j] += epsilon;  //to avoid to have any weight equal to zero
				weights[i][j] /= H+epsilon*weights[i].size();
			}
		}

	} else if(m_weight_generation == LOW_DISCREPANCY) {
		pagmo::util::discrepancy::simplex generator(prob.get_f_dimension(),1);
		for(unsigned int i = 0; i <NP; ++i) {
			weights.push_back(generator());
		}

	} else if(m_weight_generation == RANDOM) {
		pagmo::util::discrepancy::project_2_simplex projection(prob.get_f_dimension());
		for (unsigned int i = 0; i<NP; ++i) {
			fitness_vector dummy(prob.get_f_dimension()-1,0.0);
			for(unsigned int j = 0; j <prob.get_f_dimension()-1; ++j) {
				dummy[j] = r_dist();
			}
			weights.push_back(projection(dummy));
		}
	}


	// We then compute, for each weight vector, which ones are the closest ones (this will form the topology later on)
	std::vector<std::vector<int> > indices;
	pagmo::util::neighbourhood::euclidian::compute_neighbours(indices, weights);

	// Create the archipelago of NP islands:
	// each island in the archipelago solves a different single-objective problem.
	// We use here the broadcast migration model. This will force, at each migration,
	// to have individuals from all connected island to be inserted.
	pagmo::archipelago arch(pagmo::archipelago::broadcast);
	
	// Best individual will be selected for migration
	const pagmo::migration::best_s_policy  selection_policy;
	
	// As m_T neighbours are connected, we replace m_T individuals on the island
	const pagmo::migration::worst_r_policy replacement_policy(m_T);

	for(pagmo::population::size_type i=0; i<NP;++i) {
		pagmo::problem::decompose decomposed_prob(prob, m_method,weights[i],m_z);
		pagmo::population decomposed_pop(decomposed_prob);

		//Set the individuals of the new population as one individual of the original population plus m_T
		//neighbours individuals
		if(m_T < NP-1) {
			for(pagmo::population::size_type  j = 0; j <= m_T; ++j) { 
				decomposed_pop.push_back(X[indices[i][j]]);
			}
		} else { //complete topology
			for(pagmo::population::size_type  j = 0 ; j < NP; ++j) {
				decomposed_pop.push_back(X[j]);
			}
		}
		arch.push_back(pagmo::island(*m_solver,decomposed_pop, 1.0, selection_policy, replacement_policy));
	}

	if(m_T >= NP-1) {
		arch.set_topology(topology::fully_connected());
	} else {
		topology::custom topo;
		for(unsigned int i = 0; i < NP; ++i) {
			topo.push_back();
		}
		for(unsigned int i = 0; i < NP; ++i) { //connect each island with the T closest neighbours
			for(unsigned int j = 1; j <= m_T; ++j) { //start from 1 to avoid to connect with itself
				topo.add_edge(i,indices[i][j]);
			}
		}
		arch.set_topology(topo);
	}

	//Evolve the archipelago for m_gen generations
	if(m_max_parallelism == NP) { //asynchronous island evolution
		arch.evolve(m_gen);
		arch.join();
	} else {
		for(int g = 0; g < m_gen; ++g) { //batched island evolution
			arch.evolve_batch(1, m_max_parallelism);
		}
	}

	//The original population is set to contain the best individual of each island
	for(pagmo::population::size_type i=0; i<NP ;++i) {
		pop.push_back(arch.get_island(i)->get_population().champion().x);
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
	s << "max_parallelism:" << m_max_parallelism << ' ';
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
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::pade);
