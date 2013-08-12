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
#include "../problem/decompose.h"
#include "../util/discrepancy.h"
#include "../util/neighbourhood.h"
#include "../migration/worst_r_policy.h"
#include "../migration/best_s_policy.h"
#include "../types.h"
#include "base.h"
#include "mopso.h"

namespace pagmo { namespace algorithm {
/// Constructor
 /**
 * Constructs a mopso algorithm
 *
 * @param[in] gen Number of generations to evolve.
 *
 * @throws value_error if gen is negative
 */
mopso::mopso(int gen)
	:base(),m_gen(gen)
{
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
}

/// Copy constructor. Performs a deep copy. Necessary as a pointer to a base algorithm is here contained
mopso::mopso(const mopso &algo):base(algo), m_gen(algo.m_gen)
{}

/// Clone method.
base_ptr mopso::clone() const
{
	return base_ptr(new mopso(*this));
}

/// Evolve implementation.
/**
 * Run the MOPSO algorithm for the number of generations specified in the constructors.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */
void mopso::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const population::size_type NP = pop.size();

	if ( prob.get_f_dimension() < 2 ) {
		pagmo_throw(value_error, "The problem is not multiobjective, try some other algorithm than mopso");
	}

	// Get out if there is nothing to do.
	if (m_gen == 0) {
		return;
	}

	for(int g = 0; g < m_gen; ++g) {
		std::cout<<"gen: " << g << std::endl;

		population pop_best_x(pop);
		pop_best_x.clear();
		for(population::size_type idx = 0; idx < NP; ++idx) {
			pop_best_x.push_back(pop.get_individual(idx).best_x);
		}

		//Calculate the best 5% individuals of the original population according to the crowing distance
		std::vector<population::size_type> bestIndividuals = pop.get_best_idx((int)ceil(NP/20)); //best 5%
		for(population::size_type idx = 0; idx < NP; ++idx) {

			//Set as a leader for the current particle a random particle among the 5% best ones
			decision_vector bestX = pop.get_individual(
						bestIndividuals[boost::uniform_int<int>(0,bestIndividuals.size()-1)(m_drng)]).cur_x;

			//Calculate some random factors
			const double W  = 1-0.6/m_gen*g;//boost::uniform_real<double>(0.1,0.5)(m_drng);
			const double C1 = 2.0;//boost::uniform_real<double>(1.5,2)(m_drng);
			const double C2 = 2.0;//boost::uniform_real<double>(1.5,2)(m_drng);
			const double r1 = boost::uniform_real<double>(0,1)(m_drng);
			const double r2 = boost::uniform_real<double>(0,1)(m_drng);
			const double chi = 1.0;
			decision_vector maxV = pop.problem().get_ub();

			//Calculate new velocity and new position for each particle
			decision_vector newX;
			decision_vector newV;
			for(decision_vector::size_type i = 0; i < pop.get_individual(idx).cur_x.size(); ++i) {
				double v = W*pop.get_individual(idx).cur_v[i] +
								C1*r1*(pop.get_individual(idx).best_x[i] - pop.get_individual(idx).cur_x[i]) +
								C2*r2*(bestX[i] - pop.get_individual(idx).cur_x[i]);
				double x = pop.get_individual(idx).cur_x[i] + chi*v;
				if(v>maxV.at(i)){
					v = maxV.at(i);
				}
				if(x > pop.problem().get_ub()[i]) {
					x = pop.problem().get_ub()[i];
					v = 0;
				} else if (x < pop.problem().get_lb()[i]) {
					x = pop.problem().get_lb()[i];
					v = 0;
				}
				newV.push_back(v);
				newX.push_back(x);
			}

			pop.set_x(idx, newX); 
			pop.set_v(idx, newV);
		}


		//The first NP individuals contain the actual population (that will be mutated) the last NP contains
		// individuals having as X the best_x of each individual of the current population
		population newPop = population(pop);
		for(population::size_type idx = 0; idx < NP; ++idx) {
			newPop.push_back(pop_best_x.get_individual(idx).cur_x);
		}

		//Select the best NP individuals in the new population (of size 2*NP) according to the crowding distance
		std::vector<population::size_type> bestIndices = newPop.get_best_idx(NP);

		//Set the population accordingly
		pop.clear();
		for(population::size_type i = 0; i < NP; ++i) {
			pop.push_back(newPop.get_individual(bestIndices[i]).cur_x);
			pop.set_x(i, newPop.get_individual(bestIndices[i]).best_x);
			pop.set_x(i, newPop.get_individual(bestIndices[i]).cur_x); 
			pop.set_v(i, newPop.get_individual(bestIndices[i]).cur_v);
		}

	}

}

/// Algorithm name
std::string mopso::get_name() const
{
	return "Multi Objective Particle Swarm Optimizer (MOPSO)";
}

/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string mopso::human_readable_extra() const
{
	std::ostringstream s;
	s << "gen:" << m_gen << ' ';
	return s.str();
}


}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::mopso);
