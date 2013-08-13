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
 * @param[in] minW minimum particles' inertia weight (the inertia weight is decreased troughout the run between maxW and minW)
 * @param[in] maxW maximum particles' inertia weight (the inertia weight is decreased troughout the run between maxW and minW)
 * @param[in] C1 magnitude of the force, applied to the particle's velocity, in the direction of its previous best position
 * @param[in] C2 magnitude of the force, applied to the particle's velocity, in the direction of its global best (leader)
 * @param[in] CHI velocity scaling factor
 * @param[in] m_v_coeff velocity coefficient (determining the maximum allowed particle velocity)
 * @param[in] leader_selection_range the leader of each particle is selected among the best leader_selection_range% individuals
 *
 * @throws value_error if gen is negative
 */
mopso::mopso(int gen, double minW, double maxW, double C1, double C2,
	  double CHI, double v_coeff, int leader_selection_range):base(),
	m_gen(gen),
	m_minW(minW),
	m_maxW(maxW),
	m_C1(C1),
	m_C2(C2),
	m_CHI(CHI),
	m_v_coeff(v_coeff),
	m_leader_selection_range(leader_selection_range)
{
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}

	if (minW <= 0 || maxW <= 0 || C1 <=0 || C2 <= 0 || CHI <=0) {
		pagmo_throw(value_error,"minW, maxW, C1, C2 and CHI should be greater than 0");
	}

	if (minW > maxW) {
		pagmo_throw(value_error,"minW should be lower than maxW");
	}

	if (v_coeff <= 0 || v_coeff > 1) {
		pagmo_throw(value_error,"v_coeff should be in the ]0,1] range");
	}

	if (leader_selection_range <=0 || leader_selection_range > 100) {
		pagmo_throw(value_error,"leader_selection_range should be in the ]0,100] range");
	}
}

/// Copy constructor. Performs a deep copy. Necessary as a pointer to a base algorithm is here contained
mopso::mopso(const mopso &algo):base(algo), m_gen(algo.m_gen),m_minW(algo.m_minW), m_maxW(algo.m_maxW),m_C1(algo.m_C1),
	m_C2(algo.m_C2), m_CHI(algo.m_CHI), m_v_coeff(algo.m_v_coeff), m_leader_selection_range(algo.m_leader_selection_range)
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
	const problem::base             &prob = pop.problem();
	const problem::base::size_type   D = prob.get_dimension(), prob_i_dimension = prob.get_i_dimension(), prob_c_dimension = prob.get_c_dimension(), prob_f_dimension = prob.get_f_dimension();
	const problem::base::size_type   Dc = D - prob_i_dimension;
	const decision_vector           &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type      NP = pop.size();

	//We perform some checks to determine wether the problem/population are suitable for PSO
	if( Dc == 0 ){
		pagmo_throw(value_error,"There is no continuous part in the problem decision vector for MOPSO to optimise");
	}

	if( prob_c_dimension != 0 ){
		pagmo_throw(value_error,"The problem is not box constrained and MOPSO is not suitable to solve it");
	}

	if( prob_f_dimension < 2 ){
		pagmo_throw(value_error,"The problem is not multi-objective. Use a single-objectie optimization algorithm instead");
	}

	// Get out if there is nothing to do.
	if (NP == 0 || m_gen == 0) {
		return;
	}

	decision_vector minV(Dc), maxV(Dc);	// Maximum and minimum velocity allowed

	// Initialise the minimum and maximum velocity
	for(problem::base::size_type d = 0; d < Dc; d++ ){
		double v_width  = (ub[d] - lb[d]) * m_v_coeff;
		minV[d] = -1.0 * v_width;
		maxV[d] = v_width;
	}

	for(int g = 0; g < m_gen; ++g) {
		std::cout << "gen: " << g << std::endl;

		population pop_best_x(pop);
		pop_best_x.clear();
		for(population::size_type idx = 0; idx < NP; ++idx) {
			pop_best_x.push_back(pop.get_individual(idx).best_x);
			pop_best_x.set_v(idx, pop.get_individual(idx).cur_v);
		}

		//Calculate the best leader_selection_range% individuals of the original population according to the crowing distance
		std::vector<population::size_type> bestIndividuals = pop.get_best_idx((int)ceil(NP*m_leader_selection_range/100.0));

		const double W  = m_maxW - (m_maxW-m_minW)/m_gen * g; //W decreased from maxW to minW troughout the run

		for(population::size_type idx = 0; idx < NP; ++idx) {

			//Set as a leader for the current particle a random particle among the best ones in the leader_selection_range percentile
			const decision_vector &leader = pop.get_individual(
						bestIndividuals[boost::uniform_int<int>(0,bestIndividuals.size()-1)(m_drng)]).cur_x;

			//Calculate some random factors
			const double r1 = boost::uniform_real<double>(0,1)(m_drng);
			const double r2 = boost::uniform_real<double>(0,1)(m_drng);

			const decision_vector &best_X = pop.get_individual(idx).best_x;
			const decision_vector &cur_X = pop.get_individual(idx).cur_x;
			const decision_vector &cur_V = pop.get_individual(idx).cur_v;

			//Calculate new velocity and new position for each particle
			decision_vector newX(Dc);
			decision_vector newV(Dc);
			for(decision_vector::size_type i = 0; i < cur_X.size(); ++i) {
				double v = W*cur_V[i] +
								m_C1*r1*(best_X[i] - cur_X[i]) +
								m_C2*r2*(leader[i] - cur_X[i]);

				if(v > maxV[i]){
					v = maxV[i];
				}
				else if(v < minV[i]) {
					v = minV[i];
				}

				double x = cur_X[i] + m_CHI*v;
				if(x > ub[i]) {
					x = ub[i];
					v = 0;
				} else if (x < lb[i]) {
					x = lb[i];
					v = 0;
				}

				newV[i] = v;
				newX[i] = x;
			}

			pop.set_x(idx, newX); 
			pop.set_v(idx, newV);
		}


		//The first NP individuals contain the actual population (that will be mutated) the last NP contains
		// individuals having as X the best_x of each individual of the current population
		population newPop = population(pop);
		for(population::size_type idx = 0; idx < NP; ++idx) {
			newPop.push_back(pop_best_x.get_individual(idx).cur_x);
			newPop.set_v(idx+NP,pop_best_x.get_individual(idx).cur_v);
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
	s << "minW:" << m_minW << ' ';
	s << "maxW:" << m_maxW << ' ';
	s << "C1:" << m_C1 << ' ';
	s << "C2:" << m_C2 << ' ';
	s << "CHI:" << m_CHI << ' ';
	s << "v_coeff:" << m_v_coeff << ' ';
	s << "leader_selection_range:" << m_leader_selection_range << ' ';
	return s.str();
}


}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::mopso);
