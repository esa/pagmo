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

#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/special_functions/round.hpp>


#include <string>
#include <vector>
#include <algorithm>

#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../types.h"
#include "base.h"
#include "nsga2.h"

namespace pagmo { namespace algorithm {

/// Constructor
 /**
 * Constructs a NSGA II algorithm
 *
 * @param[in] gen Number of generations to evolve.
 * @param[in] cr Crossover probability
 * @param[in] eta_c Distribution index for crossover
 * @param[in] m Mutation probability
 * @param[in] eta_m Distribution index for mutation
 * @throws value_error if gen is negative, crossover probability is not \f$ \in [0,1[\f$, mutation probability or mutation width is not \f$ \in [0,1]\f$,
 */
nsga2::nsga2(int gen, double cr, double eta_c, double m, double eta_m):base(),m_gen(gen),m_cr(cr),m_eta_c(eta_c),m_m(m),m_eta_m(eta_m)
{
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if (cr >= 1 || cr < 0) {
		pagmo_throw(value_error,"crossover probability must be in the [0,1[ range");
	}
	if (m < 0 || m > 1) {
		pagmo_throw(value_error,"mutation probability must be in the [0,1] range");
	}
	if (eta_c <1 || eta_c >= 100) {
		pagmo_throw(value_error,"Distribution index for crossover must be in 1..100");
	}
	if (eta_m <1 || eta_m >= 100) {
		pagmo_throw(value_error,"Distribution index for mutation must be in 1..100");
	}
}

/// Clone method.
base_ptr nsga2::clone() const
{
	return base_ptr(new nsga2(*this));
}

pagmo::population::size_type nsga2::tournament_selection(pagmo::population::size_type idx1, pagmo::population::size_type idx2, const pagmo::population& pop) const
{
	if (pop.get_pareto_rank(idx1) < pop.get_pareto_rank(idx2)) return idx1;
	if (pop.get_pareto_rank(idx1) > pop.get_pareto_rank(idx2)) return idx2;
	if (pop.get_crowding_d(idx1) > pop.get_crowding_d(idx2)) return idx1;
	if (pop.get_crowding_d(idx1) < pop.get_crowding_d(idx2)) return idx2;
	return ((m_drng() > 0.5) ? idx1 : idx2);
}

void nsga2::crossover(decision_vector& child1, decision_vector& child2, pagmo::population::size_type parent1_idx, pagmo::population::size_type parent2_idx,const pagmo::population& pop) const
{

		problem::base::size_type D = pop.problem().get_dimension();
		problem::base::size_type Di = pop.problem().get_i_dimension();
		problem::base::size_type Dc = D - Di;
	const decision_vector &lb = pop.problem().get_lb(), &ub = pop.problem().get_ub();
	const decision_vector& parent1 = pop.get_individual(parent1_idx).cur_x;
	const decision_vector& parent2 = pop.get_individual(parent2_idx).cur_x;
	double y1,y2,yl,yu, rand, beta, alpha, betaq, c1, c2;
	child1 = parent1;
	child2 = parent2;
		int site1, site2;

		//This implements a Simulated Binary Crossover SBX
	if (m_drng() <= m_cr) {
		for (pagmo::problem::base::size_type i = 0; i < Dc; i++) {
			if ( (m_drng() <= 0.5) && (std::fabs(parent1[i]-parent2[i]) ) > 1.0e-14) {
				if (parent1[i] < parent2[i]) {
					y1 = parent1[i];
					y2 = parent2[i];
				} else {
					y1 = parent2[i];
					y2 = parent1[i];
				}
				yl = lb[i];
				yu = ub[i];
				rand = m_drng();

				beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
				alpha = 2.0 - std::pow(beta,-(m_eta_c+1.0));
				if (rand <= (1.0/alpha))
				{
					betaq = std::pow((rand*alpha),(1.0/(m_eta_c+1.0)));
				} else {
					betaq = std::pow((1.0/(2.0 - rand*alpha)),(1.0/(m_eta_c+1.0)));
				}
				c1 = 0.5*((y1+y2)-betaq*(y2-y1));

				beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
				alpha = 2.0 - std::pow(beta,-(m_eta_c+1.0));
				if (rand <= (1.0/alpha))
				{
					betaq = std::pow((rand*alpha),(1.0/(m_eta_c+1.0)));
				} else {
					betaq = std::pow((1.0/(2.0 - rand*alpha)),(1.0/(m_eta_c+1.0)));
				}
				c2 = 0.5*((y1+y2)+betaq*(y2-y1));

				if (c1<lb[i]) c1=lb[i];
				if (c2<lb[i]) c2=lb[i];
				if (c1>ub[i]) c1=ub[i];
				if (c2>ub[i]) c2=ub[i];
				if (m_drng() <= 0.5) {
					child1[i] = c1; child2[i] = c2;
				} else {
					child1[i] = c2; child2[i] = c1;
				}
			}
		}

			}

		//This implements two point binary crossover
		for (pagmo::problem::base::size_type i = Dc; i < D; i++) {
			   if (m_drng() <= m_cr) {
					 boost::uniform_int<int> in_dist(0,Di-1);
					 boost::variate_generator<boost::mt19937 &, boost::uniform_int<int> > ra_num(m_urng,in_dist);
					 site1 = ra_num();
					 site2 = ra_num();
					 if (site1 > site2) std::swap(site1,site2);
					 for(int j=0; j<site1; j++)
					 {
						 child1[j] = parent1[j];
						 child2[j] = parent2[j];
					 }
					 for(int j=site1; j<site2; j++)
					 {
						 child1[j] = parent2[j];
						 child2[j] = parent1[j];
					 }
					 for(pagmo::problem::base::size_type j=site2; j<Di; j++)
					 {
						 child1[j] = parent1[j];
						 child2[j] = parent2[j];
					 }
			  }
			  else {
				   child1[i] = parent1[i];
				   child2[i] = parent2[i];
			  }
	}
}

void nsga2::mutate(decision_vector& child, const pagmo::population& pop) const
{

	problem::base::size_type D = pop.problem().get_dimension();
		problem::base::size_type Di = pop.problem().get_i_dimension();
		problem::base::size_type Dc = D - Di;
	const decision_vector &lb = pop.problem().get_lb(), &ub = pop.problem().get_ub();
	double rnd, delta1, delta2, mut_pow, deltaq;
	double y, yl, yu, val, xy;
		int gen_num;

		//This implements the real polinomial mutation of an individual
	for (pagmo::problem::base::size_type j=0; j < Dc; ++j){
		if (m_drng() <= m_m) {
			y = child[j];
			yl = lb[j];
			yu = ub[j];
			delta1 = (y-yl)/(yu-yl);
			delta2 = (yu-y)/(yu-yl);
			rnd = m_drng();
			mut_pow = 1.0/(m_eta_m+1.0);
			if (rnd <= 0.5)
			{
				xy = 1.0-delta1;
				val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(m_eta_m+1.0)));
				deltaq =  pow(val,mut_pow) - 1.0;
			}
			else
			{
				xy = 1.0-delta2;
				val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(m_eta_m+1.0)));
				deltaq = 1.0 - (pow(val,mut_pow));
			}
			y = y + deltaq*(yu-yl);
			if (y<yl) y = yl;
			if (y>yu) y = yu;
			child[j] = y;
		}
	}

		//This implements the integer mutation for an individual
		 for (pagmo::problem::base::size_type j=Dc; j < D; ++j){
			   if (m_drng() <= m_m) {
						y = child[j];
						yl = lb[j];
						yu = ub[j];
						boost::uniform_int<int> in_dist(yl,yu-1);
						boost::variate_generator<boost::mt19937 &, boost::uniform_int<int> > ra_num(m_urng,in_dist);
						gen_num = ra_num();
						if (gen_num >= y) gen_num = gen_num + 1;
						child[j] = gen_num;
		 }
	  }
}
/// Evolve implementation.
/**
 * Run the NSGA-II algorithm for the number of generations specified in the constructors.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void nsga2::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type D = prob.get_dimension();
	const problem::base::size_type prob_c_dimension = prob.get_c_dimension();
	const population::size_type NP = pop.size();


	//We perform some checks to determine wether the problem/population are suitable for NSGA-II
	if ( prob_c_dimension != 0 ) {
			pagmo_throw(value_error, "The problem is not box constrained and NSGA-II is not suitable to solve it");
	}

	if (NP < 5 || (NP % 4 != 0) ) {
		pagmo_throw(value_error, "for NSGA-II at least 5 individuals in the population are needed and the population size must be a multiple of 4");
	}

	if ( prob.get_f_dimension() < 2 ) {
		pagmo_throw(value_error, "The problem is not multiobjective, try some other algorithm than NSGA-II");
	}

	// Get out if there is nothing to do.
	if (m_gen == 0) {
		return;
	}

	std::vector<population::size_type> best_idx(NP), shuffle1(NP),shuffle2(NP);
	population::size_type parent1_idx, parent2_idx;
	decision_vector child1(D), child2(D);

	for (pagmo::population::size_type i=0; i< NP; i++) shuffle1[i] = i;
	for (pagmo::population::size_type i=0; i< NP; i++) shuffle2[i] = i;

	boost::uniform_int<int> pop_idx(0,NP-1);
	boost::variate_generator<boost::mt19937 &, boost::uniform_int<int> > p_idx(m_urng,pop_idx);

	// Main NSGA-II loop
	for (int g = 0; g<m_gen; g++) {
		//At each generation we make a copy of the population into popnew
		// We compute the crowding distance and the pareto rank of pop
		pop.update_pareto_information();
		population popnew(pop);

		//We create some pseudo-random permutation of the poulation indexes
		std::random_shuffle(shuffle1.begin(),shuffle1.end(),p_idx);
		std::random_shuffle(shuffle2.begin(),shuffle2.end(),p_idx);

		//We then loop thorugh all individuals with increment 4 to select two pairs of parents that will
		//each create 2 new offspring
		for (pagmo::population::size_type i=0; i< NP; i+=4) {
			// We create two offsprings using the shuffled list 1
			parent1_idx = tournament_selection(shuffle1[i], shuffle1[i+1],pop);
			parent2_idx = tournament_selection(shuffle1[i+2], shuffle1[i+3],pop);
			crossover(child1, child2, parent1_idx,parent2_idx,pop);
			mutate(child1,pop);
			mutate(child2,pop);
			popnew.push_back(child1);
			popnew.push_back(child2);

			// We repeat with the shuffled list 2
			parent1_idx = tournament_selection(shuffle2[i], shuffle2[i+1],pop);
			parent2_idx = tournament_selection(shuffle2[i+2], shuffle2[i+3],pop);
			crossover(child1, child2, parent1_idx,parent2_idx,pop);
			mutate(child1,pop);
			mutate(child2,pop);
			popnew.push_back(child1);
			popnew.push_back(child2);
		} // popnew now contains 2NP individuals

		// This method returns the sorted N best individuals in the population according to the crowded comparison operator
		// defined in population.cpp
		best_idx = popnew.get_best_idx(NP);
		// We completely cancel the population (NOTE: memory of all individuals and the notion of
		// champion is thus destroyed)
		pop.clear();
		for (population::size_type i=0; i < NP; ++i) pop.push_back(popnew.get_individual(best_idx[i]).cur_x);
	} // end of main SGA loop
}

/// Algorithm name
std::string nsga2::get_name() const
{
	return "Nondominated Sorting Genetic Algorithm II (NSGA-II)";
}

/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string nsga2::human_readable_extra() const
{
	std::ostringstream s;
	s << "gen:" << m_gen << ' ';
	s << "cr:" << m_cr << ' ';
	s << "eta_c:" << m_eta_c << ' ';
	s << "m:" << m_m << ' ';
	s << "eta_m:" << m_eta_m << std::endl;

	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::nsga2)
