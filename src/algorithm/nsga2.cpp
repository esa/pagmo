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

//Used in std::sort to compare individuals according to a specific component of the fitness
class CompareFitness: std::binary_function<std::pair<decision_vector,int> , std::pair<decision_vector,int>, bool>
{
	const pagmo::problem::base *prob;
	int n;
	
	public: 
		/**
		 * @param[in] p a reference to the problem
		 * @param[in] component the component according to which we want to compare the fitness
		 */
		CompareFitness(const pagmo::problem::base &p, int component) {
			prob = &p;
			n = component;
		}

		bool operator()(const std::pair<decision_vector,int> &a, const std::pair<decision_vector,int> &b) const {
			return (prob->objfun(a.first) < prob->objfun(b.first));
		}
};

//Used in std::sort to compare individuals according to their crowding distance value
class CompareDistance: std::binary_function<int , int, bool>
{
	const std::vector<int> *I;
	
	public: 
		/**
		* @param[in] distances a reference to the crowding distances vector
		*/
		CompareDistance(const std::vector<int> &distances) {
			I = &distances;
		}

		bool operator()(int a, int b) const {
			return (I->at(a) > I->at(b));
		}
};
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] gen Number of generations to evolve.
 * @param[in] cr Crossover probability (of each allele if binomial crossover)
 * @param[in] m Mutation probability (of each allele)
 * @param[in] mut Mutation type. One of nsga2::mutation::GAUSSIAN, nsga2::mutation::RANDOM
 * @param[in] width Mutation width. When gaussian mutation is selected is the width of the mutation
 * @param[in] cro Crossover type. One of nsga2::crossover::BINOMIAL, nsga2::crossover::EXPONENTIAL
 * @throws value_error if gen is negative, crossover probability is not \f$ \in [0,1]\f$, mutation probability or mutation width is not \f$ \in [0,1]\f$,

 */
nsga2::nsga2(int gen, const double &cr, const double &m, mutation::type mut, double width, crossover::type cro)
	:base(),m_gen(gen),m_cr(cr),m_m(m),m_mut(mut,width),m_cro(cro)
{
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if (cr > 1 || cr < 0) {
		pagmo_throw(value_error,"crossover probability must be in the [0,1] range");
	}
	if (m < 0 || m > 1) {
		pagmo_throw(value_error,"mutation probability must be in the [0,1] range");
	}
	if (width <0 || width >1) {
		pagmo_throw(value_error,"mutation width must be in the [0,1] range");
	}

}

/// Clone method.
base_ptr nsga2::clone() const
{
	return base_ptr(new nsga2(*this));
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
	const problem::base::size_type D = prob.get_dimension(), Di = prob.get_i_dimension(), prob_c_dimension = prob.get_c_dimension(), prob_f_dimension = prob.get_f_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP = pop.size();
	const problem::base::size_type Dc = D - Di;


	//We perform some checks to determine wether the problem/population are suitable for NSGA-II
	if ( prob_c_dimension != 0 ) {
		pagmo_throw(value_error,"The problem is not box constrained and SGA is not suitable to solve it");
	}

	if (NP < 5) {
		pagmo_throw(value_error,"for NSGA-II at least 5 individuals in the population are needed");
	}

	// Get out if there is nothing to do.
	if (m_gen == 0) {
		return;
	}
	// Some vectors used during evolution are allocated here.
	decision_vector dummy(D,0);			//used for initialisation purposes
	std::vector<decision_vector > X(NP,dummy), Xnew(NP,dummy), P(2*NP, dummy);

	// Initialise the chromosomes and the Xnew vector.
	for (pagmo::population::size_type i = 0; i<NP; i++ ) {
		X[i]	=	pop.get_individual(i).cur_x;
		Xnew[i]	=	pop.get_individual(i).cur_x;
	}

	std::vector<std::vector<int> > S(2*NP); //Domination sets
	std::vector<std::vector<int> > F(1); //Domination fronts
	std::vector<int> n(2*NP,0); //Domination counters

	// Main SGA loop
	for (int j = 0; j<m_gen; j++) {
		
		//Crossover
		int r1,L;
		decision_vector  member1,member2;

		for (pagmo::population::size_type i=0; i< NP; i++) {
			//for each chromosome selected i.e. in Xnew
			member1 = Xnew[i];
			do {
			//we select a mating patner different from the self (i.e. no masturbation allowed)
				r1 = boost::uniform_int<int>(0,NP - 1)(m_urng);
			} while ( r1 == boost::numeric_cast<int>(i) );
			member2 = Xnew[r1];
			//and we operate crossover
			switch (m_cro) {
				//0 - binomial crossover
			case crossover::BINOMIAL: {
				size_t n = boost::uniform_int<int>(0,D-1)(m_urng);
				for (size_t L = 0; L < D; ++L) { /* perform D binomial trials */
					if ((m_drng() < m_cr) || L + 1 == D) { /* change at least one parameter */
						member1[n] = member2[n];
					}
					n = (n+1)%D;
				}
				break; }
				//1 - exponential crossover
			case crossover::EXPONENTIAL: {
				size_t n = boost::uniform_int<int>(0,D-1)(m_urng);
				L = 0;
				do {
					member1[n] = member2[n];
					n = (n+1) % D;
					L++;
				}  while ( (m_drng() < m_cr) && (L < boost::numeric_cast<int>(D)) );
				break; }
			}
			Xnew[i] = member1;

		} 

		//Mutation
		switch (m_mut.m_type) {
			case mutation::GAUSSIAN: {
				boost::normal_distribution<double> dist;
				boost::variate_generator<boost::lagged_fibonacci607 &, boost::normal_distribution<double> > delta(m_drng,dist);
				for (pagmo::problem::base::size_type k = 0; k < D;k++) { //for each continuous variable
					double std = (ub[k]-lb[k]) * m_mut.m_width;
					for (pagmo::population::size_type i = 0; i < NP;i++) { //for each individual
						if (m_drng() < m_m) {
							double mean = Xnew[i][k];
							Xnew[i][k] = (delta() * std + mean);
							if (Xnew[i][k] > ub[k]) Xnew[i][k] = ub[k];
							if (Xnew[i][k] < lb[k]) Xnew[i][k] = lb[k];
						}
					}
				}
				for (pagmo::problem::base::size_type k = Dc; k < D;k++) { //for each integer variable
					for (pagmo::population::size_type i = 0; i < NP;i++) { //for each individual
						if (m_drng() < m_m) {
							double mean = Xnew[i][k];
							Xnew[i][k] = boost::math::iround(delta() + mean);
							if (Xnew[i][k] > ub[k]) Xnew[i][k] = ub[k];
							if (Xnew[i][k] < lb[k]) Xnew[i][k] = lb[k];
						}
					}
				}
				break;
				}
			case mutation::RANDOM: {
				for (pagmo::population::size_type i = 0; i < NP;i++) {
					for (pagmo::problem::base::size_type j = 0; j < Dc;j++) { //for each continuous variable
						if (m_drng() < m_m) {
							Xnew[i][j] = boost::uniform_real<double>(lb[j],ub[j])(m_drng);
						}
					}
					for (pagmo::problem::base::size_type j = Dc; j < D;j++) {//for each integer variable
						if (m_drng() < m_m) {
							Xnew[i][j] = boost::uniform_int<int>(lb[j],ub[j])(m_urng);
						}
					}
				}
				break;
				}
		}

		//A new bigger population defined as the old one plus the crossovered/mutated one
		for(pagmo::population::size_type i = 0; i < NP; ++i) {
			P[i] = X[i];
		}
		for(pagmo::population::size_type i = NP; i < 2*NP; ++i) {
			P[i] = Xnew[i-NP];
		}

		//Non-dominated sort on population P
		for(pagmo::population::size_type p = 0; p < 2*NP; ++p) {
			for(pagmo::population::size_type q = 0; q < 2*NP; ++q) {
				if (prob.compare_x(P[p], P[q])) { //if p dominates q
					S[p].push_back(q); //add q to the set of solutions dominated by p
				}
				else if (prob.compare_x(P[q], P[p])){
					n[p]++; //increment the doination counter of p
				}
			}

			if (n[p] == 0) { //p belongs to the first front
				F[0].push_back(p);
			}
		}

		int i = 0; //Initialize the front counter
		do {
			std::vector<int> tmp;
			F.push_back(tmp); //add a new empty front
			
			int p,q;
			for(std::vector<int>::size_type j = 0; j < F[i].size(); ++j) { //for each p in F[i]
				p = F[i][j];
				for(std::vector<int>::size_type k = 0; k < S[p].size(); ++k) { //for each q in S[p]
					q = S[p][k];
					n[q]--;
					if (n[q] == 0) { //q belongs to the next front
						F[i+1].push_back(q);
					}
				}
			}
			i++;
		}
		while(F[i].size() > 0); //the next front is not empty

		int last_front = 0;
		std::vector<int>::size_type selected_fronts_size = F[0].size();
		while (selected_fronts_size < NP) {
			last_front++;
			selected_fronts_size += F[last_front].size();
		}


		//crowing-distance-assignment for P
		std::vector<int> I(2*NP,0); //initialize distance
		std::vector<std::pair<decision_vector,int> > tmpP(2*NP);

		for(problem::base::size_type i = 0; i < 2*NP; ++i) {
			tmpP[i].first = P[i];
			tmpP[i].second = i;
		}
		
		double fmax,fmin;
		for(problem::base::f_size_type i = 0; i < prob_f_dimension; ++i) {
			CompareFitness comp(prob,i);
			std::sort(tmpP.begin(), tmpP.end(), comp); //sort using each objective function
			I[tmpP[0].second] = std::numeric_limits<int>::infinity(); //so that boundary points are always selected
			I[tmpP[2*NP-1].second] = std::numeric_limits<int>::infinity();
			fmin = prob.objfun(tmpP[0].first)[i];
			fmax = prob.objfun(tmpP[2*NP-1].first)[i];
			for(problem::base::size_type j= 2; j < 2*NP-1; ++j) {
				if (I[tmpP[j].second] != std::numeric_limits<int>::infinity()) { //for all non boundary points
					I[tmpP[j].second] += (I[tmpP[j-1].second] - I[tmpP[j+1].second]) / (fmax-fmin);
				}
			}
		}
		CompareDistance comp_dist(I);
		std::sort(F[last_front].begin(), F[last_front].end(), comp_dist); //sort the last front that fits using the crowded-comparison operator

		{
		pagmo::population::size_type i = 0;
		std::vector<int>::size_type  j = 0;
		while(i < NP) { //set the new population, selecting the best fronts and the best individuals of the last front that fits
			for(std::vector<int>::size_type k=0; k < F[j].size() && i < NP; ++k) {
				pop.set_x(i,P[F[j][k]]);			
				++i;
			}
			++j;
		}
		}
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
	s << "CR:" << m_cr << ' ';
	s << "M:" << m_m << ' ';
	s << "mutation type:" << m_mut.m_type << ' ';
	s << "crossover type:" << m_cro << ' ';

	return s.str();
}

}} //namespaces
