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

#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../types.h"
#include "vega.h"
#include "../problem/base_stochastic.h"
#include "../problem/decompose.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the vega multi-objective algorithm.
 *
 * @param[in] gen Number of generations to evolve.
 * @param[in] cr Crossover probability (of each allele if binomial crossover)
 * @param[in] m Mutation probability (of each allele)
 * @param[in] elitism The best individual is reinserted in the population each elitism generations
 * @param[in] mut Mutation type. One of sga::mutation::GAUSSIAN, sga::mutation::RANDOM
 * @param[in] width Mutation width. When gaussian mutation is selected is the width of the mutation
 * @param[in] cro Crossover type. One of sga::crossover::BINOMIAL, sga::crossover::EXPONENTIAL
 * @throws value_error if gen is negative, crossover probability is not \f$ \in [0,1]\f$, mutation probability is not \f$ \in [0,1]\f$,
 * elitism is <=0
 *
 */
vega::vega(int gen, const double &cr, const double &m, int elitism, mutation::type mut, double width, crossover::type cro)
	:base(),m_gen(gen),m_cr(cr),m_m(m),m_elitism(elitism),m_mut(mut,width),m_cro(cro)
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
	if (elitism < 1) {
		pagmo_throw(value_error,"elitisim must be greater than zero");
	}
	if (width <0 || width >1) {
		pagmo_throw(value_error,"mutation width must be in the [0,1] range");
	}
}

/// Clone method.
base_ptr vega::clone() const
{
	return base_ptr(new vega(*this));
}

/// Evolve implementation.
/**
 * Run the vega algorithm for the number of generations specified in the constructors.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void vega::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type D = prob.get_dimension(), Di = prob.get_i_dimension(), prob_c_dimension = prob.get_c_dimension(), prob_f_dimension = prob.get_f_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP = pop.size();
	const problem::base::size_type Dc = D - Di;

	//We perform some checks to determine wether the problem/population are suitable for VEGA
	if( prob_c_dimension != 0 ) {
		pagmo_throw(value_error,"The problem is not box constrained and VEGA is not suitable to solve it");
	}

	if( prob_f_dimension < 2 ) {
		pagmo_throw(value_error,"The problem is not multiobjective and VEGA is not suitable to solve it");
	}

	if(NP < 5*prob_f_dimension) {
		pagmo_throw(value_error,"for VEGA at least 5 * number of objectives individuals in the population are needed");
	}

	if( (NP%prob_f_dimension != 0) && (NP != 1) ) {
		pagmo_throw(value_error,"for VEGA the population size must a multiple of the number of objectives");
	}

	// sub population size
	population::size_type sub_pop_size = NP/prob_f_dimension;

	// Get out if there is nothing to do.
	if (m_gen == 0) {
		return;
	}

	// Some vectors used during evolution are allocated here.
	decision_vector dummy(D,0);			//used for initialisation purposes
	std::vector<decision_vector > X(NP,dummy);
	std::vector<decision_vector > Xnew(NP,dummy);
	// fitness for the main problem
	std::vector<fitness_vector > fit(NP);

	// selection vectors
	std::vector<double> selectionfitness(sub_pop_size);
	std::vector<double> cumsum(sub_pop_size);
	std::vector<double> cumsumTemp(sub_pop_size);

	std::vector<population::size_type> selection(sub_pop_size);

	// Initialise the chromosomes and their fitness to that of the initial deme
	for (pagmo::population::size_type i = 0; i<NP; i++ ) {
		X[i]	=	pop.get_individual(i).cur_x;
		fit[i]	=	pop.get_individual(i).cur_f;
	}

	// Find the best member and store in bestX and bestfit
	fitness_vector bestfit;
	decision_vector bestX(D,0);
	population::size_type bestidx = pop.get_best_idx();
	bestX = pop.get_individual(bestidx).cur_x;
	bestfit = pop.get_individual(bestidx).cur_f;

	// creates new sub problems using the decompose meta-problem
	std::vector<problem::base_ptr> sub_probs;
	std::vector< std::vector<decision_vector> > sub_x(prob_f_dimension);
	std::vector< std::vector<fitness_vector> > sub_f(prob_f_dimension);

	// generating the subproblems used to evaluate fitnesses
	for(unsigned int i=0; i<prob_f_dimension; i++) {
		std::vector<double> sub_prob_weights(prob_f_dimension,0.);
		sub_prob_weights[i] = 1.;
		sub_probs.push_back(problem::decompose(prob, problem::decompose::WEIGHTED, sub_prob_weights).clone());
	}

	// Main VEGA loop
	for(int j=0; j<m_gen; j++) {

		boost::uniform_int<int> pop_idx(0,NP-1);
		boost::variate_generator<boost::mt19937 &, boost::uniform_int<int> > p_idx(m_urng,pop_idx);

		//We create some pseudo-random permutation of the poulation indexes
		std::random_shuffle(X.begin(),X.end(),p_idx);

		// Initialise the fitness with each problem
		for(unsigned int sp_idx=0; sp_idx<prob_f_dimension; sp_idx++) {
			problem::base_ptr current_sub_prob = sub_probs.at(sp_idx);
			std::vector<fitness_vector> &current_sub_f = sub_f[sp_idx];
			std::vector<decision_vector> &current_sub_x = sub_x[sp_idx];

			current_sub_f = std::vector<fitness_vector>(sub_pop_size);
			current_sub_x = std::vector<decision_vector>(sub_pop_size);

			population::size_type sub_pop_begin_idx = sp_idx*sub_pop_size;

			// for each individual of the current subpopulation
			for(pagmo::population::size_type k=0; k<sub_pop_size; k++) {
				// initialize design and fitness vectors
				current_sub_x[k] = X.at(sub_pop_begin_idx + k);
				current_sub_f[k] = current_sub_prob->objfun(current_sub_x[k]);
			}
		}

		// loop over the sub problems to do the selection
		for(unsigned int sp_idx=0; sp_idx<prob_f_dimension; sp_idx++) {

			std::vector<fitness_vector> &current_sub_f = sub_f[sp_idx];
			problem::base_ptr current_sub_prob = sub_probs.at(sp_idx);

			// ---------------------------
			// scaled roulette
			// selection of the best individuals among the current sub population
			// checked and seems to be correctly implemented
			// ---------------------------

			// We scale all fitness values from 0 (worst) to absolute value of the best fitness
			fitness_vector worstfit = current_sub_f[0];
			for(pagmo::population::size_type i=1; i < sub_pop_size; i++) {
				if(current_sub_prob->compare_fitness(worstfit, current_sub_f[i]))
					worstfit = current_sub_f[i];
			}

			for(pagmo::population::size_type i=0; i < sub_pop_size; i++) {
				selectionfitness[i] = fabs(worstfit[0] - current_sub_f[i][0]);
			}

			// We build and normalise the cumulative sum
			cumsumTemp[0] = selectionfitness[0];

			for(pagmo::population::size_type i=1; i < sub_pop_size; i++) {
				cumsumTemp[i] = cumsumTemp[i - 1] + selectionfitness[i];
			}
			for(pagmo::population::size_type i=0; i < sub_pop_size; i++) {
				cumsum[i] = cumsumTemp[i] / cumsumTemp[sub_pop_size-1];
			}

			//we throw a dice and pick up the corresponding index
			double r2;
			for(pagmo::population::size_type i=0; i < sub_pop_size; i++) {
				r2 = m_drng();
				for(pagmo::population::size_type l=0; l < sub_pop_size; l++) {
					if(cumsum[l] > r2) {
						selection[i] = l;
						break;
					}
				}
			}
			// ---------------------------
			// end of roulette selection
			// ---------------------------

			// We have the selected individuals in the sub population
			// we can now store the current individuals in the mating pool
			population::size_type sub_pop_begin_idx = sp_idx * sub_pop_size;
			std::vector<decision_vector> &current_sub_x = sub_x[sp_idx];

			// Xnew stores the new selected generation of chromosomes
			for (pagmo::population::size_type i=0; i < sub_pop_size; i++) {
				Xnew[sub_pop_begin_idx + i] = current_sub_x[selection[i]];
			}
		}

		//2 - Crossover
		{
			int r1, L;
			decision_vector member1, member2;

			for(pagmo::population::size_type i=0; i < NP; i++) {
				// for each chromosome selected i.e. in Xnew
				member1 = Xnew[i];
				// we select a mating partner different from the self (i.e. no masturbation)
				do {
					r1 = boost::uniform_int<int>(0,NP - 1)(m_urng);
				} while ( r1 == boost::numeric_cast<int>(i) );
				member2 = Xnew[r1];
				// and we operate crossover
				switch (m_cro) {
				// 0 - binomial crossover
				case crossover::BINOMIAL: {
					size_t n = boost::uniform_int<int>(0,D-1)(m_urng);
					for (size_t L = 0; L < D; ++L) { /* perform D binomial trials */
						if ((m_drng() < m_cr) || L + 1 == D) { /* change at least one parameter */
							member1[n] = member2[n];
						}
						n = (n+1)%D;
					}
					break; }
					// 1 - exponential crossover
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
		}

		//3 - Mutation
		switch (m_mut.m_type) {
		case mutation::GAUSSIAN: {
			boost::normal_distribution<double> dist;
			boost::variate_generator<boost::lagged_fibonacci607 &, boost::normal_distribution<double> > delta(m_drng,dist);
			for (pagmo::problem::base::size_type k = 0; k < Dc;k++) { //for each continuous variable
				double std = (ub[k]-lb[k]) * m_mut.m_width;
				for (pagmo::population::size_type i = 0; i < NP;i++) { //for each individual
					if (m_drng() < m_m) {
						double mean = Xnew[i][k];
						double tmp = (delta() * std + mean);
						if ( (tmp < ub[k]) &&  (tmp > lb[k]) ) Xnew[i][k] = tmp;
					}
				}
			}
			for (pagmo::problem::base::size_type k = Dc; k < D;k++) { //for each integer variable
				double std = (ub[k]-lb[k]) * m_mut.m_width;
				for (pagmo::population::size_type i = 0; i < NP;i++) { //for each individual
					if (m_drng() < m_m) {
						double mean = Xnew[i][k];
						double tmp = boost::math::iround(delta() * std + mean);
						if ( (tmp < ub[k]) &&  (tmp > lb[k]) ) Xnew[i][k] = tmp;
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

		// If the problem is a stochastic optimization change the seed and re-evaluate taking care to update also best and local bests
		try
		{
			//4 - Evaluate the new population (stochastic problem)
			dynamic_cast<const pagmo::problem::base_stochastic &>(prob).set_seed(m_urng());
			pop.clear(); // Removes memory based on different seeds (champion and best_x, best_f, best_c)

			// We re-evaluate the best individual (for elitism)
			prob.objfun(bestfit,bestX);
			// Re-evaluate wrt new seed the particle position and memory
			for (pagmo::population::size_type i=0; i < NP; i++) {
				// We evaluate here the new individual fitness
				prob.objfun(fit[i],Xnew[i]);
				// We update the velocity (in case coupling with PSO via archipelago)
				//dummy = Xnew[i];
				//std::transform(dummy.begin(), dummy.end(), pop.get_individual(i).cur_x.begin(), dummy.begin(),std::minus<double>());
				///We now set the cleared pop. cur_x is the best_x, re-evaluated with new seed.
				pop.push_back(Xnew[i]);
				//pop.set_v(i,dummy);
				if (prob.compare_fitness(fit[i], bestfit)) {
					bestfit = fit[i];
					bestX = Xnew[i];
				}
			}
		}
		catch (const std::bad_cast& e)
		{
			//4 - Evaluate the new population (deterministic problem)
			for (pagmo::population::size_type i=0; i < NP; i++) {
				prob.objfun(fit[i],Xnew[i]);
				dummy = Xnew[i];
				std::transform(dummy.begin(), dummy.end(), pop.get_individual(i).cur_x.begin(), dummy.begin(),std::minus<double>());
				//updates x and v (cache avoids to recompute the objective function and constraints)
				pop.set_x(i,Xnew[i]);
				pop.set_v(i,dummy);
				if (prob.compare_fitness(fit[i], bestfit)) {
					bestfit = fit[i];
					bestX = Xnew[i];
				}
			}
		}

		// need to check if elitism is suitable with MO

		//5 - Reinsert best individual every m_elitism generations
		if (j % m_elitism == 0) {
			int worst=0;
			for (pagmo::population::size_type i = 1; i < NP;i++) {
				if ( prob.compare_fitness(fit[worst],fit[i]) ) worst=i;
			}
			Xnew[worst] = bestX;
			fit[worst] = bestfit;
			dummy = Xnew[worst];
			std::transform(dummy.begin(), dummy.end(), pop.get_individual(worst).cur_x.begin(), dummy.begin(),std::minus<double>());
			//updates x and v (cache avoids to recompute the objective function)
			pop.set_x(worst,Xnew[worst]);
			pop.set_v(worst,dummy);
		}

		// updating the design vectors
		X = Xnew;
	} // end of main VEGA loop
}

/// Algorithm name
std::string vega::get_name() const
{
	return "Vector evaluated genetic algorithm (VEGA)";
}

/// Extra human readable algorithm info.
/**
									 * @return a formatted string displaying the parameters of the algorithm.
									 */
std::string vega::human_readable_extra() const
{
	std::ostringstream s;
	s << "gen:" << m_gen << ' ';
	s << "CR:" << m_cr << ' ';
	s << "M:" << m_m << ' ';
	s << "elitism:" << m_elitism << ' ';
	s << "mutation:";
	switch (m_mut.m_type) {
	case mutation::RANDOM: {
		s << "RANDOM ";
		break;
	}
	case mutation::GAUSSIAN: {
		s << "GAUSSIAN (" << m_mut.m_width << ") ";
		break;
	}
	}
	s << "crossover:";
	switch (m_cro) {
	case crossover::EXPONENTIAL: {
		s << "EXPONENTIAL ";
		break;
	}
	case crossover::BINOMIAL: {
		s << "BINOMIAL ";
		break;
	}
	}

	return s.str();
}

//pagmo::population::size_type vega::roulette_selection(pagmo::population::size_type idx1, pagmo::population::size_type idx2, const pagmo::population& pop) const
//{

//}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::vega)
