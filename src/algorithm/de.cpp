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
#include <string>
#include <vector>

#include "../exceptions.h"
#include "../population.h"
#include "../types.h"
#include "base.h"
#include "de.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] gen number of generations.
 * @param[in] f weight coefficient (dafault value is 0.8)
 * @param[in] cr crossover probability (dafault value is 0.9)
 * @param[in] strategy strategy (dafault strategy is 2: /rand/1/exp)
 * @param[in] ftol stopping criteria on the x tolerance
 * @param[in] xtol stopping criteria on the f tolerance
 * @throws value_error if f,cr are not in the [0,1] interval, strategy is not one of 1 .. 10, gen is negative
 */
de::de(int gen, double f, double cr, int strategy, double ftol, double xtol):base(),m_gen(gen),m_f(f),m_cr(cr),m_strategy(strategy),m_ftol(ftol),m_xtol(xtol) {
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if (strategy < 1 || strategy > 10) {
		pagmo_throw(value_error,"strategy index must be one of 1 ... 10");
	}
	if (cr < 0 || f < 0 || cr > 1 || f > 1) {
		pagmo_throw(value_error,"the f and cr parameters must be in the [0,1] range");
	}
}

/// Clone method.
base_ptr de::clone() const
{
	return base_ptr(new de(*this));
}

/// Evolve implementation.
/**
 * Run the DE algorithm for the number of generations specified in the constructors.
 * At each improvments velocity is also updated.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void de::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type D = prob.get_dimension(), prob_i_dimension = prob.get_i_dimension(), prob_c_dimension = prob.get_c_dimension(), prob_f_dimension = prob.get_f_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP = pop.size();
	const problem::base::size_type Dc = D - prob_i_dimension;


	//We perform some checks to determine wether the problem/population are suitable for DE
	if ( Dc == 0 ) {
		pagmo_throw(value_error,"There is no continuous part in the problem decision vector for DE to optimise");
	}

	if ( prob_c_dimension != 0 ) {
		pagmo_throw(value_error,"The problem is not box constrained and DE is not suitable to solve it");
	}

	if ( prob_f_dimension != 1 ) {
		pagmo_throw(value_error,"The problem is not single objective and DE is not suitable to solve it");
	}

	if (NP < 6) {
		pagmo_throw(value_error,"for DE at least 6 individuals in the population are needed");
	}

	// Get out if there is nothing to do.
	if (m_gen == 0) {
		return;
	}
	// Some vectors used during evolution are allocated here.
	decision_vector dummy(D), tmp(D); //dummy is used for initialisation purposes, tmp to contain the mutated candidate
	std::vector<decision_vector> popold(NP,dummy), popnew(NP,dummy);
	decision_vector gbX(D),gbIter(D);
	fitness_vector newfitness(prob_f_dimension);	//new fitness of the mutaded candidate
	fitness_vector gbfit(prob_f_dimension);	//global best fitness
	std::vector<fitness_vector> fit(NP,gbfit);

	//We extract from pop the chromosomes and fitness associated
	for (std::vector<double>::size_type i = 0; i < NP; ++i) {
		popold[i] = pop.get_individual(i).cur_x;
		fit[i] = pop.get_individual(i).cur_f;
	}
	popnew = popold;

	// Initialise the global bests
	gbX=pop.champion().x;
	gbfit=pop.champion().f;
	// container for the best decision vector of generation
	gbIter = gbX;

	// Main DE iterations
	size_t r1,r2,r3,r4,r5;	//indexes to the selected population members
	for (int gen = 0; gen < m_gen; ++gen) {
		//Start of the loop through the deme
		for (size_t i = 0; i < NP; ++i) {
			do {                       /* Pick a random population member */
				/* Endless loop for NP < 2 !!!     */
				r1 = boost::uniform_int<int>(0,NP-1)(m_urng);
			} while (r1==i);

			do {                       /* Pick a random population member */
				/* Endless loop for NP < 3 !!!     */
				r2 = boost::uniform_int<int>(0,NP-1)(m_urng);
			} while ((r2==i) || (r2==r1));

			do {                       /* Pick a random population member */
				/* Endless loop for NP < 4 !!!     */
				r3 = boost::uniform_int<int>(0,NP-1)(m_urng);
			} while ((r3==i) || (r3==r1) || (r3==r2));

			do {                       /* Pick a random population member */
				/* Endless loop for NP < 5 !!!     */
				r4 = boost::uniform_int<int>(0,NP-1)(m_urng);
			} while ((r4==i) || (r4==r1) || (r4==r2) || (r4==r3));

			do {                       /* Pick a random population member */
				/* Endless loop for NP < 6 !!!     */
				r5 = boost::uniform_int<int>(0,NP-1)(m_urng);
			} while ((r5==i) || (r5==r1) || (r5==r2) || (r5==r3) || (r5==r4));


			/*-------DE/best/1/exp--------------------------------------------------------------------*/
			/*-------Our oldest strategy but still not bad. However, we have found several------------*/
			/*-------optimization problems where misconvergence occurs.-------------------------------*/
			if (m_strategy == 1) { /* strategy DE0 (not in our paper) */
				tmp = popold[i];
				size_t n = boost::uniform_int<int>(0,Dc-1)(m_urng), L = 0;
				do {
					tmp[n] = gbIter[n] + m_f*(popold[r2][n]-popold[r3][n]);
					n = (n+1)%Dc;
					++L;
				} while ((m_drng() < m_cr) && (L < Dc));
			}

			/*-------DE/rand/1/exp-------------------------------------------------------------------*/
			/*-------This is one of my favourite strategies. It works especially well when the-------*/
			/*-------"gbIter[]"-schemes experience misconvergence. Try e.g. m_f=0.7 and m_cr=0.5---------*/
			/*-------as a first guess.---------------------------------------------------------------*/
			else if (m_strategy == 2) { /* strategy DE1 in the techreport */
				tmp = popold[i];
				size_t n = boost::uniform_int<int>(0,Dc-1)(m_urng), L = 0;
				do {
					tmp[n] = popold[r1][n] + m_f*(popold[r2][n]-popold[r3][n]);
					n = (n+1)%Dc;
					++L;
				} while ((m_drng() < m_cr) && (L < Dc));
			}

			/*-------DE/rand-to-best/1/exp-----------------------------------------------------------*/
			/*-------This strategy seems to be one of the best strategies. Try m_f=0.85 and m_cr=1.------*/
			/*-------If you get misconvergence try to increase NP. If this doesn't help you----------*/
			/*-------should play around with all three control variables.----------------------------*/
			else if (m_strategy == 3) { /* similiar to DE2 but generally better */
				tmp = popold[i];
				size_t n = boost::uniform_int<int>(0,Dc-1)(m_urng), L = 0;
				do {
					tmp[n] = tmp[n] + m_f*(gbIter[n] - tmp[n]) + m_f*(popold[r1][n]-popold[r2][n]);
					n = (n+1)%Dc;
					++L;
				} while ((m_drng() < m_cr) && (L < Dc));
			}
			/*-------DE/best/2/exp is another powerful strategy worth trying--------------------------*/
			else if (m_strategy == 4) {
				tmp = popold[i];
				size_t n = boost::uniform_int<int>(0,Dc-1)(m_urng), L = 0;
				do {
					tmp[n] = gbIter[n] +
						 (popold[r1][n]+popold[r2][n]-popold[r3][n]-popold[r4][n])*m_f;
					n = (n+1)%Dc;
					++L;
				} while ((m_drng() < m_cr) && (L < Dc));
			}
			/*-------DE/rand/2/exp seems to be a robust optimizer for many functions-------------------*/
			else if (m_strategy == 5) {
				tmp = popold[i];
				size_t n = boost::uniform_int<int>(0,Dc-1)(m_urng), L = 0;
				do {
					tmp[n] = popold[r5][n] +
						 (popold[r1][n]+popold[r2][n]-popold[r3][n]-popold[r4][n])*m_f;
					n = (n+1)%Dc;
					++L;
				} while ((m_drng() < m_cr) && (L < Dc));
			}

			/*=======Essentially same strategies but BINOMIAL CROSSOVER===============================*/

			/*-------DE/best/1/bin--------------------------------------------------------------------*/
			else if (m_strategy == 6) {
				tmp = popold[i];
				size_t n = boost::uniform_int<int>(0,Dc-1)(m_urng);
				for (size_t L = 0; L < Dc; ++L) { /* perform Dc binomial trials */
					if ((m_drng() < m_cr) || L + 1 == Dc) { /* change at least one parameter */
						tmp[n] = gbIter[n] + m_f*(popold[r2][n]-popold[r3][n]);
					}
					n = (n+1)%Dc;
				}
			}
			/*-------DE/rand/1/bin-------------------------------------------------------------------*/
			else if (m_strategy == 7) {
				tmp = popold[i];
				size_t n = boost::uniform_int<int>(0,Dc-1)(m_urng);
				for (size_t L = 0; L < Dc; ++L) { /* perform Dc binomial trials */
					if ((m_drng() < m_cr) || L + 1 == Dc) { /* change at least one parameter */
						tmp[n] = popold[r1][n] + m_f*(popold[r2][n]-popold[r3][n]);
					}
					n = (n+1)%Dc;
				}
			}
			/*-------DE/rand-to-best/1/bin-----------------------------------------------------------*/
			else if (m_strategy == 8) {
				tmp = popold[i];
				size_t n = boost::uniform_int<int>(0,Dc-1)(m_urng);
				for (size_t L = 0; L < Dc; ++L) { /* perform Dc binomial trials */
					if ((m_drng() < m_cr) || L + 1 == Dc) { /* change at least one parameter */
						tmp[n] = tmp[n] + m_f*(gbIter[n] - tmp[n]) + m_f*(popold[r1][n]-popold[r2][n]);
					}
					n = (n+1)%Dc;
				}
			}
			/*-------DE/best/2/bin--------------------------------------------------------------------*/
			else if (m_strategy == 9) {
				tmp = popold[i];
				size_t n = boost::uniform_int<int>(0,Dc-1)(m_urng);
				for (size_t L = 0; L < Dc; ++L) { /* perform Dc binomial trials */
					if ((m_drng() < m_cr) || L + 1 == Dc) { /* change at least one parameter */
						tmp[n] = gbIter[n] +
							 (popold[r1][n]+popold[r2][n]-popold[r3][n]-popold[r4][n])*m_f;
					}
					n = (n+1)%Dc;
				}
			}
			/*-------DE/rand/2/bin--------------------------------------------------------------------*/
			else if (m_strategy == 10) {
				tmp = popold[i];
				size_t n = boost::uniform_int<int>(0,Dc-1)(m_urng);
				for (size_t L = 0; L < Dc; ++L) { /* perform Dc binomial trials */
					if ((m_drng() < m_cr) || L + 1 == Dc) { /* change at least one parameter */
						tmp[n] = popold[r5][n] +
							 (popold[r1][n]+popold[r2][n]-popold[r3][n]-popold[r4][n])*m_f;
					}
					n = (n+1)%Dc;
				}
			}


			/*=======Trial mutation now in tmp[]. force feasibility and how good this choice really was.==================*/
			// a) feasibility
			size_t i2 = 0;
			while (i2<Dc) {
				if ((tmp[i2] < lb[i2]) || (tmp[i2] > ub[i2]))
					tmp[i2] = boost::uniform_real<double>(lb[i2],ub[i2])(m_drng);
				++i2;
			}

			//b) how good?
			prob.objfun(newfitness, tmp);    /* Evaluate new vector in tmp[] */
			if ( pop.problem().compare_fitness(newfitness,fit[i]) ) {  /* improved objective function value ? */
				fit[i]=newfitness;
				popnew[i] = tmp;
				// As a fitness improvment occured we move the point
				// and thus can evaluate a new velocity
				std::transform(tmp.begin(), tmp.end(), pop.get_individual(i).cur_x.begin(), tmp.begin(),std::minus<double>());
				//updates x and v (cache avoids to recompute the objective function)
				pop.set_x(i,popnew[i]);
				pop.set_v(i,tmp);
				if ( pop.problem().compare_fitness(newfitness,gbfit) ) {
					/* if so...*/
					gbfit=newfitness;          /* reset gbfit to new low...*/
					gbX=popnew[i];
				}
			} else {
				popnew[i] = popold[i];
			}

		}//End of the loop through the deme

		/* Save best population member of current iteration */
		gbIter = gbX;

		/* swap population arrays. New generation becomes old one */
		std::swap(popold, popnew);


		//9 - Check the exit conditions (every 40 generations)
		if (gen % 40 == 0) {
			double dx = 0;
			for (decision_vector::size_type i = 0; i < D; ++i) {
				tmp[i] = pop.get_individual(pop.get_worst_idx()).best_x[i] - pop.get_individual(pop.get_best_idx()).best_x[i];
				dx += std::fabs(tmp[i]);
			}
			
			if  ( dx < m_xtol ) {
				if (m_screen_output) { 
					std::cout << "Exit condition -- xtol < " <<  m_xtol << std::endl;
				}
				return;
			}

			double mah = std::fabs(pop.get_individual(pop.get_worst_idx()).best_f[0] - pop.get_individual(pop.get_best_idx()).best_f[0]);

			if (mah < m_ftol) {
				if (m_screen_output) {
					std::cout << "Exit condition -- ftol < " <<  m_ftol << std::endl;
				}
				return;
			}

			// outputs current values
			if (m_screen_output) {
				std::cout << "Generation " << gen << " ***" << std::endl;
				std::cout << "    Best global fitness: " << pop.champion().f << std::endl;
				std::cout << "    xtol: " << dx << ", ftol: " << mah << std::endl;
			}

		}
		

	}//end main DE iterations
	if (m_screen_output) {
		std::cout << "Exit condition -- generations > " <<  m_gen << std::endl;
	}

}

/// Algorithm name
std::string de::get_name() const
{
	return "Differential Evolution";
}

/// Sets crossover parameter.
/**
 * 
 * @param[in] cr new value for the crossover parameter
 * @throws value_error if cr not in (0,1)
 */
void de::set_cr(double cr) {
	if (cr < 0 || cr > 1 ) {
		pagmo_throw(value_error,"the cr parameter must be in the [0,1] range");
	}
	m_cr = cr;
}

/// Sets f parameter.
/**
 * 
 * @param[in] f new value for the crossover parameter
 * @throws value_error if f not in (0,1)
 */
void de::set_f(double f) {
	if (f < 0 || f > 1 ) {
		pagmo_throw(value_error,"the cr parameter must be in the [0,1] range");
	}
	m_f = f;
}

/// Gets crossover parameter.
/**
 * 
 * @return the crossover parameter
 */
double de::get_cr() const {
	return m_cr;
}


/// Gets f parameter.
/**
 * 
 * @return the f parameter
 */
double de::get_f() const {
	return m_f;
}


/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string de::human_readable_extra() const
{
	std::ostringstream s;
	s << "gen:" << m_gen << ' ';
	s << "F:" << m_f << ' ';
	s << "CR:" << m_cr << ' ';
	s << "variant:" << m_strategy << ' ';
	s << "ftol:" << m_ftol << ' ';
	s << "xtol:" << m_xtol;
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::de)
