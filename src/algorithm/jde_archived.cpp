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
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <string>
#include <vector>

#include "../exceptions.h"
#include "../population.h"
#include "../types.h"
#include "base.h"
#include "jde_archived.h"

#define KPARAM 5

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] gen number of generations.
 * @param[in] variant algorithm variant (one of 1..18)
 * @param[in] variant_adptv parameter adaptation scheme to be used (one of 1..2)
 * @param[in] ftol stopping criteria on the x tolerance
 * @param[in] xtol stopping criteria on the f tolerance
 * @param[in] memory when true the algorithm preserves its internal state (adapted parameters) through successive calls
 * @param[in] maxarchive maximum size of the archive
 * @param[in] threshold (Novelty) threshold of a solution needed to be added into the archive
 * @throws value_error if f,cr are not in the [0,1] interval, strategy is not one of 1 .. 10, gen is negative
 */
jde_archived::jde_archived(int gen, int variant, int variant_adptv, double ftol, double xtol, bool memory, int maxarchive, double threshold ):base(), m_gen(gen), m_f(0), m_cr(0),
	 m_variant(variant), m_variant_adptv(variant_adptv), m_ftol(ftol), m_xtol(xtol), m_memory(memory), m_max_archive_size(maxarchive), m_archive_threshold(threshold) {
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if (variant < 1 || variant > 18) {
		pagmo_throw(value_error,"variant index must be one of 1 ... 18");
	}
	if (variant_adptv < 1 || variant_adptv > 2) {
		pagmo_throw(value_error,"adaptive variant index must be one of 1 ... 2");
	}
    if (maxarchive < 1) {
        pagmo_throw(value_error,"maximum archive size must be at least 1!");
    }
}

/// Clone method.
base_ptr jde_archived::clone() const
{
	return base_ptr(new jde_archived(*this));
}

/// Distance metric in behaviour space (at the moment: Euclidean distance)
double jde_archived::cdist(pagmo::fitness_vector a, pagmo::fitness_vector b, problem::base::size_type dim) const
{
    double retval = 0.0;
    for (fitness_vector::size_type i = 0; i < dim; ++i) {
        retval += std::pow(a[i] - b[i], 2.0);
    }
    return std::sqrt(retval);
}

/// Evolve implementation.
/**
 * Run the jDE algorithm for the number of generations specified in the constructors.
 * At each improvments velocity is also updated.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */
void jde_archived::evolve(population &pop) const
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

// We use the fitness dimension to pass coordinates as value for the objective function
//	if ( prob_f_dimension != 1 ) {
//		pagmo_throw(value_error,"The problem is not single objective and DE is not suitable to solve it");
//	}

	if (NP < 8) {
		pagmo_throw(value_error,"for jDE at least 8 individuals in the population are needed");
	}

	// Get out if there is nothing to do.
	if (m_gen == 0) {
		return;
	}
	// Some vectors used during evolution are allocated here.
	decision_vector dummy(D), tmp(D); //dummy is used for initialisation purposes, tmp to contain the mutated candidate
	std::vector<decision_vector> popold(NP,dummy), popnew(NP,dummy);
	decision_vector gbX(D),gbIter(D);
	//fitness_vector newfitness(2);	//new fitness of the mutaded candidate
	//fitness_vector gbfit(2);	//global best fitness
	//Changed for experimenting with novelty search
	fitness_vector newfitness(4);	//new fitness of the mutaded candidate
	fitness_vector gbfit(4);	//global best fitness
	std::vector<fitness_vector> fitold(NP,gbfit);
    std::vector<fitness_vector> fitnew(NP,gbfit);
    std::vector<double> dist;   // used for saving distances in behavior space
    std::vector<double> dist_copy;   // a copy for better searching the k-nearest neighbours
    std::vector<int> ind;   // used for saving indices of k-nearest neighbours
    double minval = 0.0; // saves the minimum value
    int minind = 0; // saves th minimum index
    double acc = 0.0; // used as an accumulator
    std::vector<double> old_novelty(NP, 0.0);    // saves the novelty score of old individuals
    std::vector<double> new_novelty(NP, 0.0);    // saves the novelty score of new individuals

	//We extract from pop the chromosomes and fitness associated
	for (std::vector<double>::size_type i = 0; i < NP; ++i) {
		popold[i] = pop.get_individual(i).cur_x;
		fitold[i] = pop.get_individual(i).cur_f;
	}
	popnew = popold;

	// Initialise the global bests
	gbX=pop.champion().x;
	gbfit=pop.champion().f;
	// container for the best decision vector of generation
	gbIter = gbX;
	
	// Initializing the random number generators
	boost::normal_distribution<double> normal(0.0,1.0);
	boost::variate_generator<boost::lagged_fibonacci607 &, boost::normal_distribution<double> > n_dist(m_drng,normal);
	boost::uniform_real<double> uniform(0.0,1.0);
	boost::variate_generator<boost::lagged_fibonacci607 &, boost::uniform_real<double> > r_dist(m_drng,uniform);
	boost::uniform_int<int> r_p_idx(0,NP-1);
	boost::variate_generator<boost::mt19937 &, boost::uniform_int<int> > p_idx(m_urng,r_p_idx);
	boost::uniform_int<int> r_c_idx(0,Dc-1);
	boost::variate_generator<boost::mt19937 &, boost::uniform_int<int> > c_idx(m_urng,r_c_idx);

	
	// Initialize the F and CR vectors
	if ( (m_cr.size() != NP) || (m_f.size() != NP) || (m_memory) ) {
		m_cr.resize(NP); m_f.resize(NP);
		if (m_variant_adptv==1) {
			for (size_t i = 0; i < NP; ++i) {
				m_cr[i] = r_dist();
				m_f[i]  = r_dist() * 0.9 + 0.1;
			}
		}
		else if (m_variant_adptv==2) {
			for (size_t i = 0; i < NP; ++i) {
				m_cr[i] = n_dist() * 0.15 + 0.5;
				m_f[i]  = n_dist() * 0.15 + 0.5;
			}
		}
	}
	// We initialize the global best for F and CR as the first individual (this will soon be forgotten)
	double gbIterF = m_f[0];
	double gbIterCR = m_cr[0];

	// Main DE iterations
	size_t r1,r2,r3,r4,r5,r6,r7;	//indexes to the selected population members
	for (int gen = 0; gen < m_gen; ++gen) {
		//Start of the loop through the genom
		for (size_t i = 0; i < NP; ++i) {
			do {                       /* Pick a random population member */
				/* Endless loop for NP < 2 !!!     */
				r1 = p_idx();
			} while (r1==i);

			do {                       /* Pick a random population member */
				/* Endless loop for NP < 3 !!!     */
				r2 = p_idx();
			} while ((r2==i) || (r2==r1));

			do {                       /* Pick a random population member */
				/* Endless loop for NP < 4 !!!     */
				r3 = p_idx();
			} while ((r3==i) || (r3==r1) || (r3==r2));

			do {                       /* Pick a random population member */
				/* Endless loop for NP < 5 !!!     */
				r4 = p_idx();
			} while ((r4==i) || (r4==r1) || (r4==r2) || (r4==r3));

			do {                       /* Pick a random population member */
				/* Endless loop for NP < 6 !!!     */
				r5 = p_idx();
			} while ((r5==i) || (r5==r1) || (r5==r2) || (r5==r3) || (r5==r4));
			do {                       /* Pick a random population member */
				/* Endless loop for NP < 7 !!!     */
				r6 = p_idx();
			} while ((r6==i) || (r6==r1) || (r6==r2) || (r6==r3) || (r6==r4) || (r6==r5));
			do {                       /* Pick a random population member */
				/* Endless loop for NP < 8 !!!     */
				r7 = p_idx();
			} while ((r7==i) || (r7==r1) || (r7==r2) || (r7==r3) || (r7==r4) || (r7==r5) || (r7==r6));

			// Adapt amplification factor and crossover probability
			double F=0, CR=0;
			if (m_variant_adptv==1) {
				F =  (r_dist() < 0.9) ? m_f[i]  : r_dist() * 0.9 + 0.1;
				CR = (r_dist() < 0.9) ? m_cr[i] : r_dist();
			}
					
			/*-------DE/best/1/exp--------------------------------------------------------------------*/
			/*-------Our oldest strategy but still not bad. However, we have found several------------*/
			/*-------optimization problems where misconvergence occurs.-------------------------------*/
			if (m_variant == 1) { /* strategy DE0 (not in our paper) */
				if (m_variant_adptv==2) {
					F =  gbIterF  + n_dist() * 0.5 * (m_f[r2]-m_f[r3]);
					CR = gbIterCR + n_dist() * 0.5 * (m_cr[r2]-m_cr[r3]);
				}
				tmp = popold[i];
				size_t n = c_idx(), L = 0;
				do {
					tmp[n] = gbIter[n] + F*(popold[r2][n]-popold[r3][n]);
					n = (n+1)%Dc;
					++L;
				} while ((r_dist() < CR) && (L < Dc));
			}

			/*-------DE/rand/1/exp-------------------------------------------------------------------*/
			/*-------This is one of my favourite strategies. It works especially well when the-------*/
			/*-------"gbIter[]"-schemes experience misconvergence. Try e.g. m_f=0.7 and m_cr=0.5---------*/
			/*-------as a first guess.---------------------------------------------------------------*/
			else if (m_variant == 2) { /* strategy DE1 in the techreport */
				if (m_variant_adptv==2) {
					F =  m_f[r1]  + n_dist() * 0.5 * (m_f[r2]-m_f[r3]);
					CR = m_cr[r1] + n_dist() * 0.5 * (m_cr[r2]-m_cr[r3]);
				}
				tmp = popold[i];
				size_t n = c_idx(), L = 0;
				do {
					tmp[n] = popold[r1][n] + F*(popold[r2][n]-popold[r3][n]);
					n = (n+1)%Dc;
					++L;
				} while ((r_dist() < CR) && (L < Dc));
			}

			/*-------DE/rand-to-best/1/exp-----------------------------------------------------------*/
			/*-------This strategy seems to be one of the best strategies. Try m_f=0.85 and c=1.------*/
			/*-------If you get misconvergence try to increase NP. If this doesn't help you----------*/
			/*-------should play around with all three control variables.----------------------------*/
			else if (m_variant == 3) { /* similiar to DE2 but generally better */
				if (m_variant_adptv==2) {
					F =  m_f[i]  + n_dist() * 0.5 * (gbIterF-m_f[i]) +   n_dist() * 0.5 * (m_f[r1] - m_f[r2]);
					CR = m_cr[i] + n_dist() * 0.5 * (gbIterCR-m_cr[i]) + n_dist() * 0.5 * (m_cr[r1] - m_cr[r2]);
				}
				tmp = popold[i];
				size_t n = c_idx(), L = 0;
				do {
					tmp[n] = tmp[n] + F*(gbIter[n] - tmp[n]) + F*(popold[r1][n]-popold[r2][n]);
					n = (n+1)%Dc;
					++L;
				} while ((r_dist() < CR) && (L < Dc));
			}
			/*-------DE/best/2/exp is another powerful strategy worth trying--------------------------*/
			else if (m_variant == 4) {
				if (m_variant_adptv==2) {
					F =  gbIterF  + n_dist() * 0.5 *  (m_f[r1] - m_f[r3]) +   n_dist() * 0.5 * (m_f[r2] - m_f[r4]);
					CR = gbIterCR  + n_dist() * 0.5 * (m_cr[r1] - m_cr[r3]) + n_dist() * 0.5 * (m_cr[r2] - m_cr[r4]);
				}
				tmp = popold[i];
				size_t n = c_idx(), L = 0;
				do {
					tmp[n] = gbIter[n] +
						 (popold[r1][n]+popold[r2][n]-popold[r3][n]-popold[r4][n])*F;
					n = (n+1)%Dc;
					++L;
				} while ((r_dist() < CR) && (L < Dc));
			}
			/*-------DE/rand/2/exp seems to be a robust optimizer for many functions-------------------*/
			else if (m_variant == 5) {
				if (m_variant_adptv==2) {
					F =  m_f[r5]  + n_dist() * 0.5 * (m_f[r1] - m_f[r3]) +   n_dist() * 0.5 * (m_f[r2] - m_f[r4]);
					CR = m_cr[r5] + n_dist() * 0.5 * (m_cr[r1] - m_cr[r3]) + n_dist() * 0.5 * (m_cr[r2] - m_cr[r4]);
				}
				tmp = popold[i];
				size_t n = c_idx(), L = 0;
				do {
					tmp[n] = popold[r5][n] +
						 (popold[r1][n]+popold[r2][n]-popold[r3][n]-popold[r4][n])*F;
					n = (n+1)%Dc;
					++L;
				} while ((r_dist() < CR) && (L < Dc));
			}

			/*=======Essentially same strategies but BINOMIAL CROSSOVER===============================*/

			/*-------DE/best/1/bin--------------------------------------------------------------------*/
			else if (m_variant == 6) {
				if (m_variant_adptv==2) {
					F =  gbIterF  + n_dist() * 0.5 * (m_f[r2]-m_f[r3]);
					CR = gbIterCR + n_dist() * 0.5 * (m_cr[r2]-m_cr[r3]);
				}
				tmp = popold[i];
				size_t n = c_idx();
				for (size_t L = 0; L < Dc; ++L) { /* perform Dc binomial trials */
					if ((r_dist() < CR) || L + 1 == Dc) { /* change at least one parameter */
						tmp[n] = gbIter[n] + F*(popold[r2][n]-popold[r3][n]);
					}
					n = (n+1)%Dc;
				}
			}
			/*-------DE/rand/1/bin-------------------------------------------------------------------*/
			else if (m_variant == 7) {
				if (m_variant_adptv==2) {
					F =  m_f[r1]  + n_dist() * 0.5 * (m_f[r2]-m_f[r3]);
					CR = m_cr[r1] + n_dist() * 0.5 * (m_cr[r2]-m_cr[r3]);
				}
				tmp = popold[i];
				size_t n = c_idx();
				for (size_t L = 0; L < Dc; ++L) { /* perform Dc binomial trials */
					if ((r_dist() < CR) || L + 1 == Dc) { /* change at least one parameter */
						tmp[n] = popold[r1][n] + F*(popold[r2][n]-popold[r3][n]);
					}
					n = (n+1)%Dc;
				}
			}
			/*-------DE/rand-to-best/1/bin-----------------------------------------------------------*/
			else if (m_variant == 8) {
				if (m_variant_adptv==2) {
					F =  m_f[i]  + n_dist() * 0.5 * (gbIterF-m_f[i]) +   n_dist() * 0.5 * (m_f[r1] - m_f[r2]);
					CR = m_cr[i] + n_dist() * 0.5 * (gbIterCR-m_cr[i]) + n_dist() * 0.5 * (m_cr[r1] - m_cr[r2]);
				}
				tmp = popold[i];
				size_t n = c_idx();
				for (size_t L = 0; L < Dc; ++L) { /* perform Dc binomial trials */
					if ((r_dist() < CR) || L + 1 == Dc) { /* change at least one parameter */
						tmp[n] = tmp[n] + F*(gbIter[n] - tmp[n]) + F*(popold[r1][n]-popold[r2][n]);
					}
					n = (n+1)%Dc;
				}
			}
			/*-------DE/best/2/bin--------------------------------------------------------------------*/
			else if (m_variant == 9) {
				if (m_variant_adptv==2) {
					F =  gbIterF  + n_dist() * 0.5 *  (m_f[r1] - m_f[r3]) +   n_dist() * 0.5 * (m_f[r2] - m_f[r4]);
					CR = gbIterCR  + n_dist() * 0.5 * (m_cr[r1] - m_cr[r3]) + n_dist() * 0.5 * (m_cr[r2] - m_cr[r4]);
				}
				tmp = popold[i];
				size_t n = c_idx();
				for (size_t L = 0; L < Dc; ++L) { /* perform Dc binomial trials */
					if ((r_dist() < CR) || L + 1 == Dc) { /* change at least one parameter */
						tmp[n] = gbIter[n] +
							 (popold[r1][n]+popold[r2][n]-popold[r3][n]-popold[r4][n])*F;
					}
					n = (n+1)%Dc;
				}
			}
			/*-------DE/rand/2/bin--------------------------------------------------------------------*/
			else if (m_variant == 10) {
				if (m_variant_adptv==2) {
					F =  m_f[r5]  + n_dist() * 0.5 * (m_f[r1] - m_f[r3]) +   n_dist() * 0.5 * (m_f[r2] - m_f[r4]);
					CR = m_cr[r5] + n_dist() * 0.5 * (m_cr[r1] - m_cr[r3]) + n_dist() * 0.5 * (m_cr[r2] - m_cr[r4]);
				}
				tmp = popold[i];
				size_t n = c_idx();
				for (size_t L = 0; L < Dc; ++L) { /* perform Dc binomial trials */
					if ((r_dist() < CR) || L + 1 == Dc) { /* change at least one parameter */
						tmp[n] = popold[r5][n] +
							 (popold[r1][n]+popold[r2][n]-popold[r3][n]-popold[r4][n])*F;
					}
					n = (n+1)%Dc;
				}
			}

			/*-------DE/best/3/exp--------------------------------------------------------------------*/
			if (m_variant == 11) {
				if (m_variant_adptv==2) {
					F =  gbIterF  + n_dist() * 0.5 * (m_f[r1] - m_f[r2]) + n_dist() * 0.5 * (m_f[r3] - m_f[r4]) + n_dist() * 0.5 * (m_f[r5] - m_f[r6]);
					CR = gbIterCR + n_dist() * 0.5 * (m_cr[r1] - m_cr[r2]) + n_dist() * 0.5 * (m_cr[r3] - m_cr[r4]) + n_dist() * 0.5 * (m_cr[r5] - m_cr[r6]);
				}
				tmp = popold[i];
				size_t n = c_idx(), L = 0;
				do {
					tmp[n] = gbIter[n] + F*(popold[r1][n]-popold[r2][n]) + F*(popold[r3][n]-popold[r4][n]) + F*(popold[r5][n]-popold[r6][n]);
					n = (n+1)%Dc;
					++L;
				} while ((r_dist() < CR) && (L < Dc));
			}
			/*-------DE/best/3/bin--------------------------------------------------------------------*/
			else if (m_variant == 12) {
				if (m_variant_adptv==2) {
					F =  gbIterF  + n_dist() * 0.5 * (m_f[r1] - m_f[r2]) + n_dist() * 0.5 * (m_f[r3] - m_f[r4]) + n_dist() * 0.5 * (m_f[r5] - m_f[r6]);
					CR = gbIterCR + n_dist() * 0.5 * (m_cr[r1] - m_cr[r2]) + n_dist() * 0.5 * (m_cr[r3] - m_cr[r4]) + n_dist() * 0.5 * (m_cr[r5] - m_cr[r6]);
				}
				tmp = popold[i];
				size_t n = c_idx();
				for (size_t L = 0; L < Dc; ++L) { /* perform Dc binomial trials */
					if ((r_dist() < CR) || L + 1 == Dc) { /* change at least one parameter */
						tmp[n] = gbIter[n] + F*(popold[r1][n]-popold[r2][n]) + F*(popold[r3][n]-popold[r4][n]) + F*(popold[r5][n]-popold[r6][n]);
					}
					n = (n+1)%Dc;
				}
			}
			/*-------DE/rand/3/exp--------------------------------------------------------------------*/
			if (m_variant == 13) {
				if (m_variant_adptv==2) {
					F =  m_f[r7]  + n_dist() * 0.5 * (m_f[r1] - m_f[r2]) + n_dist() * 0.5 * (m_f[r3] - m_f[r4]) + n_dist() * 0.5 * (m_f[r5] - m_f[r6]);
					CR = m_cr[r7] + n_dist() * 0.5 * (m_cr[r1] - m_cr[r2]) + n_dist() * 0.5 * (m_cr[r3] - m_cr[r4]) + n_dist() * 0.5 * (m_cr[r5] - m_cr[r6]);
				}
				tmp = popold[i];
				size_t n = c_idx(), L = 0;
				do {
					tmp[n] = popold[r7][n] + F*(popold[r1][n]-popold[r2][n]) + F*(popold[r3][n]-popold[r4][n]) + F*(popold[r5][n]-popold[r6][n]);
					n = (n+1)%Dc;
					++L;
				} while ((r_dist() < CR) && (L < Dc));
			}
			/*-------DE/rand/3/bin--------------------------------------------------------------------*/
			else if (m_variant == 14) {
				if (m_variant_adptv==2) {
					F =  m_f[r7]  + n_dist() * 0.5 * (m_f[r1] - m_f[r2]) + n_dist() * 0.5 * (m_f[r3] - m_f[r4]) + n_dist() * 0.5 * (m_f[r5] - m_f[r6]);
					CR = m_cr[r7] + n_dist() * 0.5 * (m_cr[r1] - m_cr[r2]) + n_dist() * 0.5 * (m_cr[r3] - m_cr[r4]) + n_dist() * 0.5 * (m_cr[r5] - m_cr[r6]);
				}
				tmp = popold[i];
				size_t n = c_idx();
				for (size_t L = 0; L < Dc; ++L) { /* perform Dc binomial trials */
					if ((r_dist() < CR) || L + 1 == Dc) { /* change at least one parameter */
						tmp[n] = popold[r7][n] + F*(popold[r1][n]-popold[r2][n]) + F*(popold[r3][n]-popold[r4][n]) + F*(popold[r5][n]-popold[r6][n]);
					}
					n = (n+1)%Dc;
				}
			}
			/*-------DE/rand-to-current/2/exp---------------------------------------------------------*/
			if (m_variant == 15) {
				if (m_variant_adptv==2) {
					F =  m_f[r7]  + n_dist() * 0.5 * (m_f[r1] - m_f[i]) + n_dist() * 0.5 * (m_f[r3] - m_f[r4]) + n_dist() * 0.5 * (m_f[r5] - m_f[r6]);
					CR = m_cr[r7] + n_dist() * 0.5 * (m_cr[r1] - m_cr[i]) + n_dist() * 0.5 * (m_cr[r3] - m_cr[r4]) + n_dist() * 0.5 * (m_cr[r5] - m_cr[r6]);
				}
				tmp = popold[i];
				size_t n = c_idx(), L = 0;
				do {
					tmp[n] = popold[r7][n] + F*(popold[r1][n]-popold[i][n]) + F*(popold[r3][n]-popold[r4][n]);
					n = (n+1)%Dc;
					++L;
				} while ((r_dist() < CR) && (L < Dc));
			}
			/*-------DE/rand-to-current/2/bin---------------------------------------------------------*/
			else if (m_variant == 16) {
				if (m_variant_adptv==2) {
					F =  m_f[r7]  + n_dist() * 0.5 * (m_f[r1] - m_f[i]) + n_dist() * 0.5 * (m_f[r3] - m_f[r4]) + n_dist() * 0.5 * (m_f[r5] - m_f[r6]);
					CR = m_cr[r7] + n_dist() * 0.5 * (m_cr[r1] - m_cr[i]) + n_dist() * 0.5 * (m_cr[r3] - m_cr[r4]) + n_dist() * 0.5 * (m_cr[r5] - m_cr[r6]);
				}
				tmp = popold[i];
				size_t n = c_idx();
				for (size_t L = 0; L < Dc; ++L) { /* perform Dc binomial trials */
					if ((r_dist() < CR) || L + 1 == Dc) { /* change at least one parameter */
						tmp[n] = popold[r7][n] + F*(popold[r1][n]-popold[i][n]) + F*(popold[r3][n]-popold[r4][n]);
					}
					n = (n+1)%Dc;
				}
			}
			/*-------DE/rand-to-best-and-current/2/exp------------------------------------------------*/
			if (m_variant == 17) {
				if (m_variant_adptv==2) {
					F =  m_f[r7]  + n_dist() * 0.5 * (m_f[r1] - m_f[i]) + n_dist() * 0.5 * (gbIterF - m_f[r4]);
					CR = m_cr[r7] + n_dist() * 0.5 * (m_cr[r1] - m_cr[i]) + n_dist() * 0.5 * (gbIterCR - m_cr[r4]);
				}
				tmp = popold[i];
				size_t n = c_idx(), L = 0;
				do {
					tmp[n] = popold[r7][n] + F*(popold[r1][n]-popold[i][n]) + F*(gbIter[n]-popold[r4][n]);
					n = (n+1)%Dc;
					++L;
				} while ((r_dist() < CR) && (L < Dc));
			}
			/*-------DE/rand-to-best-and-current/2/bin------------------------------------------------*/
			else if (m_variant == 18) {
				if (m_variant_adptv==2) {
					F =  m_f[r7]  + n_dist() * 0.5 * (m_f[r1] - m_f[i]) + n_dist() * 0.5 * (gbIterF - m_f[r4]);
					CR = m_cr[r7] + n_dist() * 0.5 * (m_cr[r1] - m_cr[i]) + n_dist() * 0.5 * (gbIterCR - m_cr[r4]);
				}
				tmp = popold[i];
				size_t n = c_idx();
				for (size_t L = 0; L < Dc; ++L) { /* perform Dc binomial trials */
					if ((r_dist() < CR) || L + 1 == Dc) { /* change at least one parameter */
						tmp[n] = popold[r7][n] + F*(popold[r1][n]-popold[i][n]) + F*(gbIter[n]-popold[r4][n]);
					}
					n = (n+1)%Dc;
				}
			}

			/*=======Trial mutation now in tmp[]. force feasibility and how good this choice really was.==================*/
			// a) feasibility
			size_t i2 = 0;
			while (i2<Dc) {
				if ((tmp[i2] < lb[i2]) || (tmp[i2] > ub[i2]))
					tmp[i2] = r_dist() * (ub[i2]-lb[i2]) + lb[i2];
				++i2;
			}

            // b) evaluate coordinates in behavior space
            popnew[i] = tmp;
            prob.objfun(fitnew[i], tmp);
            

		}//End of the loop through the deme

        // loop through old population and compute distances
        dist.clear();
        ind.clear();
        for (pagmo::population::size_type i = 0; i < NP; ++i) {
            for (pagmo::population::size_type j = 0; j < NP; ++j) {
                dist.push_back(cdist(fitold[i], fitold[j], prob_f_dimension));
            }
            for (pagmo::population::size_type j = 0; j < NP; ++j) {
                dist.push_back(cdist(fitold[i], fitnew[j], prob_f_dimension));
            }
            for (std::vector<fitness_vector>::size_type j = 0; j < m_archive.size(); ++j) {
                dist.push_back(cdist(fitold[i], m_archive[j], prob_f_dimension));
            }

            // copy distance vector
            dist_copy = dist;

            // determine k nearest neighbours
            for (size_t k = 0; k < KPARAM; ++k) {
                minind = 0;
                minval = dist[0];
                for (size_t j = 1; j < dist_copy.size(); ++j) {
                    if (minval < dist_copy[j]) {
                        minval = dist_copy[j];
                        minind = j;
                    }
                }
                // poison vector with infinity
                dist_copy[minind] = 1.7e+300;
                ind.push_back(minind);
            }
             
            // compute novelty metric
            acc = 0.0;
            for (size_t k = 0; k < KPARAM; ++k) {
                for (size_t j = 0; j < ind.size(); ++j) {
                    acc += dist[ind[j]];
                }
            }
            old_novelty[i] = acc / KPARAM;
        }
            
        // loop through candidate soluations and compute distances
        dist.clear();
        ind.clear();
        for (pagmo::population::size_type i = 0; i < NP; ++i) {
            for (pagmo::population::size_type j = 0; j < NP; ++j) {
                dist.push_back(cdist(fitnew[i], fitold[j], prob_f_dimension));
            }
            for (pagmo::population::size_type j = 0; j < NP; ++j) {
                dist.push_back(cdist(fitnew[i], fitnew[j], prob_f_dimension));
            }
            for (std::vector<fitness_vector>::size_type j = 0; j < m_archive.size(); ++j) {
                dist.push_back(cdist(fitnew[i], m_archive[j], prob_f_dimension));
            }

            // copy distance vector
            dist_copy = dist;

            // determine k nearest neighbours
            for (size_t k = 0; k < KPARAM; ++k) {
                minind = 0;
                minval = dist[0];
                for (size_t j = 1; j < dist_copy.size(); ++j) {
                    if (minval < dist_copy[j]) {
                        minval = dist_copy[j];
                        minind = j;
                    }
                }
                // poison vector with infinity
                dist_copy[minind] = 1.7e+300;
                ind.push_back(minind);
            }
             
            // compute novelty metric
            acc = 0.0;
            for (size_t k = 0; k < KPARAM; ++k) {
                for (size_t j = 0; j < ind.size(); ++j) {
                    acc += dist[ind[j]];
                }
            }
            new_novelty[i] = acc / KPARAM;
        }
        
        // compare old and new via novelty metric
		for (size_t i = 0; i < NP; ++i) {
            if ( old_novelty[i] > new_novelty[i] ) {
//                popnew[i] = popold[i];  // old individual has stronger novelty
				pop.set_x(i, popold[i]);
            } else {
				// adapt parameters of a succesful solution
            //    m_cr[i] = CR;    // WHY YOU NOT WORK?!?
			//	m_f[i] = F;

                // individual change due to higher novelty of new solution
                pop.set_x(i, popnew[i]);
            }
            // check with archive
            if (m_archive.size() < m_max_archive_size) {
                if (new_novelty[i] > m_archive_threshold) {
                    m_archive.push_back(popnew[i]);
                }
            }
        }

		/* Save best population member of current iteration */
		gbIter = gbX;       // today: give solution with highest novelty

		/* swap population arrays. New generation becomes old one */
		//std::swap(popold, popnew);


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
		}
		

	}//end main DE iterations
	if (m_screen_output) {
		std::cout << "Exit condition -- generations > " <<  m_gen << std::endl;
	}

}

/// Algorithm name
std::string jde_archived::get_name() const
{
	return "jDE_archived";
}

/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string jde_archived::human_readable_extra() const
{
	std::ostringstream s;
	s << "gen:" << m_gen << ' ';
	s << "variant:" << m_variant << ' ';
	s << "self_adaptation:" << m_variant_adptv << ' ';
	s << "memory:" << m_memory << ' ';
	s << "ftol:" << m_ftol << ' ';
	s << "xtol:" << m_xtol << ' ';
    s << "maxarchive:" << m_max_archive_size << ' ';
    s << "threshold:" << m_archive_threshold;

	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::jde_archived);
