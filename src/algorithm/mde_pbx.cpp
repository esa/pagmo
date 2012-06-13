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
#include "mde_pbx.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] gen number of generations.
 * @param[in] qperc percentage of population to choose the best vector
 * @param[in] nexp exponent for the powermean
 * @param[in] ftol stopping criteria on the x tolerance
 * @param[in] xtol stopping criteria on the f tolerance
 * @param[in] restart when true the algorithm re-initialize randomly the parameters at each call
 * @throws value_error if f,cr are not in the [0,1] interval, strategy is not one of 1 .. 10, gen is negative
 */

mde_pbx::mde_pbx(int gen, double qperc, double nexp, double ftol, double xtol):base(), m_gen(gen), 
	 m_f(0), m_fsuccess(0), m_fm(0.5), m_cr(0), m_crsuccess(0), m_crm(0.6), m_qperc(qperc), m_nexp(nexp), 
	 m_ftol(ftol), m_xtol(xtol) {
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if (qperc < 0.0 || qperc > 1.0) {
		pagmo_throw(value_error,"percentage of population to choose from must be between 0.0 and 1.0");
	}
	if (nexp == 0.0) {
		pagmo_throw(value_error,"the exponent of the powermean must not be 0.0!");
	}
}

/// Clone method.
base_ptr mde_pbx::clone() const
{
	return base_ptr(new mde_pbx(*this));
}

/// Evolve implementation.
/**
 * Run the jDE algorithm for the number of generations specified in the constructors.
 * At each improvments velocity is also updated.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void mde_pbx::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type D = prob.get_dimension(), 
	      prob_i_dimension = prob.get_i_dimension(), 
	      prob_c_dimension = prob.get_c_dimension(), 
	      prob_f_dimension = prob.get_f_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP = pop.size();
	const problem::base::size_type Dc = D - prob_i_dimension;
	const population::size_type NP_Part = double_to_int::convert(m_qperc * NP);
	
	size_t trials; // used to prevent infinite loops
	const size_t maxtrials = 5000; // maximum number of trials before throwing an error

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

	if (NP < 3) {
		pagmo_throw(value_error,"for this algorithm, at least 3 individuals in the population are needed");
	}

	if (NP_Part < 1) {
		pagmo_throw(value_error,"You need a higher number of individuals in your population for sampling");
	}

	// Get out if there is nothing to do.
	if (m_gen == 0) {
		return;
	}
	// Some vectors used during evolution are allocated here.
	decision_vector dummy(D), tmp(D); 		//dummy is used for initialisation purposes, tmp to contain the mutated candidate
	std::vector<decision_vector> popold(NP,dummy), popnew(NP,dummy);
	std::vector<decision_vector> poppart(NP_Part, dummy);	// vector used for the q%-subset of population
	decision_vector gbX(D),gbIter(D);
	fitness_vector newfitness(1);			//new fitness of the mutaded candidate
	fitness_vector gbfit(1);			//global best fitness
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

	// reserve space for saving successful values for f and cr
	m_fsuccess.reserve(NP);
	m_crsuccess.reserve(NP);
	m_cr.reserve(NP);
	m_f.reserve(NP);
	
	// Initializing the random number generators
	boost::uniform_real<double> uniform(0.0,1.0);
	boost::variate_generator<boost::lagged_fibonacci607 &, boost::uniform_real<double> > r_dist(m_drng,uniform);

	boost::uniform_int<int> r_p_idx(0,NP-1);
	boost::variate_generator<boost::mt19937 &, boost::uniform_int<int> > p_idx(m_urng,r_p_idx);

	boost::uniform_int<int> r_c_idx(0,Dc-1);
	boost::variate_generator<boost::mt19937 &, boost::uniform_int<int> > c_idx(m_urng,r_c_idx);


	// We initialize the global best for F and CR as the first individual (this will soon be forgotten)
	double gbIterF = m_f[0];
	double gbIterCR = m_cr[0];

	// Main DE loop
	size_t r1,r2;	//indexes to the selected population members
	double p;
		
	for (int gen = 0; gen < m_gen; ++gen) {
	  
	 //Check the exit conditions (every 40 generations)
	    if (gen % 40) {
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
	}

	if (m_screen_output) {
		std::cout << "Exit condition -- generations > " <<  m_gen << std::endl;
	}

}

/// Algorithm name
std::string mde_pbx::get_name() const
{
	return "MDE_pBX";
}

/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string mde_pbx::human_readable_extra() const
{
	std::ostringstream s;
	s << "gen:" << m_gen << ' ';
	s << "q percentage:" << m_qperc << ' ';
	s << "power mean exponent:" << m_nexp << ' ';
	s << "ftol:" << m_ftol << ' ';
	s << "xtol:" << m_xtol;

	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::mde_pbx);
