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
#include <boost/random/cauchy_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <string>
#include <vector>
#include <cmath>

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

mde_pbx::mde_pbx(int gen, double qperc, double nexp, double ftol, double xtol):base(), m_gen(gen), m_fsuccess(0), m_fm(0.5), 
		 m_crsuccess(0), m_crm(0.6), m_qperc(qperc), m_nexp(nexp), m_ftol(ftol), m_xtol(xtol) {
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
* Run the MDE_pBX algorithm for the number of generations specified in the constructors.
* At each improvments velocity is also updated.
*
* @param[in,out] pop input/output pagmo::population to be evolved.
*/
void mde_pbx::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type D = prob.get_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP = pop.size(), NP_Part = double_to_int::convert(m_qperc * NP);

	//We perform some checks to determine wether the problem/population are suitable for DE
	if ( D == 0 ) {
		pagmo_throw(value_error,"Dimension of the problem shall not be zero!");
	}

	if ( prob.get_c_dimension() != 0 ) {
		pagmo_throw(value_error,"The problem is not box constrained and MDE_pBX is not suitable to solve it");
	}

	if ( prob.get_f_dimension() != 1 ) {
		pagmo_throw(value_error,"The problem is not single objective and MDE_pBX is not suitable to solve it");
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
	fitness_vector newfitness(1);			//new fitness of the mutated candidate
	fitness_vector gbfit(1);			//global best fitness
	std::vector<fitness_vector> fit(NP,gbfit);	//saves current fitness of all individuals
	
	//We extract from pop the chromosomes and fitness associated
	for (std::vector<double>::size_type i = 0; i < NP; ++i) {
		popnew[i] = pop.get_individual(i).cur_x;
		fit[i] = pop.get_individual(i).cur_f;
	}
	
	// Initialize the global bests
// 	gbfit=pop.champion().f;

	// reserve space for saving successful values for f and cr
	m_fsuccess.reserve(NP);
	m_crsuccess.reserve(NP);
	
	// Initializing the random number generators
	boost::uniform_real<double> uniform(0.0,1.0);
	boost::variate_generator<boost::lagged_fibonacci607 &, boost::uniform_real<double> > r_dist(m_drng,uniform);

	boost::uniform_int<int> r_c_idx(0,D-1);
	boost::variate_generator<boost::mt19937 &, boost::uniform_int<int> > c_idx(m_urng,r_c_idx);
	
	// Declaring temporary variables used by the main-loop
	pagmo::population::size_type p;
	size_t r1, r2, bestq_idx, bestp_idx, j_rand, trials;
	size_t a[NP];
	double cri, fi, wcr, wf;
	
	// **** Main Loop of differential evolution ****
	for (int gen = 0; gen < m_gen; ++gen) {
	    // starting to evolve a new generation: latest generation becomes the old one
	    std::swap(popold, popnew);
	  
	    // clear the sets of successful scale factors and crossover probabilities
	    m_fsuccess.clear();
	    m_crsuccess.clear();
	    
	    // adjust parameter p controlling the elite of individuals to choose from 
	    // as we start counting with gen=0 and not with gen=1 we use (gen) instead of (gen - 1) here
	    p = ceil((NP / 2.0) * ( 1.0 - (static_cast<double>(gen) / m_gen)));

	    // update random distributions
	    boost::normal_distribution<double> nd(m_crm, 0.1);
	    boost::variate_generator<boost::lagged_fibonacci607 &, boost::normal_distribution<double> > gauss(m_drng,nd);	    
	    
    	    boost::cauchy_distribution<double> cd(m_fm, 0.1);
	    boost::variate_generator<boost::lagged_fibonacci607 &, boost::cauchy_distribution<double> > cauchy(m_drng,cd);	    
	    
	    // loop through all individuals
	    for (pagmo::population::size_type i = 0; i < NP; ++i) {

		// get q% random indices
		for (pagmo::population::size_type k = 0; k < NP; ++k) {
		    a[k] = k;
		}
		
		// We only swap the first q% of the indices
		for (pagmo::population::size_type k = 0; k < NP_Part; ++k) {
		    r1 = double_to_int::convert( boost::uniform_int<int>(k,NP-1)(m_urng) );
		    std::swap(a[k], a[r1]);
		}
		
		// find index of individual from q% sample with best fitness
		bestq_idx = a[0];
		for (pagmo::population::size_type k = 1; k < NP_Part; ++k) {
		    if ( prob.compare_fitness(fit[a[k]], fit[bestq_idx]) ) {
			    bestq_idx = a[k];
		    }
		}

		// choose two random distinct pop members
		do {    
			/* Endless loop for NP < 2 !!!     */
			r1 = boost::uniform_int<pagmo::population::size_type>(0,NP-1)(m_urng);
		} while (r1==i);

		do {            
			/* Endless loop for NP < 3 !!!     */
			r2 = boost::uniform_int<pagmo::population::size_type>(0,NP-1)(m_urng);
		} while ((r2==i) || (r2==r1));
		
		// get a random vector out of the p-best
		std::vector<population::size_type> pbest = pop.get_best_idx(p);
		bestp_idx = pbest[boost::uniform_int<int>(0,p-1)(m_urng)];
		
		// sample scale factors
		trials = 0;
		do {
		  cri = gauss();
		  trials++;
		} while ( ((cri < 0.0) || (cri > 1.0)) && (trials < 20) );
		
		if (trials > 20) {
		    pagmo_throw(value_error,"Random number sampling for Cr is no longer efficient. Evolution aborted.");
		}

		trials = 0;
		do {
		  fi = cauchy();
		  trials++;
		} while ( ((fi <= 0.0) || (fi > 1.0)) && (trials < 20) );
		
		if (trials > 20) {
		    pagmo_throw(value_error,"Random number sampling for F is no longer efficient. Evolution aborted.");
		}

		// fix a random dimension index
		j_rand = c_idx();
		
		// Mutation + Crossover
		tmp = popold[i];
		for (size_t j = 0; j < D; ++j) {
		    if (j == j_rand || r_dist() < cri) {
			tmp[j] = popold[i][j] + fi * (popold[bestq_idx][j] - popold[i][j] + popold[r1][j] - popold[r2][j]);
		    } else {
			tmp[j] = popold[bestp_idx][j];
		    }
		}
		 
		/*=======Trial mutation now in tmp[]. Force feasibility and test how good this choice really was.==========*/
		// a) feasibility (throw in some random value if we fall out of the borders)
		size_t i2 = 0;
		while (i2<D) {
			if ((tmp[i2] < lb[i2]) || (tmp[i2] > ub[i2]))
				tmp[i2] = boost::uniform_real<double>(lb[i2],ub[i2])(m_drng);
			++i2;
		}

		// b) Compare with the objective function
		prob.objfun(newfitness, tmp);    /* Evaluate new vector in tmp[] */
		if ( pop.problem().compare_fitness(newfitness,fit[i]) ) {  /* improved objective function value ? */
			fit[i]=newfitness;
			popnew[i] = tmp;
			// As a fitness improvment occured we move the point
			// and thus can evaluate a new velocity
			std::transform(tmp.begin(), tmp.end(), pop.get_individual(i).cur_x.begin(), tmp.begin(),std::minus<double>());
			// updates x and v (cache avoids to recompute the objective function)
			pop.set_x(i,popnew[i]);
			pop.set_v(i,tmp);
			// remember our successful scale factors
			m_crsuccess.push_back(cri);
			m_fsuccess.push_back(fi);
		} else {
			popnew[i] = popold[i];
		}
		
	    }
		    
	    // Update Crossover Probability
	    if (!m_crsuccess.empty()) {
	      wcr = 0.9  + (0.1 * r_dist());
	      m_crm = (wcr * m_crm) + ((1.0 - wcr) * powermean(m_crsuccess, m_nexp));
	    }
	    
// 	    Update Fitness Scale Factor
	    if (!m_fsuccess.empty()) {
	      wf = 0.8 + (0.2 * r_dist());
	      m_fm = (wf * m_fm) + ((1.0 - wf) * powermean(m_fsuccess, m_nexp));
	    }

	    //Check the exit conditions (every 40 generations)
	    if (gen % 30 == 0) {
		double dx = 0;
		
		for (decision_vector::size_type k = 0; k < D; ++k) {
			tmp[k] = pop.get_individual(pop.get_worst_idx()).best_x[k] - pop.get_individual(pop.get_best_idx()).best_x[k];
			dx += std::fabs(tmp[k]);
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
		    std::cout << "    Fm: " << m_fm << ", Crm: " << m_crm << std::endl;
		}
	    }
	} // End of Generation main iteration

	if (m_screen_output) {
		std::cout << "Exit condition -- generations > " <<  m_gen << std::endl;
	}
}

/// Algorithm name
std::string mde_pbx::get_name() const
{
	return "MDE_pBX";
}

/// Computes the powermean of a set given as a vector
double mde_pbx::powermean(std::vector<double> v, double exp) const
{
	// correct implementation of the power mean
	double sum = 0.0;
	size_t vsize = v.size();

	if (vsize == 0) return 0;
	
	for (size_t i = 0; i < vsize; ++i) {
	    sum += std::pow(v[i], exp) ;
	}
	return std::pow((sum / vsize), (1.0 / exp));

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
