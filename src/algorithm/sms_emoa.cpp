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
#include "sms_emoa.h"

namespace pagmo { namespace algorithm {

/// Constructor
 /**
 * Constructs the SMS-EMOA algorithm.
 *
 * @param[in] gen Number of generations to evolve.
 * @param[in] sel_m Selection method of the algorithm. If sel_m=1, least_contributor is used for the computation regardless of the number of fronts.
 * If sel_m=2, domination count is used in case there are more than one front.
 * @param[in] cr Crossover probability
 * @param[in] eta_c Distribution index for crossover
 * @param[in] m Mutation probability
 * @param[in] eta_m Distribution index for mutation
 * @throws value_error if gen is negative, crossover probability is not \f$ \in [0,1[\f$, mutation probability or mutation width is not \f$ \in [0,1]\f$,
 */
sms_emoa::sms_emoa(int gen, int sel_m, double cr, double eta_c, double m, double eta_m):base(),m_gen(gen),
	m_sel_m(sel_m),m_cr(cr),m_eta_c(eta_c),m_m(m),m_eta_m(eta_m)
{
	validate_parameters();
}

/// Copy constructor
 /**
  * Copy constructor of SMS-EMOA.
  *
  * @param[in] orig Original instance of SMS-EMOA to make a copy from
  */
sms_emoa::sms_emoa(const sms_emoa &orig) : base(),m_gen(orig.m_gen),m_sel_m(orig.m_sel_m),
	m_cr(orig.m_cr),m_eta_c(orig.m_eta_c),m_m(orig.m_m),m_eta_m(orig.m_eta_m)
{
	if (orig.m_hv_algorithm) {
		m_hv_algorithm = orig.m_hv_algorithm->clone();
	}
}

/// Constructor
 /**
 * Constructs the SMS-EMOA algorithm with the preferred hypervolume algorithm for computation.
 *
 * @param[in] gen Number of generations to evolve.
 * @param[in] sel_m Selection method of the algorithm. If sel_m=1, least_contributor is used for the computation regardless of the number of fronts.
 * If sel_m=2, domination count is used in case there are more than one front.
 * @param[in] cr Crossover probability
 * @param[in] eta_c Distribution index for crossover
 * @param[in] m Mutation probability
 * @param[in] eta_m Distribution index for mutation
 * @param[in] hv_algorithm Hypervolume algorithm used for the computation of the least contributor
 * @throws value_error if gen is negative, crossover probability is not \f$ \in [0,1[\f$, mutation probability or mutation width is not \f$ \in [0,1]\f$,
 */
sms_emoa::sms_emoa(pagmo::util::hv_algorithm::base_ptr hv_algorithm, int gen, int sel_m, double cr, double eta_c, double m, double eta_m):base(),
	m_gen(gen),m_sel_m(sel_m), m_cr(cr),m_eta_c(eta_c),m_m(m),m_eta_m(eta_m)
{
	m_hv_algorithm = hv_algorithm;
	validate_parameters();
}

/// Get the hypervolume algorithm used for the computation
/**
 * Returns an instance of the algorithm currently set for the computation of the least_contributor.
 *
 * @return shared pointer to the hv_algorithm instance
 */
pagmo::util::hv_algorithm::base_ptr sms_emoa::get_hv_algorithm() const
{
	return m_hv_algorithm;
}

/// Clone method.
base_ptr sms_emoa::clone() const
{
	return base_ptr(new sms_emoa(*this));
}

void sms_emoa::validate_parameters()
{
	if (m_gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if (m_sel_m < 1 || m_sel_m > 2) {
		pagmo_throw(value_error,"selection method must be equal to 1 (least contributor) or 2 (domination count)");
	}
	if (m_cr >= 1 || m_cr < 0) {
		pagmo_throw(value_error,"crossover probability must be in the [0,1] range");
	}
	if (m_m < 0 || m_m > 1) {
		pagmo_throw(value_error,"mutation probability must be in the [0,1] range");
	}
	if (m_eta_c <1 || m_eta_c >= 100) {
		pagmo_throw(value_error,"Distribution index for crossover must be in 1..100");
	}
	if (m_eta_m <1 || m_eta_m >= 100) {
		pagmo_throw(value_error,"Distribution index for mutation must be in 1..100");
	}
}

void sms_emoa::crossover(decision_vector& child1, decision_vector& child2, pagmo::population::size_type parent1_idx, pagmo::population::size_type parent2_idx,const pagmo::population& pop) const
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

void sms_emoa::mutate(decision_vector& child, const pagmo::population& pop) const
{

	problem::base::size_type D = pop.problem().get_dimension();
	problem::base::size_type Di = pop.problem().get_i_dimension();
	problem::base::size_type Dc = D - Di;
	const decision_vector &lb = pop.problem().get_lb(), &ub = pop.problem().get_ub();
	double rnd, delta1, delta2, mut_pow, deltaq;
	double y, yl, yu, val, xy;
	int gen_num;
	
	//This implements the real polynomial mutation of an individual
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
				deltaq = pow(val,mut_pow) - 1.0;
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
	for (pagmo::problem::base::size_type j=Dc; j < D; ++j) {
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

// Find the index of the least contributing individual
population::size_type sms_emoa::evaluate_s_metric_selection(const population & pop) const
{

	std::vector< std::vector< population::size_type> > fronts = pop.compute_pareto_fronts();

	const std::vector< population::size_type> &last_front = fronts.back();

	if (last_front.size() == 1) {
		return last_front[0];
	}

	// if the chosen method is to always to pick the least contributor, or when working solely on the first front
	if (m_sel_m == 1 || fronts.size() == 1) {
		std::vector<fitness_vector> points;
		points.resize(last_front.size());

		for (population::size_type idx = 0 ; idx < last_front.size() ; ++idx) {
			points[idx] = fitness_vector(pop.get_individual(last_front[idx]).cur_f);
		}

		pagmo::util::hypervolume hypvol(points);
		fitness_vector r = hypvol.get_nadir_point(1.0);

		population::size_type least_idx;

		if (m_hv_algorithm) {
			least_idx = hypvol.least_contributor(r, m_hv_algorithm);
		} else {
			least_idx = hypvol.least_contributor(r);
		}

		return last_front[least_idx];
	} else { // if m_sel_m == 2 && fronts.size() > 1
		population::size_type max_dom_count = 0;
		population::size_type individual_idx = 0;

		for (population::size_type idx = 0 ; idx < last_front.size() ; ++idx) {
			population::size_type current_dom_count = pop.get_domination_count(last_front[idx]);
			if (current_dom_count > max_dom_count) {
				max_dom_count = current_dom_count;
				individual_idx = idx;
			}
		}
		return last_front[individual_idx];
	}
}


/// Evolve implementation.
/**
 * Run the SMS-EMOA evolve routine.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void sms_emoa::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type D = prob.get_dimension();
	const problem::base::size_type prob_c_dimension = prob.get_c_dimension();
	const population::size_type NP = pop.size();

	if ( prob_c_dimension != 0 ) {
		pagmo_throw(value_error, "The problem is not box constrained and SMS-EMOA is not suitable to solve it");
	}
	
	if ( prob.get_f_dimension() < 2 ) {
		pagmo_throw(value_error, "The problem is not multiobjective, try some other algorithm than SMS-EMOA");
	}

	// Get out if there is nothing to do.
	if (m_gen == 0) {
		return;
	}
	
	population::size_type parent1_idx, parent2_idx;
	decision_vector child1(D), child2(D);
	
	// Main SMS-EMOA loop
	for (int g = 0; g < m_gen; g++) {
		// select two different parent indices from the population
		parent1_idx = m_urng() % NP;
		parent2_idx = ((m_urng() % (NP-1)) + parent1_idx) % NP;

		crossover(child1, child2, parent1_idx, parent2_idx, pop);
		++m_fevals;
		mutate(child1, pop);
		pop.push_back(child1);
		pop.erase(evaluate_s_metric_selection(pop));
	}
}

/// Algorithm name
std::string sms_emoa::get_name() const
{
	return "S-Metric Selection Evolutionary Multiobjective Optimisation Algorithm (SMS-EMOA)";
}

/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string sms_emoa::human_readable_extra() const
{
	std::ostringstream s;
	s << "gen:" << m_gen << ' ';
	s << "cr:" << m_cr << ' ';	
	s << "eta_c:" << m_eta_c << ' ';
	s << "m:" << m_m << ' ';
	s << "eta_m:" << m_eta_m << ' ';
	s << "hv_algorithm:";
	if (m_hv_algorithm) {
		s << m_hv_algorithm->get_name();
	} else {
		s << "Chosen dynamically";
	}
	s << std::endl;

	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::sms_emoa)
