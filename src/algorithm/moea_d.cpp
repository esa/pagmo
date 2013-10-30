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
#include<limits>

#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/math/special_functions/binomial.hpp>

#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../island.h"
#include "../population.h"
#include "../problem/decompose.h"
#include "../util/discrepancy.h"
#include "../util/neighbourhood.h"
#include "../types.h"
#include "base.h"
#include "moea_d.h"

namespace pagmo { namespace algorithm {

moead::moead(int gen,
		   pagmo::problem::decompose::method_type method,
		   population::size_type T, 
		   weight_generation_type weight_generation,
		   double realb,
		   unsigned int limit
		   )
	  :base(),
	  m_gen(gen),
	  m_method(method),
	  m_T(T),
	  m_weight_generation(weight_generation),
	  m_realb(realb),
	  m_limit(limit)
{
	// Sanity checks
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}

	if(m_weight_generation != RANDOM && m_weight_generation != GRID && m_weight_generation != LOW_DISCREPANCY) {
		pagmo_throw(value_error,"non existing weight generation method");
	}
}

/// Clone method.
base_ptr moead::clone() const
{
	return base_ptr(new moead(*this));
}

//Recursive function building all m-ple of elements of X summing to s
void moead::reksum(std::vector<std::vector<double> > &retval,
		   const std::vector<unsigned int>& X,
		   unsigned int m,
		   unsigned int s,
		   std::vector<double> eggs) const {

	if (m==1) {
		if (std::find(X.begin(),X.end(),s) == X.end()) { //not found
			return;
		} else {
			eggs.push_back(s);
			retval.push_back(eggs);
		}
	} else {
		for (unsigned int i=0; i<X.size(); ++i) {
			eggs.push_back(X[i]);
			reksum(retval,X,m-1,s-X[i],eggs);
			eggs.pop_back();
		}
	}
}

/// Generates the weights used in the problem decomposition
 std::vector<fitness_vector> moead::generate_weights(const unsigned int n_f, const unsigned int n_w) const {

	// Definition of useful probability distributions
	boost::uniform_real<double> uniform(0.0,1.0);
	boost::variate_generator<boost::lagged_fibonacci607 &, boost::uniform_real<double> > r_dist(m_drng,uniform);

	std::vector<fitness_vector> retval;
	if(m_weight_generation == GRID) {
			//find the largest H resulting in a population smaller or equal to NP
			unsigned int H;
			if (n_f == 2) {
				H = n_w-1;
			} else if (n_f == 3) {
				H = floor(0.5 * (sqrt(8*n_w + 1) - 3));
			} else {
				std::cout << "Fitness dimension is " << n_f << std::endl;
				H = 1;
				while(boost::math::binomial_coefficient<double>(H+n_f-1, n_f-1) <= n_w) {
					++H;
				}
				H--;
			}

			// We check that NP equals the population size resulting from H
			if (fabs(n_w-(boost::math::binomial_coefficient<double>(H+n_f-1, n_f-1))) > 1E-8) {
				std::ostringstream error_message;
				error_message << "Invalid population size. Select " << boost::math::binomial_coefficient<double>(H+n_f-1, n_f-1)
						<< " or " << boost::math::binomial_coefficient<double>(H+1+n_f-1, n_f-1)
						<< ".";
				pagmo_throw(value_error,error_message.str());
			}
	
			// We generate the weights
			std::vector<unsigned int> range;
			for (unsigned int i=0; i<H+1;++i) {
				range.push_back(i);
			}
			double epsilon = 1E-6;
			reksum(retval, range, n_f, H);
			for(unsigned int i=0; i< retval.size(); ++i) {
				for(unsigned int j=0; j< retval[i].size(); ++j) {
					retval[i][j] += epsilon;  //NOTE: to avoid to have any weight exactly equal to zero
					retval[i][j] /= H+epsilon*retval[i].size();
				}
			}
	
		} else if(m_weight_generation == LOW_DISCREPANCY) {
			pagmo::util::discrepancy::simplex generator(n_f,1);
			for(unsigned int i = 0; i <n_w; ++i) {
				retval.push_back(generator());
			}
	
		} else if(m_weight_generation == RANDOM) {
			pagmo::util::discrepancy::project_2_simplex projection(n_f);
			for (unsigned int i = 0; i<n_w; ++i) {
				fitness_vector dummy(n_f-1,0.0);
				for(unsigned int j = 0; j <n_f-1; ++j) {
					dummy[j] = r_dist();
				}
				retval.push_back(projection(dummy));
			}
		}
		return retval;
 }

// Performs polynomial mutation (code from nsgaII)
void moead::mutation(decision_vector& child, const pagmo::population& pop, double rate) const
{

	problem::base::size_type D = pop.problem().get_dimension();
	const decision_vector &lb = pop.problem().get_lb(), &ub = pop.problem().get_ub();
	double rnd, delta1, delta2, mut_pow, deltaq;
	double y, yl, yu, val, xy;
	double eta_m = 20;

	//This implements the real polinomial mutation of an individual
	for (pagmo::problem::base::size_type j=0; j < D; ++j){
		if (m_drng() <= rate) {
			y = child[j];
			yl = lb[j];
			yu = ub[j];
			delta1 = (y-yl)/(yu-yl);
			delta2 = (yu-y)/(yu-yl);
			rnd = m_drng();
			mut_pow = 1.0/(eta_m+1.0);
			if (rnd <= 0.5)
			{
				xy = 1.0-delta1;
				val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(eta_m+1.0)));
				deltaq =  pow(val,mut_pow) - 1.0;
			}
			else
			{
				xy = 1.0-delta2;
				val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(eta_m+1.0)));
				deltaq = 1.0 - (pow(val,mut_pow));
			}
			y = y + deltaq*(yu-yl);
			if (y<yl) y = yl;
			if (y>yu) y = yu;
			child[j] = y;
		}
	}
}

void moead::mating_selection(std::vector<population::size_type> &list, int n, int type,const std::vector<std::vector<population::size_type> >& neigh_idx) const {
		list.clear();
		std::vector<population::size_type>::size_type ss   = neigh_idx[n].size(), p;
		
		boost::uniform_int<int> idx(0,neigh_idx.size()-1);
		boost::variate_generator<boost::mt19937 &, boost::uniform_int<int> > p_idx(m_urng,idx);
		
		while(list.size()<2)
		{
			if(type==1){
				p = neigh_idx[n][p_idx()%ss];
			} else {
				p = p_idx();
			}
	
			bool flag = true;
			for(std::vector<population::size_type>::size_type i=0; i<list.size(); i++)
			{
				if(list[i]==p) // p is in the list
				{ 
					flag = false;
					break;
				}
			}
			if(flag) list.push_back(p);
		}
}

/// Evolve implementation.

void moead::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const population::size_type NP = pop.size();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();

	// And make some sanity checks
	if ( prob.get_f_dimension() < 2 ) {
		pagmo_throw(value_error, "The problem is not multiobjective, try some other algorithm than MOEA/D");
	}
	if ( m_T > NP-1 ) {
		pagmo_throw(value_error, "Too many neighbours specified");
	}

	// Get out if there is nothing to do.
	if (m_gen == 0) {
		return;
	}
	
	// Variate generators
	boost::uniform_int<int> pop_idx(0,NP-1);
	boost::variate_generator<boost::mt19937 &, boost::uniform_int<int> > p_idx(m_urng,pop_idx);
	
	// Initialize the candidate chromosome
	decision_vector candidate(prob.get_dimension()); 

	// Compute the ideal point
	fitness_vector ideal_point = pop.compute_ideal();

	// Generate the weights for NP decomposed problems
	std::vector<fitness_vector> weights = generate_weights(prob.get_f_dimension(), NP);
	
	// We compute, for each weight vector, the m_T neighbouring ones
	std::vector<std::vector<population::size_type> > neigh_idx;
	pagmo::util::neighbourhood::euclidian::compute_neighbours(neigh_idx, weights);
	for (unsigned int i=0; i < neigh_idx.size();++i) {
		neigh_idx[i].erase(neigh_idx[i].begin()+m_T, neigh_idx[i].end());
	}

	// We create the decomposed problem which we will use to compute the decomposed fitness values according
	// to weights and ideal points (the construction parameter weights[0] is thus irrelevant)
	pagmo::problem::decompose prob_decomposed(prob, m_method, weights[0], ideal_point);

	// We create a pseudo-random permutation of the indexes 1..NP
	std::vector<population::size_type> shuffle(NP);
	for(pagmo::population::size_type i=0; i < shuffle.size(); ++i) shuffle[i] = i;

	// Main MOEA/D loop
	for (int g = 0; g<m_gen; ++g) {
	//Shuffle the indexes
	std::random_shuffle(shuffle.begin(), shuffle.end(), p_idx);
		for (population::size_type i = 0; i<NP;++i) {
			// We consider the subproblem with index n
			int n = shuffle[i];
			// We select at random between a neighborhood and the whole pop
			int type;
			if(m_drng()<m_realb)	type = 1;   // neighborhood
			else					type = 2;   // whole population

			// 1 - We select two mating partners (not n) in the neighbourhood
			std::vector<population::size_type> p(2);
			mating_selection(p,n,type,neigh_idx);
			
			// 2 - We produce and evaluate an offspring using a DE operator
			double rate =0.5;
			for(decision_vector::size_type kk=0;kk<prob.get_dimension(); ++kk)
			{
				/*Selected Two Parents*/
				candidate[kk] = pop.get_individual(n).cur_x[kk] + rate*(pop.get_individual(p[0]).cur_x[kk] - pop.get_individual(p[1]).cur_x[kk]);
				
				// Fix the bounds
				if(candidate[kk]<lb[kk]){
					candidate[kk] = lb[kk] + m_drng()*(pop.get_individual(n).cur_x[kk] - lb[kk]);
				}
				if(candidate[kk]>ub[kk]){ 
					candidate[kk] = ub[kk] - m_drng()*(ub[kk] - pop.get_individual(n).cur_x[kk]);
				}
			}
			mutation(candidate, pop, 1.0 / prob.get_dimension());
			fitness_vector new_f = prob.objfun(candidate);
			m_fevals++;
			
			// 3 - We update the ideal point
			for (fitness_vector::size_type j=0; j<prob.get_f_dimension(); ++j){
				if (new_f[j] < ideal_point[j]) ideal_point[j] = new_f[j];
			}
			prob_decomposed.set_ideal_point(ideal_point);
			
			// 4-  We insert the newly found solution into the population
			unsigned int size, time = 0;
			// First try on problem n
			fitness_vector f1(1),f2(1);
			prob_decomposed.compute_decomposed_fitness(f1,pop.get_individual(n).cur_f,weights[n]);
			prob_decomposed.compute_decomposed_fitness(f2,new_f,weights[n]);
			if(f2[0]<f1[0])
			{
				pop.set_x(n,candidate);
				time++;
			}
			// Then on neighbouring problems up to m_limit (to preserve diversity)
			if(type==1)	size = neigh_idx[n].size();	// neighborhood
			else		size = NP;					// whole population
			std::vector<population::size_type> shuffle2(size);
			for(pagmo::population::size_type k=0; k < shuffle2.size(); ++k) shuffle2[k] = k;
			std::random_shuffle(shuffle2.begin(), shuffle2.end(), p_idx);
			for (unsigned int k=0; k<shuffle2.size(); ++k) {
				population::size_type pick;
				if(type==1)	pick = neigh_idx[n][shuffle2[k]];		// neighborhood
				else		pick = shuffle2[k];					// whole population

				prob_decomposed.compute_decomposed_fitness(f1,pop.get_individual(pick).cur_f,weights[pick]);
				prob_decomposed.compute_decomposed_fitness(f2,new_f,weights[pick]);
				if(f2[0]<f1[0])
				{
					pop.set_x(pick,candidate);
					time++;
				}
				// the maximal number of solutions updated is not allowed to exceed 'limit'
				if(time>=m_limit) break;
			}
		}
	}
}

/// Algorithm name
std::string moead::get_name() const
{
	return "MOEA-D-DE (MOEA-D-DE)";
}

/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string moead::human_readable_extra() const
{
	std::ostringstream s;
	s << "gen:" << m_gen << ' ';
	s << "neighbours:" << m_T << ' ';
	s << "decomposition:";
	switch (m_method)
	{
		case pagmo::problem::decompose::BI : s << "BI" << ' ';
			break;
		case pagmo::problem::decompose::WEIGHTED : s << "WEIGHTED" << ' ';
			break;
		case pagmo::problem::decompose::TCHEBYCHEFF : s << "TCHEBYCHEFF" << ' ';
			break;
	}
	s << "weights:";
	switch (m_weight_generation)
	{
		case RANDOM : s << "RANDOM" << ' ';
			break;
		case LOW_DISCREPANCY : s << "LOW_DISCREPANCY" << ' ';
			break;
		case GRID : s << "GRID" << ' ';
			break;
	}
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::moead)
