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
#include "../population.h"
#include "../util/neighbourhood.h"
#include "base.h"
#include "spea2.h"


namespace pagmo { namespace algorithm {
/// Constructor
 /**
 * Constructs a SPEA2 algorithm
 *
 * @param[in] gen Number of generations to evolve.
 * @param[in] cr Crossover probability
 * @param[in] eta_c Distribution index for crossover
 * @param[in] m Mutation probability
 * @param[in] eta_m Distribution index for mutation
 *
 * @throws value_error if gen is negative
 */
spea2::spea2(int gen, double cr, double eta_c, double m, double eta_m):base(),
	m_gen(gen),m_cr(cr),m_eta_c(eta_c),m_m(m),m_eta_m(eta_m)
{
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
}

/// Clone method.
base_ptr spea2::clone() const
{
	return base_ptr(new spea2(*this));
}

/// Evolve implementation.
/**
 * Run the SPEA2 algorithm for the number of generations specified in the constructors.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */
void spea2::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base             &prob = pop.problem();
	const problem::base::size_type   D = prob.get_dimension(), prob_i_dimension = prob.get_i_dimension(), prob_c_dimension = prob.get_c_dimension(), prob_f_dimension = prob.get_f_dimension();
	const problem::base::size_type   Dc = D - prob_i_dimension;
	const population::size_type      NP = pop.size();

	//We perform some checks to determine wether the problem/population are suitable for PSO
	if( Dc == 0 ){
		pagmo_throw(value_error,"There is no continuous part in the problem decision vector for spea2 to optimise");
	}

	if( prob_c_dimension != 0 ){
		pagmo_throw(value_error,"The problem is not box constrained and spea2 is not suitable to solve it");
	}

	if( prob_f_dimension < 2 ){
		pagmo_throw(value_error,"The problem is not multi-objective. Use a single-objectie optimization algorithm instead");
	}

	if (NP < 5 || (NP % 4 != 0) ) {
		pagmo_throw(value_error, "for NSGA-II at least 5 individuals in the population are needed and the population size must be a multiple of 4");
	}

	// Get out if there is nothing to do.
	if (NP == 0 || m_gen == 0) {
		return;
	}

	const population::size_type archive_size = NP;

	std::vector<population::size_type> bestIndices_archive(archive_size, 0);

	population archive(prob,0);
	population new_archive(prob,0);

	for(int g = 0; g < m_gen; ++g) {
		std::cout << "gen: " << g << std::endl;


		new_archive.clear();
		for(unsigned int i=0; i < NP; ++i) {
			new_archive.push_back(pop.get_individual(i).cur_x);
		}
		for(unsigned int i=0; i < archive.size(); ++i) {
			new_archive.push_back(archive.get_individual(i).cur_x);
		}


		std::vector<fitness_vector> fit(new_archive.size());
		for ( population::size_type i = 0; i<new_archive.size(); i++ ) {
			fit[i]	=	new_archive.get_individual(i).cur_f;
		}
		std::vector<std::vector<pagmo::population::size_type> > neighbours;
		pagmo::util::neighbourhood::euclidian::compute_neighbours(neighbours, fit);

		std::vector<double> F_archive(new_archive.size(),0);// individuals' fitness (according to raw fitness and density)
		compute_spea2_fitness(F_archive, neighbours, fit, new_archive);

		std::vector<population::size_type> order_by_fitness = pagmo::util::neighbourhood::order(F_archive);

		unsigned int n_non_dominated;
		for(n_non_dominated = 0; n_non_dominated < new_archive.size() && F_archive[order_by_fitness[n_non_dominated]] < 1; ++n_non_dominated);

		if(n_non_dominated > archive_size) { //truncate according to delta

			std::vector<fitness_vector> fit_nd(n_non_dominated);
			for ( population::size_type i = 0; i<n_non_dominated; i++ ) {
				fit_nd[i]	=	fit[order_by_fitness[i]];
			}
			std::vector<std::vector<pagmo::population::size_type> > neighbours_nd;
			pagmo::util::neighbourhood::euclidian::compute_neighbours(neighbours_nd, fit_nd);

			std::vector<population::size_type> rv(n_non_dominated);
			for(unsigned int i=0; i<n_non_dominated; ++i) rv[i] = i;

			order_by_delta(rv, neighbours_nd, fit_nd, static_cast<int>(sqrt(2*n_non_dominated)));

			for(unsigned int i=0; i<archive_size; ++i) {
				bestIndices_archive[i] = order_by_fitness[rv[i]];
			}
		} else { //fill with dominated individuals
			bestIndices_archive = std::vector<population::size_type>(order_by_fitness.begin(), order_by_fitness.begin()+archive_size);
		}

		//Update the new archive
		archive.clear();
		for(unsigned int i=0; i<archive_size; ++i) {
			archive.push_back(new_archive.get_individual(bestIndices_archive[i]).cur_x);
		}

		std::vector<population::size_type> shuffle1(archive_size), shuffle2(archive_size);
		for (pagmo::population::size_type i=0; i< archive_size; i++) {
			shuffle1[i] = i;
			shuffle2[i] = i;
		}
		boost::uniform_int<int> pop_idx(0,archive_size-1);
		boost::variate_generator<boost::mt19937 &, boost::uniform_int<int> > p_idx(m_urng,pop_idx);

		population::size_type parent1_idx, parent2_idx;
		decision_vector child1(D), child2(D);

		//We create some pseudo-random permutation of the poulation indexes
		std::random_shuffle(shuffle1.begin(),shuffle1.end(), p_idx);
		std::random_shuffle(shuffle2.begin(),shuffle2.end(), p_idx);

		// We completely cancel the population (NOTE: memory of all individuals and the notion of
		// champion is thus destroyed)
		pop.clear();

		//We then loop thorugh all individuals with increment 4 to select two pairs of parents that will
		//each create 2 new offspring
		for (pagmo::population::size_type i=0; i<archive_size; i+=4) {
			// We create two offsprings using the shuffled list 1
			parent1_idx = tournament_selection(bestIndices_archive[shuffle1[i]], bestIndices_archive[shuffle1[i+1]],F_archive);
			parent2_idx = tournament_selection(bestIndices_archive[shuffle1[i+2]], bestIndices_archive[shuffle1[i+3]],F_archive);
			crossover(child1, child2, parent1_idx,parent2_idx,new_archive);
			mutate(child1,new_archive);
			mutate(child2,new_archive);
			pop.push_back(child1);
			pop.push_back(child2);
		}

		//same with shuffle2
		for (pagmo::population::size_type i=0; i<archive_size; i+=4) {
			// We create two offsprings using the shuffled list 1
			parent1_idx = tournament_selection(bestIndices_archive[shuffle2[i]], bestIndices_archive[shuffle2[i+1]],F_archive);
			parent2_idx = tournament_selection(bestIndices_archive[shuffle2[i+2]], bestIndices_archive[shuffle2[i+3]],F_archive);
			crossover(child1, child2, parent1_idx,parent2_idx,new_archive);
			mutate(child1,new_archive);
			mutate(child2,new_archive);
			pop.push_back(child1);
			pop.push_back(child2);
		}

	}

}

void spea2::order_by_delta(std::vector<population::size_type> &rv, const std::vector<std::vector<pagmo::population::size_type>  > &neighbours, const std::vector<fitness_vector> &fit, int K) const
{
	if(rv.size() != fit.size() ) {
		std::cout << "order_by_delta: rv size should be equal to fit size, returning doing nothing";
		return;
	} else {
		std::sort(rv.begin(), rv.end(), distance_sorter(neighbours, fit, K));
	}
}

void spea2::compute_spea2_fitness(std::vector<double> &F,
			const std::vector<std::vector<pagmo::population::size_type> > &neighbours,
			const std::vector<fitness_vector> &fit,
			const pagmo::population &pop) const {
	std::vector<std::vector<population::size_type> > domination_list;
	const population::size_type NP = pop.size();
	std::vector<double> S(NP, 0.0);
	const int K = sqrt(2*NP); //TODO: this should be global, check where else it is used



	for(unsigned i=0; i<NP; ++i) {
		domination_list.push_back(pop.get_domination_list(i));
	}

	for(unsigned int i=0; i<NP; ++i) {
		S[i] = static_cast<double>(domination_list[i].size()) / (NP+1);
	}

	std::fill(F.begin(), F.end(), 0);

	for(unsigned int i=0; i<NP; ++i) {
		for(unsigned int j=0; j<domination_list[i].size(); ++j) {
			F[domination_list[i][j]] += S[i];
		}
	}

	for(unsigned int i=0; i<NP; ++i) {
		F[i] = S[i] +
				(1.0 / (pagmo::util::neighbourhood::euclidian::distance(fit[i], fit[neighbours[i][K]]) + 2));
	}
}

pagmo::population::size_type spea2::tournament_selection(pagmo::population::size_type idx1, pagmo::population::size_type idx2, const std::vector<double> &F) const
{
	if (F[idx1] < F[idx2]) {
		return idx1;
	}
	else if (F[idx1] > F[idx2]) {
		return idx2;
	}
	else {
		return ((m_drng() > 0.5) ? idx1 : idx2);
	}
}

void spea2::crossover(decision_vector& child1, decision_vector& child2, pagmo::population::size_type parent1_idx, pagmo::population::size_type parent2_idx,const pagmo::population& pop) const
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

void spea2::mutate(decision_vector& child, const pagmo::population& pop) const
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

/// Algorithm name
std::string spea2::get_name() const
{
	return "Non-dominated Sorting Particle Swarm Optimizer (spea2)";
}

/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string spea2::human_readable_extra() const
{
	std::ostringstream s;
	s << "gen:" << m_gen << ' ';
	return s.str();
}


}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::spea2);
