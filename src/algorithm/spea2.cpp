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
spea2::spea2(int gen, double cr, double eta_c, double m, double eta_m, int archive_size):base(),
	m_gen(gen),m_cr(cr),m_eta_c(eta_c),m_m(m),m_eta_m(eta_m),m_archive_size(archive_size)
{
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if (archive_size < -1) {
		pagmo_throw(value_error,"archive_size must be positive or -1");
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
	population::size_type max_archive_size, archive_size=0;

	if(m_archive_size == -1) {
		max_archive_size = NP;
	} else {
		max_archive_size = static_cast<population::size_type>(m_archive_size);
	}

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

	if (max_archive_size < 5 || (max_archive_size % 4 != 0) ) {
		pagmo_throw(value_error, "for SPEA2 at least 5 individuals in the archive are needed and the archive size must be a multiple of 4");
	}

	if (NP < 5 || (NP % 4 != 0) ) {
		pagmo_throw(value_error, "for SPEA2 at least 5 individuals in the population are needed and the population size must be a multiple of 4");
	}

	// Get out if there is nothing to do.
	if (NP == 0 || m_gen == 0) {
		return;
	}

	//population new_pop(pop);

	std::vector<spea2_individual> new_pop;
	for(unsigned int i=0; i < NP; ++i) {
		new_pop.push_back(spea2_individual());
		new_pop[i].x = pop.get_individual(i).cur_x;
		new_pop[i].f = pop.get_individual(i).cur_f;
		new_pop[i].c = pop.get_individual(i).cur_c;
	}

	//the cycle is until m_gen+1 because on the last iteration is necessary to do the truncation of the archive according to the dimension of the population to return pop
	for(int g = 0; g <= m_gen; ++g) {

		std::vector<double> F(NP+archive_size,0); 
		
		//computation of individuals' fitness (according to raw fitness and density)
		compute_spea2_fitness(F, sqrt(NP + archive_size), new_pop, prob);

		std::vector<population::size_type> ordered_by_fitness = pagmo::util::neighbourhood::order(F);

		unsigned int n_non_dominated;
		for(n_non_dominated = 0; n_non_dominated < NP+archive_size && F[ordered_by_fitness[n_non_dominated]] < 1; ++n_non_dominated);

		//if it is the last iteration store in pop the best solutions between the archive and the current population so the dimension of the archive is force to be the same as the population
		if(g==m_gen)
			max_archive_size = NP;

		if(n_non_dominated > max_archive_size) { //truncate according to delta
			archive_size = max_archive_size;

			while(n_non_dominated>max_archive_size){

				//fitness vector fo the non-dominated individuals
				std::vector<fitness_vector> fit_nd(n_non_dominated);
				for ( population::size_type i = 0; i<n_non_dominated; i++ ) {
					fit_nd[i]	=	new_pop[ordered_by_fitness[i]].f;
				}

				//computing K-NN distances for K=0,...,n_non_dominated
				std::vector<std::vector<pagmo::population::size_type> > neighbours_nd;
				pagmo::util::neighbourhood::euclidian::compute_neighbours(neighbours_nd, fit_nd);

				//sorting the individuals according to the density estimation given by the K-NN
				std::vector<population::size_type> ordered_by_delta = order_by_delta(neighbours_nd, fit_nd);

				//identifying the last individual of the ordered vector (the one to be eliminated)
				population::size_type idx = ordered_by_delta.size() - 1 ;
				population::size_type pop_idx = ordered_by_fitness[ordered_by_delta[idx]];

				new_pop.erase(new_pop.begin()+pop_idx);
				ordered_by_fitness.erase(ordered_by_fitness.begin()+ordered_by_delta[idx]);
				n_non_dominated--;

				//update the vector of indeces after the removal of the po_idx element
				for(population::size_type j=0; j<ordered_by_fitness.size(); j++){
					if(pop_idx<ordered_by_fitness[j])
						ordered_by_fitness[j]--;
				}

			}

			//erasing all the dominated solutions
			for(population::size_type i=archive_size; i<ordered_by_fitness.size(); i++){
				new_pop.erase(new_pop.begin() + ordered_by_fitness[i]);
				for(population::size_type j=0; j<ordered_by_fitness.size(); j++){
					if(ordered_by_fitness[i]<ordered_by_fitness[j])
						ordered_by_fitness[j]--;
				}
			}

		} else { //fill with dominated individuals
			archive_size = std::min(max_archive_size, ordered_by_fitness.size());

			//erasing all the dominated solutions
			for(population::size_type i=archive_size; i<ordered_by_fitness.size(); i++){
				new_pop.erase(new_pop.begin() + ordered_by_fitness[i]);
				for(population::size_type j=0; j<ordered_by_fitness.size(); j++){
					if(ordered_by_fitness[i]<ordered_by_fitness[j])
						ordered_by_fitness[j]--;
				}
			}
		}

		if(g==m_gen){
			//return the population stored into the archive that now has dimension NP
			pop.clear();
			for(unsigned int i=0; i<NP; ++i) {
				pop.push_back(new_pop[i].x);
			}
			break;//exiting from cycle over the generations 
		}

		//Now new_pop contains just the archive. We perform the genetic operations on the archive individuals and fill the population

		std::vector<population::size_type> shuffle1(archive_size), shuffle2(archive_size);
		for (pagmo::population::size_type i=0; i< archive_size; i++) {
			shuffle1[i] = i;
			shuffle2[i] = i;
		}
		boost::uniform_int<int> pop_idx(0,archive_size-1);
		boost::variate_generator<boost::mt19937 &, boost::uniform_int<int> > p_idx(m_urng,pop_idx);

		population::size_type parent1_idx, parent2_idx;
		decision_vector child1_x(D), child2_x(D);

		//We create some pseudo-random permutation of the poulation indexes
		std::random_shuffle(shuffle1.begin(),shuffle1.end(), p_idx);
		std::random_shuffle(shuffle2.begin(),shuffle2.end(), p_idx);

		//new_pop.update_pareto_information(); ///TODO: check why this was needed when using pagmo::population and if we still need

		std::vector<fitness_vector> new_pop_fit(n_non_dominated);
		std::vector<constraint_vector> new_pop_cons(n_non_dominated);
		for ( population::size_type i = 0; i<n_non_dominated; i++ ) {
			new_pop_fit[i]	=	new_pop[i].f;
			new_pop_cons[i]	=	new_pop[i].c;
		}

		std::vector<std::vector<population::size_type> > domination_list = compute_domination_list(prob, new_pop_fit,new_pop_cons);
		std::vector<population::size_type> pareto_rank = compute_pareto_rank(domination_list);

		int j = 0;
		//We then loop thorugh all individuals with increment 4 to select two pairs of parents that will
		//each create 2 new offspring
		for (pagmo::population::size_type i=0; i<archive_size; i+=4) {
			// We create two offsprings using the shuffled list 1
			parent1_idx = tournament_selection(shuffle1[i], shuffle1[i+1], pareto_rank);
			parent2_idx = tournament_selection(shuffle1[i+2], shuffle1[i+3], pareto_rank);
			crossover(child1_x, child2_x, parent1_idx,parent2_idx,new_pop, prob);
			mutate(child1_x,prob);
			mutate(child2_x,prob);

			spea2_individual child1;
			child1.x = child1_x;
			child1.f = prob.objfun(child1_x);
			child1.c = prob.compute_constraints(child1_x);

			spea2_individual child2;
			child2.x = child2_x;
			child2.f = prob.objfun(child2_x);
			child2.c = prob.compute_constraints(child2_x);

			new_pop.push_back(child1);
			new_pop.push_back(child2);
		}

		//Same with shuffle2
		for (pagmo::population::size_type i=0; i<archive_size; i+=4) {
			// We create two offsprings using the shuffled list 1
			parent1_idx = tournament_selection(shuffle2[i], shuffle2[i+1], pareto_rank);
			parent2_idx = tournament_selection(shuffle2[i+2], shuffle2[i+3], pareto_rank);
			crossover(child1_x, child2_x, parent1_idx,parent2_idx,new_pop, prob);
			mutate(child1_x,prob);
			mutate(child2_x,prob);

			spea2_individual child1;
			child1.x = child1_x;
			child1.f = prob.objfun(child1_x);
			child1.c = prob.compute_constraints(child1_x);

			spea2_individual child2;
			child2.x = child2_x;
			child2.f = prob.objfun(child2_x);
			child2.c = prob.compute_constraints(child2_x);

			new_pop.push_back(child1);
			new_pop.push_back(child2);
		}

	}

}

std::vector<std::vector<population::size_type> > spea2::compute_domination_list(const pagmo::problem::base &prob,
																   const std::vector<fitness_vector> &fit,
																   const std::vector<constraint_vector> &cons) const
{
	std::vector<population::size_type> dummy;
	std::vector<std::vector<population::size_type> > domination_list(fit.size(), dummy);

	for(unsigned int i=0; i<fit.size();++i) {
		for(unsigned int j=0; j<fit.size(); ++j) {
			// Check if individual in position i dominates individual in position n.
			if(prob.compare_fc(fit[i],cons[i],fit[j],cons[j])) {
				domination_list[i].push_back(j);
			}
		}
	}

	return domination_list;
}

std::vector<population::size_type> spea2::complement(std::vector<population::size_type> v, population::size_type N) const
{
	std::vector<population::size_type> rv(N - v.size());
	std::sort(v.begin(), v.end());
	population::size_type n = 0;
	std::vector<population::size_type>::size_type i = 0;
	std::vector<population::size_type>::size_type j = 0;
	while(i < v.size()) {
		while(n < v[i]) {
			rv[j] = n;
			j++;
			n++;
		}
		i++;
		n++;
	}
	while(n < N) {
		rv[j] = n;
		j++;
		n++;
	}
	return rv;
}

std::vector<population::size_type> spea2::order_by_delta(const std::vector<std::vector<pagmo::population::size_type>  > &neighbours, const std::vector<fitness_vector> &fit) const
{
	std::vector<population::size_type> rv(fit.size());
	for(unsigned int i=0; i<rv.size(); ++i) rv[i] = i;

	std::sort(rv.begin(), rv.end(), distance_sorter(neighbours, fit));

	return rv;
}

void spea2::compute_spea2_fitness(std::vector<double> &F,
			int K,
			const std::vector<spea2_individual> &pop,
			const pagmo::problem::base &prob) const
{

	const population::size_type NP = pop.size();
	std::vector<int> S(NP, 0);


	std::vector<fitness_vector> fit(NP);
	std::vector<fitness_vector> cons(NP);
	for ( population::size_type i = 0; i<NP; i++ ) {
		fit[i]	=	pop[i].f;
		cons[i]	=	pop[i].c;
	}
	std::vector<std::vector<pagmo::population::size_type> > neighbours;
	pagmo::util::neighbourhood::euclidian::compute_neighbours(neighbours, fit);

	/*for(unsigned i=0; i<NP; ++i) {
		domination_list.push_back(pop.get_domination_list(i));
	}*/

	std::vector<std::vector<population::size_type> > domination_list = compute_domination_list(prob, fit,cons);


	for(unsigned int i=0; i<NP; ++i) {
		S[i] = domination_list[i].size();
	}

	std::fill(F.begin(), F.end(), 0);

	for(unsigned int i=0; i<NP; ++i) {
		for(unsigned int j=0; j<domination_list[i].size(); ++j) {
			F[domination_list[i][j]] += S[i];
		}
	}

	for(unsigned int i=0; i<NP; ++i) {
		F[i] = F[i] +
				(1.0 / (pagmo::util::neighbourhood::euclidian::distance(fit[i], fit[neighbours[i][K]]) + 2));
	}
}

pagmo::population::size_type spea2::tournament_selection(pagmo::population::size_type idx1, pagmo::population::size_type idx2, const std::vector<population::size_type> &pareto_rank) const
{
	if (pareto_rank[idx1] < pareto_rank[idx2]) return idx1;
	if (pareto_rank[idx1] > pareto_rank[idx2]) return idx2;
	return ((m_drng() > 0.5) ? idx1 : idx2);
	
	/*fitness_vector f1 = pop.get_individual(idx1).cur_f;
	fitness_vector f2 = pop.get_individual(idx2).cur_f; 
	if (f1 < f2) {
		return idx1;
	}
	else if (f1 > f2) {
		return idx2;
	}
	else {
		return ((m_drng() > 0.5) ? idx1 : idx2);
	}*/
}

std::vector<population::size_type> spea2::compute_domination_count(const std::vector<std::vector<population::size_type> > &dom_list) const
{
	std::vector<population::size_type> domination_count(dom_list.size(),0);
	for(unsigned int i=0; i<dom_list.size(); ++i) {
		for(unsigned int j = 0; j < dom_list[i].size(); ++j) {
			domination_count[dom_list[i][j]]++;
		}
	}

	return domination_count;
}

std::vector<population::size_type> spea2::compute_pareto_rank(const std::vector<std::vector<population::size_type> > &dom_list) const
{
	std::vector<population::size_type> pareto_rank(dom_list.size(),0);

	// We define some utility vectors .....
	std::vector<population::size_type> F,S;

	// And make a copy of the domination count (number of individuals that dominating one individual)
	std::vector<population::size_type> dom_count = compute_domination_count(dom_list);

	// 1 - Find the first Pareto Front
	for (population::size_type idx = 0; idx < dom_count.size(); ++idx){
		if (dom_count[idx] == 0) {
			F.push_back(idx);
		}
	}

	unsigned int irank = 1;

	// We loop to find subsequent fronts
	while (F.size()!=0) {
		//For each individual F in the current front
		for (population::size_type i=0; i < F.size(); ++i) {
			//For each individual dominated by F
			for (population::size_type j=0; j<dom_list[F[i]].size(); ++j) {
				dom_count[dom_list[F[i]][j]]--;
				if (dom_count[dom_list[F[i]][j]] == 0){
					S.push_back(dom_list[F[i]][j]);
					pareto_rank[dom_list[F[i]][j]] = irank;
				}
			}
		}
		F = S;
		S.clear();
		irank++;
	}

	//std::cout << "pareto rank before return " << pareto_rank << std::endl;
	return pareto_rank;
}


void spea2::crossover(decision_vector& child1, decision_vector& child2, pagmo::population::size_type parent1_idx, pagmo::population::size_type parent2_idx,
					  const std::vector<spea2_individual> &pop, const pagmo::problem::base &prob) const
{

	problem::base::size_type D = prob.get_dimension();
	problem::base::size_type Di = prob.get_i_dimension();
	problem::base::size_type Dc = D - Di;
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const decision_vector& parent1 = pop[parent1_idx].x;
	const decision_vector& parent2 = pop[parent2_idx].x;
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

void spea2::mutate(decision_vector& child, const pagmo::problem::base& prob) const
{

	problem::base::size_type D = prob.get_dimension();
	problem::base::size_type Di = prob.get_i_dimension();
	problem::base::size_type Dc = D - Di;
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
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
