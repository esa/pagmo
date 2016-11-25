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


#include <string>
#include <vector>
#include <algorithm>
#include <fstream>

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
#include "nspso.h"


namespace pagmo { namespace algorithm {
/// Constructor
 /**
 * Constructs a NSPSO algorithm (multi objective PSO)
 *
 * @param[in] gen Number of generations to evolve.
 * @param[in] minW minimum particles' inertia weight (the inertia weight is decreased troughout the run between maxW and minW)
 * @param[in] maxW maximum particles' inertia weight (the inertia weight is decreased troughout the run between maxW and minW)
 * @param[in] C1 magnitude of the force, applied to the particle's velocity, in the direction of its previous best position
 * @param[in] C2 magnitude of the force, applied to the particle's velocity, in the direction of its global best (leader)
 * @param[in] CHI velocity scaling factor
 * @param[in] v_coeff velocity coefficient (determining the maximum allowed particle velocity)
 * @param[in] leader_selection_range the leader of each particle is selected among the best leader_selection_range% individuals
 * @param[in] diversity_mechanism the diversity mechanism to use to mantein diversity on the pareto front
 *
 * @throws value_error if gen is negative
 */
nspso::nspso(int gen, double minW, double maxW, double C1, double C2,
	  double CHI, double v_coeff, int leader_selection_range,
	   diversity_mechanism_type diversity_mechanism):base(),
	m_gen(gen),
	m_minW(minW),
	m_maxW(maxW),
	m_C1(C1),
	m_C2(C2),
	m_CHI(CHI),
	m_v_coeff(v_coeff),
	m_leader_selection_range(leader_selection_range),
	m_diversity_mechanism(diversity_mechanism)
{
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}

	if (minW <= 0 || maxW <= 0 || C1 <=0 || C2 <= 0 || CHI <=0) {
		pagmo_throw(value_error,"minW, maxW, C1, C2 and CHI should be greater than 0");
	}

	if (minW > maxW) {
		pagmo_throw(value_error,"minW should be lower than maxW");
	}

	if (v_coeff <= 0 || v_coeff > 1) {
		pagmo_throw(value_error,"v_coeff should be in the ]0,1] range");
	}

	if (leader_selection_range <=0 || leader_selection_range > 100) {
		pagmo_throw(value_error,"leader_selection_range should be in the ]0,100] range");
	}
}


/// Clone method.
base_ptr nspso::clone() const
{
	return base_ptr(new nspso(*this));
}

/// Evolve implementation.
/**
 * Run the NSPSO algorithm for the number of generations specified in the constructors.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */
void nspso::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base					&prob = pop.problem();
	const problem::base::size_type		D = prob.get_dimension(), prob_i_dimension = prob.get_i_dimension(), prob_c_dimension = prob.get_c_dimension(), prob_f_dimension = prob.get_f_dimension();
	const problem::base::size_type		Dc = D - prob_i_dimension;
	const decision_vector				&lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type			NP = pop.size();

	//We perform some checks to determine wether the problem/population are suitable for PSO
	if( Dc == 0 ){
		pagmo_throw(value_error,"There is no continuous part in the problem decision vector for NSPSO to optimise");
	}

	if( prob_c_dimension != 0 ){
		pagmo_throw(value_error,"The problem is not box constrained and NSPSO is not suitable to solve it");
	}

	if( prob_f_dimension < 2 ){
		pagmo_throw(value_error,"The problem is not multi-objective. Use a single-objectie optimization algorithm instead");
	}

	if(m_diversity_mechanism != CROWDING_DISTANCE && m_diversity_mechanism != NICHE_COUNT && m_diversity_mechanism != MAXMIN) {
		pagmo_throw(value_error,"non existing diversity mechanism method");
	}

	// Get out if there is nothing to do.
	if (NP == 0 || m_gen == 0) {
		return;
	}

	decision_vector minV(Dc), maxV(Dc);	// Maximum and minimum velocity allowed

	// Initialise the minimum and maximum velocity
	for(problem::base::size_type d = 0; d < Dc; d++ ){
		double v_width  = (ub[d] - lb[d]) * m_v_coeff;
		minV[d] = -1.0 * v_width;
		maxV[d] = v_width;
	}

	// 0 - Copy the population into nextPopList
	std::vector<nspso_individual> nextPopList;
	for(unsigned int i=0; i < NP; ++i) {
		nextPopList.push_back(nspso_individual());
		nextPopList[i].cur_x = pop.get_individual(i).cur_x;
		nextPopList[i].best_x = pop.get_individual(i).best_x;
		nextPopList[i].cur_v = pop.get_individual(i).cur_v;
		nextPopList[i].cur_f = pop.get_individual(i).cur_f;
		nextPopList[i].best_f = pop.get_individual(i).best_f;
		nextPopList[i].cur_c = pop.get_individual(i).cur_c;
		nextPopList[i].best_c = pop.get_individual(i).best_c;
	}

	for(int g = 0; g < m_gen; ++g) {

		std::vector<population::size_type> bestNonDomIndices;
		std::vector<fitness_vector> fit(NP);// particles' current fitness values
		std::vector<constraint_vector> cons(NP);
		for ( population::size_type i = 0; i<NP; i++ ) {
			fit[i]	=	nextPopList[i].best_f;
			cons[i] = nextPopList[i].best_c;
		}

		// 1 - Calculate non-dominated population
		if(m_diversity_mechanism == CROWDING_DISTANCE) {
			std::vector<std::vector<population::size_type> > dom_list = compute_domination_list(prob,fit,cons);
			std::vector<population::size_type> pareto_rank = compute_pareto_rank(dom_list);
			std::vector<std::vector<population::size_type> > pareto_fronts = compute_pareto_fronts(pareto_rank);
			std::vector<double> crowding_d = compute_crowding_d(fit, pareto_fronts);

			crowding_pareto_comp comp(pareto_rank, crowding_d);
			std::vector<population::size_type> dummy(NP);
			for(unsigned int i=0; i<NP; ++i) dummy[i] = i;
			std::sort(dummy.begin(), dummy.end(),comp);
			if(pareto_fronts[0].size() > 1) {
				bestNonDomIndices = std::vector<population::size_type>(dummy.begin(), dummy.begin() + pareto_fronts[0].size());
			} else { //ensure the non-dom population has at least 2 individuals (to avoid convergence to a point)
				bestNonDomIndices = std::vector<population::size_type>(dummy.begin(), dummy.begin() + 2);
			}

		} else if(m_diversity_mechanism == NICHE_COUNT) {
			std::vector<decision_vector> nonDomChromosomes;

			std::vector<std::vector<population::size_type> > dom_list = compute_domination_list(prob,fit,cons);
			std::vector<population::size_type> pareto_rank = compute_pareto_rank(dom_list);
			std::vector<std::vector<population::size_type> > pareto_fronts = compute_pareto_fronts(pareto_rank);

			for(unsigned int i = 0; i < pareto_fronts[0].size(); ++i) {
				nonDomChromosomes.push_back(nextPopList[pareto_fronts[0][i]].best_x);
			}

			std::vector<double> nadir = compute_nadir(fit, pareto_rank);
			std::vector<double> ideal = compute_ideal(fit, pareto_rank);

			//Fonseca-Fleming setting for delta
			double delta = 1.0;
			if (prob_f_dimension == 2) {
				delta = ( (nadir[0] - ideal[0]) + (nadir[1] - ideal[1]) ) / (nonDomChromosomes.size()-1);
			} else if (prob_f_dimension == 3) {
				const double d1 = nadir[0] - ideal[0];
				const double d2 = nadir[1] - ideal[1];
				const double d3 = nadir[2] - ideal[2];
				const double N = nonDomChromosomes.size();
				delta = sqrt(4*d2*d1*N + 4*d3*d1*N + 4*d2*d3*N + pow(d1,2) + pow(d2,2) + pow(d3,2)
							 - 2*d2*d1 - 2*d3*d1 - 2*d2*d3 + d1 + d2 + d3) / (2*(N-1));
			} else { //for higher dimension we just divide equally the entire volume containing the pareto front
				for(unsigned int i=0; i<nadir.size(); ++i) {
					delta *= nadir[i] - ideal[i];
				}
				delta = pow(delta, 1.0/nadir.size())/nonDomChromosomes.size();
			}

			std::vector<int> count(nonDomChromosomes.size(),0);
			compute_niche_count(count, nonDomChromosomes, delta);
			std::vector<population::size_type> tmp = pagmo::util::neighbourhood::order(count);

			if(pareto_fronts[0].size() > 1) {
				for(unsigned int i=0; i < tmp.size(); ++i) {
					bestNonDomIndices.push_back(pareto_fronts[0][tmp[i]]);
				}
			} else { //ensure the non-dom population has at least 2 individuals (to avoid convergence to a point)
				unsigned int minPopSize = 2;
				for(unsigned f = 0; minPopSize > 0 && f < pareto_fronts.size(); ++f) {
						for(unsigned int i=0; minPopSize > 0 && i < pareto_fronts[f].size(); ++i) {
							bestNonDomIndices.push_back(pareto_fronts[f][i]);
							minPopSize--;
						}
					}
			}
		} else { // m_diversity_method == MAXMIN
			std::vector<double> maxmin(NP,0);
			compute_maxmin(maxmin, fit);
			bestNonDomIndices = pagmo::util::neighbourhood::order(maxmin);
			unsigned int i = 1;
			for(;i < bestNonDomIndices.size() && maxmin[bestNonDomIndices[i]] < 0; ++i);
			if(i<2) i=2; //ensure the non-dom population has at least 2 individuals (to avoid convergence to a point)
			bestNonDomIndices = std::vector<population::size_type>(bestNonDomIndices.begin(), bestNonDomIndices.begin() + i); //keep just the non dominated
		}

		//decrease W from maxW to minW troughout the run
		const double W  = m_maxW - (m_maxW-m_minW)/m_gen * g;

		//2 - Move the points
		for(population::size_type idx = 0; idx < NP; ++idx) {
			const decision_vector &best_X = nextPopList[idx].best_x;
			const decision_vector &cur_X = nextPopList[idx].cur_x;
			const decision_vector &cur_V = nextPopList[idx].cur_v;

			// 2.1 - Calculate the leader
			int ext = ceil(bestNonDomIndices.size()*m_leader_selection_range/100.0)-1;
			if (ext < 1) ext = 1;
			unsigned int leaderIdx;
			do {
				leaderIdx = boost::uniform_int<int>(0,ext)(m_drng);
			} while (bestNonDomIndices[leaderIdx] == idx); //never pick yourself as a leader
			decision_vector leader = nextPopList[bestNonDomIndices[leaderIdx]].best_x;

			//Calculate some random factors
			const double r1 = boost::uniform_real<double>(0,1)(m_drng);
			const double r2 = boost::uniform_real<double>(0,1)(m_drng);

			//Calculate new velocity and new position for each particle
			decision_vector newX(Dc);
			decision_vector newV(Dc);
			for(decision_vector::size_type i = 0; i < cur_X.size(); ++i) {
				double v = W*cur_V[i] +
								m_C1*r1*(best_X[i] - cur_X[i]) +
								m_C2*r2*(leader[i] - cur_X[i]);

				if(v > maxV[i]){
					v = maxV[i];
				}
				else if(v < minV[i]) {
					v = minV[i];
				}

				double x = cur_X[i] + m_CHI*v;
				if(x > ub[i]) {
					x = ub[i];
					v = 0;
				} else if (x < lb[i]) {
					x = lb[i];
					v = 0;
				}

				newV[i] = v;
				newX[i] = x;
			}

			//Add the moved particle to the population
			nextPopList.push_back(nspso_individual());
			nextPopList[idx+NP].cur_x = newX;
			nextPopList[idx+NP].best_x = newX;
			nextPopList[idx+NP].cur_v = newV;
			nextPopList[idx+NP].cur_f = prob.objfun(newX);
			nextPopList[idx+NP].best_f = nextPopList[idx+NP].cur_f;
			nextPopList[idx+NP].cur_c = prob.compute_constraints(newX);
			nextPopList[idx+NP].best_c = nextPopList[idx+NP].cur_c;
		}

		//3 - Select the best NP individuals in the new population (of size 2*NP) according to pareto dominance
		std::vector<fitness_vector> nextPop_fit(nextPopList.size());
		std::vector<constraint_vector> nextPop_cons(nextPopList.size());
		for(unsigned int i=0; i<nextPopList.size(); ++i) {
			nextPop_fit[i] = nextPopList[i].best_f;
			nextPop_cons[i] = nextPopList[i].best_c;
		}

		std::vector<population::size_type> bestNextPopIndices(NP,0);

		if(m_diversity_mechanism != MAXMIN) {
			std::vector<std::vector<population::size_type> > nextPop_pareto_fronts = compute_pareto_fronts(prob, nextPop_fit, nextPop_cons);
			for(unsigned int f = 0, i=0; i<NP && f < nextPop_pareto_fronts.size(); ++f) {
				if(nextPop_pareto_fronts[f].size() < NP-i) { //then push the whole front in the population
					for(unsigned int j = 0; j < nextPop_pareto_fronts[f].size(); ++j) {
						bestNextPopIndices[i] = nextPop_pareto_fronts[f][j];
						++i;
					}
				} else {
					boost::uniform_int<int> pop_idx(0,nextPop_pareto_fronts[f].size());
					boost::variate_generator<boost::mt19937 &, boost::uniform_int<int> > p_idx(m_urng,pop_idx);
					std::random_shuffle(nextPop_pareto_fronts[f].begin(), nextPop_pareto_fronts[f].end(), p_idx);
					for(unsigned int j = 0; i<NP; ++j) {
						bestNextPopIndices[i] = nextPop_pareto_fronts[f][j];
						++i;
					}
				}
			}
		} else {
			std::vector<double> maxmin(2*NP,0);
			compute_maxmin(maxmin, nextPop_fit);
			std::vector<population::size_type> tmp = pagmo::util::neighbourhood::order(maxmin);
			bestNextPopIndices = std::vector<population::size_type>(tmp.begin(), tmp.begin() + NP);
		}

		// The nextPopList for the next generation will contain the best NP individuals out of 2NP according to pareto dominance
		std::vector<nspso_individual> tmpPop(NP);
		for(unsigned int i=0; i < NP; ++i) {
			tmpPop[i] = nextPopList[bestNextPopIndices[i]];
		}

		nextPopList = tmpPop;
	}

	//4 - Evolution is over, copy the last population back to the orginal population
	pop.clear();
	for(population::size_type i = 0; i < NP; ++i){
		pop.push_back(nextPopList[i].best_x);
		pop.set_x(i, nextPopList[i].cur_x);
		pop.set_v(i, nextPopList[i].cur_v);
	}

}

double nspso::minfit(unsigned int i, unsigned int j, const std::vector<fitness_vector> &fit) const
{
	double min = fit[i][0] - fit[j][0];
	for(unsigned int f=0; f<fit[i].size(); ++f) {
		double tmp = fit[i][f] - fit[j][f];
		if(tmp < min) {
			min = tmp;
		}
	}
	return min;
}

void nspso::compute_maxmin(std::vector<double> &maxmin, const std::vector<fitness_vector> &fit) const
{
	const unsigned int NP = fit.size();
	for(unsigned int i=0; i<NP; ++i) {
		maxmin[i] = minfit(i, (i+1)%NP, fit);
		for(unsigned j=0; j<NP; ++j) {
			if(i != j) {
				double tmp = minfit(i, j, fit);
				if(tmp > maxmin[i]) {
					maxmin[i] = tmp;
				}
			}
		}
	}
}



double nspso::euclidian_distance(const std::vector<double> &x, const std::vector<double> &y) const
{
	pagmo_assert(x.size() == y.size());
	double sum = 0;
	for(unsigned int i = 0; i < x.size(); ++i) {
		sum+= pow(x[i]-y[i],2);
	}
	return sqrt(sum);
}

void nspso::compute_niche_count(std::vector<int> &count, const std::vector<std::vector<double> > &chromosomes, double delta) const
{
	std::fill(count.begin(), count.end(),0);
	for(unsigned int i=0; i<chromosomes.size(); ++i) {
		for(unsigned int j=0; j<chromosomes.size(); ++j) {
			if(euclidian_distance(chromosomes[i], chromosomes[j]) < delta) {
				count[i]++;
			}
		}
	}

}

std::vector<std::vector<population::size_type> > nspso::compute_domination_list(const pagmo::problem::base &prob,
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

std::vector<population::size_type> nspso::compute_domination_count(const std::vector<std::vector<population::size_type> > &dom_list) const
{
	std::vector<population::size_type> domination_count(dom_list.size(),0);
	for(unsigned int i=0; i<dom_list.size(); ++i) {
		for(unsigned int j = 0; j < dom_list[i].size(); ++j) {
			domination_count[dom_list[i][j]]++;
		}
	}

	return domination_count;
}

std::vector<population::size_type> nspso::compute_pareto_rank(const std::vector<std::vector<population::size_type> > &dom_list) const
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

std::vector<double> nspso::compute_crowding_d(const std::vector<fitness_vector> &fit, const std::vector<std::vector<population::size_type> > &pareto_fronts) const {

	std::vector<double> crowding_d(fit.size(),0);

	for(unsigned f=0; f < pareto_fronts.size(); ++f) {
		std::vector<fitness_vector::size_type> I(pareto_fronts[f]);

		fitness_vector::size_type lastidx = I.size() - 1;

		// we construct the comparison functor along the first fitness component
		one_dim_fit_comp funct(fit,0);

		// we loop along fitness components
		for (fitness_vector::size_type i = 0; i < fit[0].size(); ++i) {
			funct.m_dim = i;
			// we sort I along the fitness_dimension i
			std::sort(I.begin(),I.end(), funct);
			// assign Inf to the boundaries
			crowding_d[I[0]] = std::numeric_limits<double>::max();
			crowding_d[I[lastidx]] = std::numeric_limits<double>::max();
			//and compute the crowding distance
			double df = fit[I[lastidx]][i] - fit[I[0]][i];
			for (population::size_type j = 1; j < lastidx; ++j) {
				if (df == 0.0) { 						// handles the case in which the pareto front collapses to one single point
					crowding_d[I[j]] += 0.0;			// avoiding creation of nans that can't be serialized
				} else {
					crowding_d[I[j]] += (fit[I[j+1]][i] -fit[I[j-1]][i])/df;
				}
			}
		}
	}

	return crowding_d;
}

std::vector<std::vector<population::size_type> > nspso::compute_pareto_fronts(const pagmo::problem::base &prob,
																			  const std::vector<fitness_vector> &fit,
																			  const std::vector<constraint_vector> &cons) const
{
	std::vector<std::vector<population::size_type> > dom_list = compute_domination_list(prob,fit,cons);
	std::vector<population::size_type> pareto_rank = compute_pareto_rank(dom_list);

	return compute_pareto_fronts(pareto_rank);
}

std::vector<std::vector<population::size_type> > nspso::compute_pareto_fronts(const std::vector<population::size_type> &pareto_rank) const
{
	std::vector<std::vector<population::size_type> > retval;

	for (population::size_type idx = 0; idx < pareto_rank.size(); ++idx) {
		if (pareto_rank[idx] >= retval.size()) {
			retval.resize(pareto_rank[idx] + 1);
		}
		retval[pareto_rank[idx]].push_back(idx);
	}

	return retval;
}

fitness_vector nspso::compute_ideal(const std::vector<fitness_vector> &fit, const std::vector<population::size_type> &pareto_rank) const {

	unsigned int firstFrontIdx = 0;
	for(;pareto_rank[firstFrontIdx] > 0; ++firstFrontIdx);

	fitness_vector ideal(fit[firstFrontIdx]);
	for (population::size_type idx = 0; idx < fit.size(); ++idx) {
		if (pareto_rank[idx] == 0) { //it is in the first pareto front
			for(fitness_vector::size_type i = 0; i < ideal.size(); ++i) {
				if (fit[idx][i] < ideal[i]) {
					ideal[i] = fit[idx][i];
				}
			}
		}
	}
	return ideal;
}

fitness_vector nspso::compute_nadir(const std::vector<fitness_vector> &fit, const std::vector<population::size_type> &pareto_rank) const {
	unsigned int firstFrontIdx = 0;
	for(;pareto_rank[firstFrontIdx] > 0; ++firstFrontIdx);

	fitness_vector nadir(fit[firstFrontIdx]);
		for (population::size_type idx = 0; idx < fit.size(); ++idx) {
		if (pareto_rank[idx] == 0) { //it is in the first pareto front
			for(fitness_vector::size_type i = 0; i < nadir.size(); ++i) {
				if (fit[idx][i] > nadir[i]) {
					nadir[i] = fit[idx][i];
				}
			}
		}
	}
	return nadir;
}


/// Algorithm name
std::string nspso::get_name() const
{
	return "Non-dominated Sorting Particle Swarm Optimizer (NSPSO)";
}

/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string nspso::human_readable_extra() const
{
	std::ostringstream s;
	s << "gen:" << m_gen << ' ';
	s << "minW:" << m_minW << ' ';
	s << "maxW:" << m_maxW << ' ';
	s << "C1:" << m_C1 << ' ';
	s << "C2:" << m_C2 << ' ';
	s << "CHI:" << m_CHI << ' ';
	s << "v_coeff:" << m_v_coeff << ' ';
	s << "leader_selection_range:" << m_leader_selection_range << ' ';
	s << "diversity_mechanism:";
	switch (m_diversity_mechanism)
	{
		case CROWDING_DISTANCE : s << "CROWDING_DISTANCE" << ' ';
			break;
		case NICHE_COUNT : s << "NICHE_COUNT" << ' ';
			break;
		case MAXMIN : s << "MAXMIN" << ' ';
			break;
	}
	return s.str();
}


}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::nspso)
