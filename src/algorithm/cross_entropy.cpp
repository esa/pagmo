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

#include <string>
#include <vector>
#include <boost/random/normal_distribution.hpp>
#include <boost/numeric/conversion/cast.hpp>


#include "cross_entropy.h"
#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../types.h"




namespace pagmo { namespace algorithm {

//This is a comparator used in the sort function to sort individuals by the number of other individuals dominated
class CompareFitness: std::binary_function<std::pair<population::individual_type,int> , std::pair<population::individual_type,int> , bool>
{
	population *pop;
	
	public: 
		CompareFitness(population &p) {
			pop = &p;
		}

		bool operator()(const std::pair<population::individual_type,int> &a, const std::pair<population::individual_type,int> &b) const {
			return (pop->n_dominated(a.first) > pop->n_dominated(b.first));
		}
};

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] iter number of iterations
 * @param[in] fraction_elite the fraction of samples to be considered elite
 * @param[in] alpha mean smothing factor
 * @param[in] beta standard deviation smothing factor
 * 
 * */
cross_entropy::cross_entropy(int iter, double fraction_elite, double alpha, double beta):base(),m_iter(iter),m_fraction_elite(fraction_elite),m_alpha(alpha),m_beta(beta){
	//prob = NULL;
}

/// Clone method.
base_ptr cross_entropy::clone() const
{
	return base_ptr(new cross_entropy(*this));
}

/// Evolve implementation.
/**
 * Run Cross Entropy
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void cross_entropy::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type prob_i_dimension = prob.get_i_dimension(), D = prob.get_dimension(), Dc = D - prob_i_dimension, prob_c_dimension = prob.get_c_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP = pop.size();
	const population::size_type Nelite = boost::numeric_cast<population::size_type>(ceil(m_fraction_elite * NP));

	//We perform some checks to determine whether the problem/population are suitable for Cross Entropy
	if ( Dc == 0 ) {
		pagmo_throw(value_error,"There is no continuous part in the problem decision vector for CE to optimise");
	}

	if ( prob.get_f_dimension() != 1 ) {
		pagmo_throw(value_error,"The problem is not single objective and CE is not suitable to solve it");
	}

	if ( prob_c_dimension != 0 ) {
		pagmo_throw(value_error,"The problem is not box constrained and CE is not suitable to solve it");
	}

	if (NP < 2) {
		pagmo_throw(value_error,"for CE at least 2 individuals in the population are needed");
	}

	// Get out if there is nothing to do.
	if (m_iter == 0) {
		return;
	}

	// Some vectors used during evolution are allocated here.
	decision_vector dummy(D,0);			    				    //used for initialisation purposes
	std::vector<decision_vector> X(NP,dummy);	    				    //set of solutions
	std::vector<decision_vector> temp_X(Nelite,dummy);				    //set of elite solutions
	std::vector<std::pair<population::individual_type,int> > individuals(NP); 	    //set of individuals
	decision_vector population_mean(Dc,0);   					    //population's mean
	decision_vector tmp_population_mean(Dc,0);   	    				    //tmp population's mean
	decision_vector population_std(Dc,0);               				    //population's standard deviation
	decision_vector tmp_population_std(Dc,0);               			    //tmp population's standard deviation

	// Get population
	for ( population::size_type i = 0; i<NP; i++ ) {
		X[i]	   	       =	pop.get_individual(i).cur_x;
		individuals[i].first   =	pop.get_individual(i);
		individuals[i].second  =	i;
	}

	population_mean = calculate_mean(X);
	population_std = calculate_std(X, population_mean);
	
	//Main CE loop
	for(int t = 0; t < m_iter; ++t) {
		for(population::size_type i = 0; i < NP; ++i) {
			for(problem::base::size_type k = 0; k < Dc; ++k) {
				X[i][k] = boost::normal_distribution<double>(population_mean[k], population_std[k])(m_drng); //random gaussian sample

				//check constraints
				if (X[i][k] < lb[k]) {
					X[i][k] = lb[k];
				}			
				if (X[i][k] > ub[k]) {
					X[i][k] = ub[k];
				}
			}
			pop.set_x(i,X[i]);
		}
		//sort the individuals vector
		CompareFitness comp(pop);
		std::sort(individuals.begin(), individuals.end(), comp);

		//get the best Nelite individuals according to the number of dominated individuals
		for(population::size_type i = 0; i < Nelite; ++i) {
			temp_X[i] = X[individuals[i].second];
		}

		//calculate mean and std on these best individuals
		tmp_population_mean = calculate_mean(temp_X);
		tmp_population_std = calculate_std(temp_X, tmp_population_mean);
		
		double beta = m_beta - m_beta * pow((1.0 - 1.0/(t+1)),7);
		
		//smooth mean and standard deviation with old ones
		for(problem::base::size_type k=0; k < Dc; ++k) {
			population_mean[k] = m_alpha * tmp_population_mean[k] + (1-m_alpha) * population_mean[k];
			population_std[k] = beta * tmp_population_std[k] + (1-beta) * population_std[k];
		}
	}
	
}

//Calculate the mean vector of a vector calculating the mean of each component
decision_vector cross_entropy::calculate_mean(std::vector<decision_vector> X) {
	decision_vector mean_vector(X[0].size(), 0);

	for(decision_vector::size_type k = 0; k < X[0].size(); ++k) {
		for(std::vector<decision_vector>::size_type i = 0; i < X.size(); ++i) {
			mean_vector[k] += X[i][k];
		}
		mean_vector[k] /= X.size();
	}
	return mean_vector;
}

//Calculate the standard deviation vector of a vector calculating the standard deviation of each component
decision_vector cross_entropy::calculate_std(std::vector<decision_vector> X, decision_vector mean_vector) {
	decision_vector std_vector(X[0].size(), 0);
	for(decision_vector::size_type k = 0; k < X[0].size(); ++k) {
		for(std::vector<decision_vector>::size_type i = 0; i < X.size(); ++i) {
			std_vector[k] += (X[i][k] - mean_vector[k]) * (X[i][k] - mean_vector[k]);
		}
		std_vector[k] /= X.size();
		std_vector[k]  = sqrt(std_vector[k]);
	}
	return std_vector;
}

/// Algorithm name
std::string cross_entropy::get_name() const
{
	return "Cross Entropy method";
}


/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string cross_entropy::human_readable_extra() const
{
	std::ostringstream s;
	s << "iter:" << m_iter << ' ';
	s << "fraction_elite:" << m_fraction_elite << ' ';
	return s.str();
}



}} //namespaces
