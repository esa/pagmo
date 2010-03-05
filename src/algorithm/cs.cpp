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

#include "cs.h"
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>


namespace pagmo { namespace algorithm {

cs::cs(const double &stop_range, const double &start_range = 0.25, const double &reduction_coeff = 0.5):base(),m_stop_range(stop_range),m_start_range(start_range),
		m_reduction_coeff(reduction_coeff)
{
	if (reduction_coeff >= 1 || reduction_coeff <=0) {
		pagmo_throw(value_error,"the reduction coefficient must be smaller than one and positive, You Fool!!");
	}
	if (start_range > 1 || start_range <= 0) {
		pagmo_throw(value_error,"the starting range must be smaller than one and positive, You Fool!!");
	}
	if (stop_range > 1 || stop_range <= 0 || stop_range>range_) {
		pagmo_throw(value_error,"the minimum range must be smaller than one, positive and smaller than the starting range, (o44portebat studuisse)!!");
	}
}

/// Clone method.
base_ptr cs::clone() const
{
	return base_ptr(new cs(*this));
}

/// Evolve implementation.
/**
 * Run the compass search algorithm for the number of iterations specified in the constructors.
 * At each accepted point velocity is also updated.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved. The population champion is considered.
 */

void cs::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type D = prob.get_dimension(), prob_i_dimension = prob.get_i_dimension(), prob_c_dimension = prob.get_c_dimension(), prob_f_dimension = prob.get_f_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP = pop.size();
	const problem::base::size_type Dc = D - prob_i_dimension;


	//We perform some checks to determine wether the problem/population are suitable for sa_corana
	if ( Dc == 0 ) {
		pagmo_throw(value_error,"There is no continuous part in the problem decision vector for sa_corana to optimise");
	}

	if ( prob_c_dimension != 0 ) {
		pagmo_throw(value_error,"The problem is not box constrained and sa_corana is not suitable to solve it");
	}

	if ( prob_f_dimension != 1 ) {
		pagmo_throw(value_error,"The problem is not single objective and sa_corana is not suitable to solve it");
	}

	//Determines the number of temperature adjustment for the annealing procedure
	const size_t n_T = m_niter / (m_step_adj * m_bin_size * Dc);
	if (n_T == 0) {
		pagmo_throw(value_error,"n_T is zero, increase niter");
	}

	// Get out if there is nothing to do.
	if (NP == 0 ) {
		return;
	}

	//Starting point is the best individual
	const int bestidx = pop.get_best_idx();
	const decision_vector &x0 = pop.get_individual(bestidx).cur_x;
	const fitness_vector &fit0 = pop.get_individual(bestidx).cur_f;

	decision_vector x=x0,newx;
	fitness_vector = f=fit0,newf;
	size_t D = problem.getDimension();
	bool flag = false;

	double newrange=range;

	while (newrange > minRange) {
		flag = false;
		for (unsigned int i=0; i<Dc; i++) {
			newx=x;
			newx[i] = x[i] + newrange * (ub[i]-lb[i]);
			//feasibility correction
			if (newx[i] > ub [i]) newx[i]=ub[i];

			prob.objfun(newf,newx);
			if (prob.compare_fitness(newf,f)) {
				f = newf;
				x = newx;
				flag=true;
				break; //accept
			}

			newx[i] = x[i] - newrange * (UB[i]-LB[i]);
			//feasibility correction
			if (newx[i] < LB [i]) newx[i]=LB[i];

			newf = problem.objfun(newx);
			if (newf < f) {  //accept
				f = newf;
				x = newx;
				flag=true;
				break;
			}
		}
		if (!flag) {
			newrange *= reduxCoeff;
		}
	} //end while
}




/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string sa_corana::human_readable_extra() const
{
	std::ostringstream s;
	s << "\tIteration:\t" << m_niter << '\n';
	s << "\tStarting Temperature:\t" << m_Ts << '\n';
	s << "\tFinal Temperature:\t" << m_Tf << '\n';
	s << "\tRatio of neighbourhood over temperature adjustments:\t" << m_step_adj << '\n';
	s << "\tSize of the bin to evaluate the acceptance rate:\t" << m_bin_size << '\n';
	s << "\tSize of the strating range:\t" << m_range << '\n';
	return s.str();
}

}} //namespaces
