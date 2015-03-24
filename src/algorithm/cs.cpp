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

#include "cs.h"
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>


namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] max_eval Maximum number of function evaluations. The actual number might be much lower.
 * @param[in] stop_range Stopping criteron based on the perturbation size
 * @param[in] start_range Starting perturbation size
 * @param[in] reduction_coeff Size reduction of the perturbation size
 * @throws value_error if start and stop range not \f$ \in [0,1[ \f$ and not decreasing. max_eval negative
 * reduction_coeff not \f$ \in ]0,1[\f$
 */

cs::cs(const int& max_eval, const double &stop_range, const double &start_range, const double &reduction_coeff)
	:base(),m_stop_range(stop_range),m_start_range(start_range),m_reduction_coeff(reduction_coeff),m_max_eval(max_eval)
{
	if (reduction_coeff >= 1 || reduction_coeff <=0) {
		pagmo_throw(value_error,"the reduction coefficient must be smaller than one and positive, You Fool!!");
	}
	if (start_range > 1 || start_range <= 0) {
		pagmo_throw(value_error,"the starting range must be smaller than one and positive, You Fool!!");
	}
	if (stop_range > 1 || stop_range <= 0 || stop_range>start_range) {
		pagmo_throw(value_error,"the minimum range must be smaller than one, positive and smaller than the starting range, (o44portebat studuisse)!!");
	}
	if (max_eval < 0) {
		pagmo_throw(value_error,"Maximum number of function evaluations needs to be positive");
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
 * @param[in,out] pop input/output pagmo::population to be evolved. Best member only is evolved.
 * Velocity is evaluated at the end as difference between decision vector before and after evolution
 */

void cs::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type D = prob.get_dimension(), prob_i_dimension = prob.get_i_dimension(), prob_c_dimension = prob.get_c_dimension(), prob_f_dimension = prob.get_f_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP = pop.size();
	const problem::base::size_type Dc = D - prob_i_dimension;


	//We perform some checks to determine whether the problem/population are suitable for compass search
	if ( Dc == 0 ) {
		pagmo_throw(value_error,"There is no continuous part in the problem decision vector for compass search to optimise");
	}

	if ( prob_c_dimension != 0 ) {
		pagmo_throw(value_error,"The problem is not box constrained and compass search is not suitable to solve it");
	}

	if ( prob_f_dimension != 1 ) {
		pagmo_throw(value_error,"The problem is not single objective and compass search is not suitable to solve it");
	}

	// Get out if there is nothing to do.
	if (NP == 0 || m_max_eval == 0) {
		return;
	}

	//Starting point is the best individual
	const int bestidx = pop.get_best_idx();
	const decision_vector &x0 = pop.get_individual(bestidx).cur_x;
	const fitness_vector &fit0 = pop.get_individual(bestidx).cur_f;

	decision_vector x=x0,newx;
	fitness_vector f=fit0,newf=fit0;
	bool flag = false;
	int eval=0;

	double newrange=m_start_range;

	while (newrange > m_stop_range && eval <= m_max_eval) {
		flag = false;
		for (unsigned int i=0; i<Dc; i++) {
			newx=x;

			//move up
			newx[i] = x[i] + newrange * (ub[i]-lb[i]);
			//feasibility correction
			if (newx[i] > ub [i]) newx[i]=ub[i];

			prob.objfun(newf,newx); eval++;
			if (prob.compare_fitness(newf,f)) {
				f = newf;
				x = newx;
				flag=true;
				break; //accept
			}

			//move down
			newx[i] = x[i] - newrange * (ub[i]-lb[i]);
			//feasibility correction
			if (newx[i] < lb [i]) newx[i]=lb[i];

			prob.objfun(newf,newx); eval++;
			if (prob.compare_fitness(newf,f)) {  //accept
				f = newf;
				x = newx;
				flag=true;
				break;
			}
		}
		if (!flag) {
			newrange *= m_reduction_coeff;
		}
	} //end while
	std::transform(x.begin(), x.end(), pop.get_individual(bestidx).cur_x.begin(), newx.begin(),std::minus<double>()); // newx is now velocity
	pop.set_x(bestidx,x); //new evaluation is possible here......
	pop.set_v(bestidx,newx);
}

/// Algorithm name
std::string cs::get_name() const
{
	return "Compass Search";
}


/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string cs::human_readable_extra() const
{
	std::ostringstream s;
	s << "max_eval:" << m_max_eval << ' ';
	s << "stop_range:" << m_stop_range << ' ';
	s << "start_range:" << m_start_range << ' ';
	s << "reduction_coeff:" << m_reduction_coeff;
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::cs)
