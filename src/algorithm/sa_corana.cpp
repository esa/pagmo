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

#include "sa_corana.h"
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>


namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] niter number of total iterations.
 * @param[in] Ts starting temperature
 * @param[in] Tf final temperature
 * @param[in] niterT
 * @param[in] niterR
 * @param[in] range
 * @throws value_error niter is non positive, Ts is greater than Tf, Ts is non positive, Tf is non positive,
 * niterT or niterR are negative, range is not in the [0,1] interval
 */
sa_corana::sa_corana(int niter, const double &Ts, const double &Tf, int niterT, int niterR, const double &range):
		base(),m_niter(niter),m_Ts(Ts),m_Tf(Tf),m_step_adj(niterT),m_bin_size(niterR),m_range(range)
{
	if (niter < 0) {
		pagmo_throw(value_error,"number of iterations must be nonnegative");
	}
	if (Ts <= 0 || Tf <= 0 || Ts <= Tf) {
		pagmo_throw(value_error,"temperatures must be positive and Ts must be greater than Tf");
	}
	if (niterT < 0) {
		pagmo_throw(value_error,"number of iteration before adjusting the temperature must be positive");
	}
	if (niterR < 0) {
		pagmo_throw(value_error,"number of iteration before adjusting the neighbourhood must be positive");
	}
	if (range < 0 || range >1) {
		pagmo_throw(value_error,"Initial range must be between 0 and 1");
	}
}
/// Clone method.
base_ptr sa_corana::clone() const
{
	return base_ptr(new sa_corana(*this));
}

/// Evolve implementation.
/**
 * Run the sa_corana algorithm for the number of iterations specified in the constructors.
 * At each accepted point velocity is also updated.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved. Best member only is evolved.
 * Velocity is evaluated at the end as difference between decision vector before and after evolution
 */

void sa_corana::evolve(population &pop) const {

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

	// Get out if there is nothing to do.
	if (NP == 0 || m_niter == 0) {
		return;
	}
	if (n_T == 0) {
		pagmo_throw(value_error,"n_T is zero, increase niter");
	}

	//Starting point is the best individual
	const int bestidx = pop.get_best_idx();
	const decision_vector &x0 = pop.get_individual(bestidx).cur_x;
	const fitness_vector &fit0 = pop.get_individual(bestidx).cur_f;
	//Determines the coefficient to dcrease the temperature
	const double Tcoeff = std::pow(m_Tf/m_Ts,1.0/(double)(n_T));
	//Stores the current and new points
	decision_vector xNEW = x0, xOLD = xNEW;
	fitness_vector fNEW = fit0, fOLD = fNEW;
	//Stores the adaptive steps of each component (integer part included but not used)
	decision_vector step(D,m_range);

	//Stores the number of accepted points per component (integer part included but not used)
	std::vector<int> acp(D,0) ;
	double ratio = 0, currentT = m_Ts, probab = 0;

	//Main SA loops
	for (size_t jter = 0; jter < n_T; ++jter) {
		for (int mter = 0; mter < m_step_adj; ++mter) {
			for (int kter = 0; kter < m_bin_size; ++kter) {
				size_t nter = boost::uniform_int<int>(0,Dc-1)(m_urng);
				for (size_t numb = 0; numb < Dc ; ++numb) {
					nter = (nter + 1) % Dc;
					//We modify the current point actsol by mutating its nter component within
					//a step that we will later adapt
					xNEW[nter] = xOLD[nter] + boost::uniform_real<double>(-1,1)(m_drng) * step[nter] * (ub[nter]-lb[nter]);

					// If new solution produced is infeasible ignore it
					if ((xNEW[nter] > ub[nter]) || (xNEW[nter] < lb[nter])) {
						xNEW[nter]=xOLD[nter];
						continue;
					}
					//And we valuate the objective function for the new point
					prob.objfun(fNEW,xNEW);

					// We decide wether to accept or discard the point
					if (prob.compare_fitness(fNEW,fOLD) ) {
						//accept
						xOLD[nter] = xNEW[nter];
						fOLD = fNEW;
						acp[nter]++;	//Increase the number of accepted values
					} else {
						//test it with Boltzmann to decide the acceptance
						probab = exp ( - fabs(fOLD[0] - fNEW[0] ) / currentT );

						// we compare prob with a random probability.
						if (probab > m_drng()) {
							xOLD[nter] = xNEW[nter];
							fOLD = fNEW;
							acp[nter]++;	//Increase the number of accepted values
						} else {
							xNEW[nter] = xOLD[nter];
						}
					} // end if
				} // end for(nter = 0; ...
			} // end for(kter = 0; ...
			// adjust the step (adaptively)
			for (size_t iter = 0; iter < Dc; ++iter) {
				ratio = (double)acp[iter]/(double)m_bin_size;
				acp[iter] = 0;  //reset the counter
				if (ratio > .6) {
					//too many acceptances, increase the step by a factor 3 maximum
					step[iter] = step [iter] * (1 + 2 *(ratio - .6)/.4);
				} else {
					if (ratio < .4) {
						//too few acceptance, decrease the step by a factor 3 maximum
						step [iter]= step [iter] / (1 + 2 * ((.4 - ratio)/.4));
					};
				};
				//And if it becomes too large, reset it to its initial value
				if ( step[iter] > m_range ) {
					step [iter] = m_range;
				};
			}
		}
		// Cooling schedule
		currentT *= Tcoeff;
	}
	if ( prob.compare_fitness(fOLD,fit0) ){
		pop.set_x(bestidx,xOLD); //new evaluation is possible here......
		std::transform(xOLD.begin(), xOLD.end(), pop.get_individual(bestidx).cur_x.begin(), xOLD.begin(),std::minus<double>());
		pop.set_v(bestidx,xOLD);
	}
}


/// Algorithm name
std::string sa_corana::get_name() const
{
	return "Simulated Annealing (Corana's)";
}

/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string sa_corana::human_readable_extra() const
{
	std::ostringstream s;
	s << "iter:" << m_niter << ' ';
	s << "Ts:" << m_Ts << ' ';
	s << "Tf:" << m_Tf << ' ';
	s << "steps:" << m_step_adj << ' ';
	s << "bin_size:" << m_bin_size << ' ';
	s << "range:" << m_range << ' ';
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::sa_corana)
