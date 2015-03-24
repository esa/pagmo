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

#include <cmath>
#include <iostream>
#include <boost/functional/hash.hpp>

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "robust.h"

using namespace std;

namespace pagmo { namespace problem {

/**
 * Construct by specifying a problem to be transformed and the parameter rho.
 *
 * @param[in] p pagmo::problem::base to be made robust
 * @param[in] trials sample size to average over in objective / constraint function
 * @param[in] param_rho neighbourhood size
 * @param[in] seed seed of the underlying rng
 *
 * @see problem::base_stochastic constructors.
 */

robust::robust(const base & p, unsigned int trials, const double param_rho, unsigned int seed):
	base_stochastic((int)p.get_dimension(),
		 p.get_i_dimension(),
		 p.get_f_dimension(),
		 p.get_c_dimension(),
		 p.get_ic_dimension(),
		 p.get_c_tol(), seed),
	m_original_problem(p.clone()),
	m_normal_dist(0, 1),
	m_uniform_dist(0, 1),
	m_trials(trials),
	m_rho(param_rho)
{
	if(param_rho < 0){
		pagmo_throw(value_error, "Rho should be greater than 0");
	}
	set_bounds(p.get_lb(),p.get_ub());
}

/// Copy Constructor. Performs a deep copy
robust::robust(const robust &prob):
	 base_stochastic(prob),
	 m_original_problem(prob.m_original_problem->clone()),
	 m_uniform_dist(0, prob.m_rho),
	 m_trials(prob.m_trials),
	 m_rho(prob.m_rho) {}

/// Clone method.
base_ptr robust::clone() const
{
	return base_ptr(new robust(*this));
}

/**
 * Configure parameter to control the noise
 *
 * param[in] rho The parameter rho.
 *
 */
void robust::set_rho(double rho)
{
	m_rho = rho;
}

/**
 * Returns the parameter rho.
*/
double robust::get_rho() const
{
	return m_rho;
}

/// Implementation of the objective function.
/// Add noises to the decision vector before calling the actual objective function.
void robust::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	// Temporary storage used for averaging
	fitness_vector tmp(f.size(),0.0);
	f = tmp;

	// Set the seed
	m_drng.seed(m_seed);

	// Perturb decision vector and evaluate
	decision_vector x_perturbed(x);
	for(unsigned int i = 0; i < m_trials; ++i){
		inject_noise_x(x_perturbed);
		m_original_problem->objfun(tmp, x_perturbed);
		for(fitness_vector::size_type j = 0; j < f.size(); ++j){
			f[j] += tmp[j] / (double)m_trials;
		}
	}
}

/// Implementation of the constraints computation.
/// Add noises to the decision vector before calling the actual constraint function.
void robust::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	// Temporary storage used for averaging
	constraint_vector tmp(c.size(), 0.0);
	c = tmp;

	// Set the seed
	m_drng.seed(m_seed);

	// Perturb decision vector and evaluate
	decision_vector x_perturbed(x);
	for(unsigned int i = 0; i < m_trials; ++i){
		inject_noise_x(x_perturbed);
		m_original_problem->compute_constraints(tmp, x_perturbed);
		for(constraint_vector::size_type j = 0; j < c.size(); ++j){
			c[j] = tmp[j] / (double)m_trials;
		}
	}
}

/// Apply noise on the decision vector based on rho
void robust::inject_noise_x(decision_vector &x) const
{
	// We follow the algorithm at
	// http://math.stackexchange.com/questions/87230/picking-random-points-in-the-volume-of-sphere-with-uniform-probability

	// 0. Define the radius
	double radius = m_rho * pow(m_uniform_dist(m_drng),1.0/x.size());

	// 1. Sampling N(0,1) on each dimension
	std::vector<double> perturbation(x.size(), 0.0);
	double c2=0;
	for(size_type i = 0; i < perturbation.size(); i++){
		perturbation[i] = m_normal_dist(m_drng);
		c2 += perturbation[i]*perturbation[i];
	}

	// 2. Normalize the vector
	for(size_type i = 0; i < perturbation.size(); i++){
		perturbation[i] *= (radius / sqrt(c2) );
		x[i] += perturbation[i];
	}

	// 3. Clip the variables to the valid bounds
	for(base::size_type i = 0; i < x.size(); i++){
		x[i] = std::max(x[i], get_lb()[i]);
		x[i] = std::min(x[i], get_ub()[i]);
	}
}

std::string robust::get_name() const
{
	return m_original_problem->get_name() + " [Robust]";
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the translation vector
 */
std::string robust::human_readable_extra() const
{
	std::ostringstream oss;
	oss << m_original_problem->human_readable_extra() << std::endl;
	oss << "\tNeighbourhood radius = " << m_rho;
	oss << "\n\ttrials: "<<m_trials;
	oss << "\n\tseed: "<<m_seed<< std::endl;
	return oss.str();
}
}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::robust)
