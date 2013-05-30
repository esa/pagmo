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

#include <cmath>

#include "../exceptions.h"
#include "../types.h"
#include "../population.h"
#include "base.h"
#include "noisy.h"
#include <iostream>
using namespace std;

namespace pagmo { namespace problem {

/**
 * Construct by using the default noise ~N(0,0.1)
 *
 * @param[in] p pagmo::problem::base to be noisy
 *
 * @see problem::base_stochastic constructors.
 */

noisy::noisy(const base & p):
	base_stochastic((int)p.get_dimension(),
		 p.get_i_dimension(),
		 p.get_f_dimension(),
		 p.get_c_dimension(),
		 p.get_ic_dimension(),
		 p.get_c_tol(), 0),
	m_original_problem(p.clone()),
	m_normal_dist(0, 0.1),
	m_param_mu(0),
	m_param_sigma(0.1)
	
{
	set_bounds(p.get_lb(),p.get_ub());
}

/**
 * Construct by specifying a normally distributed noise
 *
 * @param[in] p pagmo::problem::base to be noisy
 * @param[in] mu Mean of the normal distribution underlying noise
 * @param[in] std Standard deviation of the normal distribution underlyig noise
 *
 * @see problem::base_stochastic constructors.
 */

noisy::noisy(const base & p, const double mu, const double sigma, unsigned int seed):
	base_stochastic((int)p.get_dimension(),
		 p.get_i_dimension(),
		 p.get_f_dimension(),
		 p.get_c_dimension(),
		 p.get_ic_dimension(),
		 p.get_c_tol(), seed),
	m_original_problem(p.clone()),
	m_normal_dist(mu, sigma),
	m_param_mu(mu),
	m_param_sigma(sigma)
	
{
	set_bounds(p.get_lb(),p.get_ub());
}

/// Copy Constructor. Performs a deep copy
noisy::noisy(const noisy &prob):
	base_stochastic((int)prob.get_dimension(), 
		 prob.get_i_dimension(),
		 prob.get_f_dimension(),
		 prob.get_c_dimension(),
		 prob.get_ic_dimension(),
		 prob.get_c_tol(),
		 prob.m_seed),
	m_original_problem(prob.m_original_problem->clone()),
	m_normal_dist(prob.m_normal_dist),
	m_param_mu(prob.m_param_mu),
	m_param_sigma(prob.m_param_sigma)
		 
{
	set_bounds(prob.get_lb(),prob.get_ub());
}

/// Clone method.
base_ptr noisy::clone() const
{
	return base_ptr(new noisy(*this));
}

/**
 * Configure parameters for the noise distribution
 * 
 * param[in] mu Mean of the noise
 * param[in] sigma Standard deviation of the noise
 *
 */ 
void noisy::set_noise_param(double mu, double sigma)
{
	boost::random::normal_distribution<double>::param_type params(mu, sigma);
	m_normal_dist.param(params);
}

/**
 * Returns the noise parameter: Mean
*/
double noisy::get_param_mu()
{
	return m_param_mu;
}

/**
 * Return the noise parameter: Standard deviation
 */
double noisy::get_param_sigma()
{
	return m_param_sigma;
}

/// Implementation of the objective function.
/// Add noises to the computed fitness vector.
void noisy::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	m_original_problem->objfun(f, x);
	inject_noise_f(f);
}

/// Implementation of the constraints computation.
/// Add noises to the computed constraint vector.
void noisy::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	m_original_problem->compute_constraints(c, x);
	inject_noise_c(c);
}

/// Apply noise on a fitness vector
void noisy::inject_noise_f(fitness_vector& f) const
{
	for(f_size_type i = 0; i < f.size(); i++){
		f[i] += m_normal_dist(m_drng);
	}
}

/// Apply noise on a constraint vector
void noisy::inject_noise_c(constraint_vector& c) const
{
	for(c_size_type i = 0; i < c.size(); i++){
		c[i] += m_normal_dist(m_drng);
	}
}

std::string noisy::get_name() const
{
	return m_original_problem->get_name() + " [Noisy]"; 
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the translation vector
 */
std::string noisy::human_readable_extra() const
{
	std::ostringstream oss;
	oss << m_original_problem->human_readable_extra() << std::endl;
	oss << "\n\tNoise type: ";
	oss << "\n\tNormal distribution of ("<<m_param_mu<<","<<m_param_sigma<<")";
	//oss << "\n\tRNG state: "<<m_drng;
	//oss << "\n\tDistribution state: "<<m_normal_dist;
	return oss.str();
}
}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::noisy);
