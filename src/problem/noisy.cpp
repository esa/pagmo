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
 * Construct by specifying a problem to be transformed and the noise distribution,
 * controlled by a flag and two params. Currently two types of noise distribution is
 * supported, namely the normally distributed noise (NORMAL) and uniformly distributed 
 * noise (UNIFORM).
 *
 * @param[in] p pagmo::problem::base to be noisy
 * @param[in] param_first Mean of the Gaussian noise / Lower bound of the uniform noise
 * @param[in] param_second Standard deviation of the Gaussian noise / Upper bound of the uniform noise
 * @param[in] noise_type Two types of noise is currently supported: NORMAL or UNIFORM
 * @param[in] seed seed for the underlying rng
 * @see problem::base_stochastic constructors.
 */

noisy::noisy(const base & p, const double param_first, const double param_second, noise_distribution::type noise_type, unsigned int seed):
	base_stochastic((int)p.get_dimension(),
		 p.get_i_dimension(),
		 p.get_f_dimension(),
		 p.get_c_dimension(),
		 p.get_ic_dimension(),
		 p.get_c_tol(), seed),
	m_original_problem(p.clone()),
	m_normal_dist(param_first, param_second),
	m_uniform_dist(param_first, param_second),
	m_param_first(param_first),
	m_param_second(param_second),
	m_noise_type(noise_type)
	
{
	if(noise_type == noise_distribution::UNIFORM && param_first > param_second){
		pagmo_throw(value_error, "Bounds specified for the uniform noise is not valid.");
	}
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
	m_param_first(prob.m_param_first),
	m_param_second(prob.m_param_second)
		 
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
 * param[in] param_first Mean of the Gaussian noise / Lower bound of the uniform noise
 * param[in] param_second Standard deviation of the Gaussian noise / Upper bound of the uniform noise
 *
 */ 
void noisy::set_noise_param(double param_first, double param_second)
{
	if(m_noise_type == noise_distribution::NORMAL){
		boost::random::normal_distribution<double>::param_type params(param_first, param_second);
		m_normal_dist.param(params);
	}
	else if(m_noise_type == noise_distribution::UNIFORM){
		boost::random::uniform_real_distribution<double>::param_type params(param_first, param_second);
		m_uniform_dist.param(params);
	}
}

/**
 * Returns the noise parameter: Mean
*/
double noisy::get_param_first()
{
	return m_param_first;
}

/**
 * Return the noise parameter: Standard deviation
 */
double noisy::get_param_second()
{
	return m_param_second;
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
		if(m_noise_type == noise_distribution::NORMAL){
			f[i] += m_normal_dist(m_drng);
		}
		else if(m_noise_type == noise_distribution::UNIFORM){
			f[i] += m_uniform_dist(m_drng);
		}
	}
}

/// Apply noise on a constraint vector
void noisy::inject_noise_c(constraint_vector& c) const
{
	for(c_size_type i = 0; i < c.size(); i++){
		if(m_noise_type == noise_distribution::NORMAL){
			c[i] += m_normal_dist(m_drng);
		}
		else if(m_noise_type == noise_distribution::UNIFORM){
			c[i] += m_uniform_dist(m_drng);
		}
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
	if(m_noise_type == noise_distribution::NORMAL){
		oss << "\n\tNormal distribution of ("<<m_param_first<<","<<m_param_second<<")";
	}
	else if(m_noise_type == noise_distribution::UNIFORM){
		oss << "\n\tUniform distribution over ("<<m_param_first<<","<<m_param_second<<")";
	}
	else{
		oss << "\n\t Unknown????";
	}
	//oss << "\n\tRNG state: "<<m_drng;
	//oss << "\n\tDistribution state: "<<m_normal_dist;
	return oss.str();
}
}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::noisy);
