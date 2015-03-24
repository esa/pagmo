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
#include "../population.h"
#include "base.h"
#include "noisy.h"

using namespace std;

namespace pagmo { namespace problem {

/**
 * Construct by specifying a problem to be transformed and the noise distribution,
 * controlled by a flag and two params. Currently two types of noise distribution are
 * supported, namely the normally distributed noise (NORMAL) and uniformly distributed 
 * noise (UNIFORM).
 *
 * @param[in] p pagmo::problem::base to be noisy
 * @param[in] trials number of samples to average upon
 * @param[in] param_first Mean of the Gaussian noise / Lower bound of the uniform noise
 * @param[in] param_second Standard deviation of the Gaussian noise / Upper bound of the uniform noise
 * @param[in] distribution Two types of distributions for the noise are currently supported: NORMAL or UNIFORM
 * @param[in] seed seed for the underlying rng
 * @see problem::base_stochastic constructors.
 */

noisy::noisy(const base & p, unsigned int trials, const double param_first, const double param_second, noise_type distribution, unsigned int seed):
	base_stochastic((int)p.get_dimension(),
		 p.get_i_dimension(),
		 p.get_f_dimension(),
		 p.get_c_dimension(),
		 p.get_ic_dimension(),
		 p.get_c_tol(), seed),
	m_original_problem(p.clone()),
	m_trials(trials),
	m_normal_dist(0.0,1.0),
	m_uniform_dist(0.0,1.0),
	m_decision_vector_hash(),
	m_param_first(param_first),
	m_param_second(param_second),
	m_noise_type(distribution)
{
	if(distribution == UNIFORM && param_first > param_second){
		pagmo_throw(value_error, "Bounds specified for the uniform noise are not valid.");
	}
	set_bounds(p.get_lb(),p.get_ub());
}

/// Copy Constructor. Performs a deep copy
noisy::noisy(const noisy &prob):
	base_stochastic(prob),
	m_original_problem(prob.m_original_problem->clone()),
	m_trials(prob.m_trials),
	m_normal_dist(0.0,1.0),
	m_uniform_dist(0.0,1.0),
	m_decision_vector_hash(),
	m_param_first(prob.m_param_first),
	m_param_second(prob.m_param_second),
	m_noise_type(prob.m_noise_type) {}

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
	if(m_noise_type == UNIFORM && param_first > param_second){
		pagmo_throw(value_error, "Bounds specified for the uniform noise are not valid.");
	}
	m_param_first = param_first;
	m_param_second = param_second;
}

/**
 * Returns the first parameter.
 * Interpretation depends on the noise specified, see constructor.
*/
double noisy::get_param_first() const
{
	return m_param_first;
}

/**
 * Return the second parameter.
 * Interpretation depends on the noise specified, see constructor.
 */
double noisy::get_param_second() const
{
	return m_param_second;
}

/// Implementation of the objective function.
/// Add noises to the computed fitness vector.
void noisy::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	//1 - Initialize a temporary fitness vector storing one trial result
	//and we use it also to init the return value 
	fitness_vector tmp(f.size(),0.0);
	f=tmp;
	//2 - We set the seed
	m_drng.seed(m_seed+m_decision_vector_hash(x));
	//3 - We average upon multiple runs
	for (unsigned int j=0; j< m_trials; ++j) {
		m_original_problem->objfun(tmp, x);
		inject_noise_f(tmp);
		for (fitness_vector::size_type i=0; i<f.size();++i) {
			f[i] = f[i] + tmp[i] / (double)m_trials;
		}
	}
}

/// Implementation of the constraints computation.
/// Add noises to the computed constraint vector.
void noisy::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	//1 - Initialize a temporary constraint vector storing one trial result
	//and we use it also to init the return value 
	constraint_vector tmp(c.size(),0.0);
	c=tmp;
	//2 - We set the seed
	m_drng.seed(m_seed+m_decision_vector_hash(x));
	//3 - We average upon multiple runs
	for (unsigned int j=0; j< m_trials; ++j) {
		m_original_problem->compute_constraints(tmp, x);
		inject_noise_c(tmp);
		for (constraint_vector::size_type i=0; i<c.size();++i) {
			c[i] = c[i] + tmp[i] / (double)m_trials;
		}
	}
}

/// Apply noise on a fitness vector
void noisy::inject_noise_f(fitness_vector& f) const
{
	for(f_size_type i = 0; i < f.size(); i++){
		if(m_noise_type == NORMAL){
			f[i] += m_normal_dist(m_drng)*m_param_second+m_param_first;
		}
		else if(m_noise_type == UNIFORM){
			f[i] += m_uniform_dist(m_drng)*(m_param_second-m_param_first)+m_param_first;
		}
	}
}

/// Apply noise on a constraint vector
void noisy::inject_noise_c(constraint_vector& c) const
{
	for(c_size_type i = 0; i < c.size(); i++){
		if(m_noise_type == NORMAL){
			c[i] += m_normal_dist(m_drng)*m_param_second+m_param_first;
		}
		else if(m_noise_type == UNIFORM){
			c[i] += m_uniform_dist(m_drng)*(m_param_second-m_param_first)+m_param_first;
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
	oss << "\tNoise type: ";
	if(m_noise_type == NORMAL){
		oss << "Normal distribution of ("<<m_param_first<<","<<m_param_second<<")";
	}
	else if(m_noise_type == UNIFORM){
		oss << "Uniform distribution over ("<<m_param_first<<","<<m_param_second<<")";
	}
	else{
		oss << "\n\t Unknown????";
	}
	oss << "\n\ttrials: "<<m_trials;
	oss << "\n\tseed: "<<m_seed << std::endl;
	//oss << "\n\tDistribution state: "<<m_normal_dist;
	return oss.str();
}
}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::noisy)
