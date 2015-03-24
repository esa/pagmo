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
#include <algorithm>

#include "../exceptions.h"
#include "../types.h"
#include "../population.h"
#include "death_penalty.h"

namespace pagmo { namespace problem {

/// Constructor from a constrained problem
/**
 * Will construct, using the selected method, the passed problem
 *
 * @param[in] problem base::problem to be modified to use a death penalty
 * as constraints handling technique.
 * @param[in] method death_penalty::method_type
 * @param[in] penalty_factors std::vector containing 
 *
 */
death_penalty::death_penalty(const base &problem, const method_type method, const std::vector<double> & penalty_factors):
	base_meta(
		 problem,
		 problem.get_dimension(),
		 problem.get_i_dimension(),
		 problem.get_f_dimension(),
		 0,
		 0,
		 std::vector<double>()),
	m_method(method),
	m_penalty_factors(penalty_factors)
{
	if(m_original_problem->get_c_dimension() <= 0){
		pagmo_throw(value_error,"The original problem is unconstrained.");
	}

	if( method > 2 || method < 0) {
		pagmo_throw(value_error, "the death penalty method must be set to 0 for simple death, 1 for Kuri death and 2 for static penalty with predefined penalty coefficients.");
	}

	if(method == WEIGHTED && m_penalty_factors.size() != m_original_problem->get_c_dimension()){
		pagmo_throw(value_error, "the vector of penalties factors is missing or needs to match constraints size");
	}
	
	if(method != WEIGHTED && m_penalty_factors.size() > 0){
		pagmo_throw(value_error, "Methods which are not WEIGHTED do not make use of weight vectors");
	}

}

/// Clone method.
base_ptr death_penalty::clone() const
{
	return base_ptr(new death_penalty(*this));
}

/// Implementation of the objective function.
/// (Wraps over the original implementation)
void death_penalty::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	constraint_vector c(m_original_problem->get_c_dimension(),0);
	m_original_problem->compute_constraints(c,x);

	if(m_original_problem->feasibility_c(c)) {
		m_original_problem->objfun(f, x);
	} else {
		double high_value = boost::numeric::bounds<double>::highest();

		switch(m_method)
		{
		case SIMPLE:
		{
			std::fill(f.begin(),f.end(),high_value);
			break;
		}
		case KURI:
		{
			constraint_vector::size_type number_of_constraints = c.size();
			constraint_vector::size_type number_of_satisfied_constraints = 0;

			// computes the number of satisfied constraints
			for(c_size_type i=0; i<number_of_constraints; i++){
				if(m_original_problem->test_constraint(c,i))
					number_of_satisfied_constraints += 1;
			}

			// sets the Kuri penalization
			double penalization = high_value * (1. - (double)number_of_satisfied_constraints / (double)number_of_constraints);

			std::fill(f.begin(),f.end(),penalization);
			break;
		}
		case WEIGHTED:
		{
			m_original_problem->objfun(f, x);
			const std::vector<double> &c_tol = m_original_problem->get_c_tol();

			// modify equality constraints to behave as inequality constraints:
			c_size_type number_of_constraints = m_original_problem->get_c_dimension();
			c_size_type number_of_eq_constraints = number_of_constraints - m_original_problem->get_ic_dimension();
			double penalization = 0;

			for(c_size_type i=0; i<number_of_constraints; i++) {
					if(i<number_of_eq_constraints){
						c[i] = std::abs(c[i]) - c_tol[i];
					}
					else{
						c[i] = c[i] - c_tol[i];
					}
			}

			for(c_size_type i=0; i<number_of_constraints; i++) {
				if(c[i] > 0.) {
					penalization += m_penalty_factors[i]*c[i];
				}
			}

			// penalizing the objective with the sum of the violation, weighted with the given factors
			for(f_size_type i=0; i<f.size(); i++){
				f[i] += penalization;
			}

			break;
		}
		default:
			pagmo_throw(value_error, "Error: There are only 2 methods for the death penalty!");
			break;
		}
	}
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the type of constraint handling
 */
std::string death_penalty::human_readable_extra() const
{
	std::ostringstream oss;
	oss << m_original_problem->human_readable_extra() << std::endl;
	oss << "\n\tConstraints handled with death penalty, method ";
	switch(m_method){
		case SIMPLE: {
			oss << "SIMPLE ";
			break;
		}
		case KURI: {
			oss << "KURI ";
			break;
			}
		case WEIGHTED: {
			oss << "WEIGHTED ";
			break;
			}
	}
	oss << std::endl;
	return oss.str();
}

std::string death_penalty::get_name() const
{
	std::string method;

	switch(m_method){
		case SIMPLE: {
			method = "SIMPLE ";
			break;
		}
		case KURI: {
			method = "KURI ";
			break;
			}
		case WEIGHTED: {
			method = "WEIGHTED ";
			break;
			}
	}
	return m_original_problem->get_name() + " [death_penalty, method_" + method + "]";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::death_penalty)

