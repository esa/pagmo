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
#include <algorithm>

#include "../exceptions.h"
#include "../types.h"
#include "../population.h"
#include "base.h"
#include "death_penalty.h"

namespace pagmo { namespace problem {

/**
 * Constructor using initial constrained problem
 *
 * @param[in] problem base::problem to be modified to use a death penalty
 * as constraints handling technique.
 * @param[in] death_penalty_methos int to be modified to use a simple death penalty
 * if defined with SIMPLE and a Kuri death penalty with KURI.
 *
 */
death_penalty::death_penalty(const base &problem, const method_type method):
	base((int)problem.get_dimension(),
		 problem.get_i_dimension(),
		 problem.get_f_dimension(),
		 0,
		 0,
		 0.),
	m_original_problem(problem.clone()),
	m_method(method)
{
	if(m_original_problem->get_c_dimension() <= 0){
		pagmo_throw(value_error,"The original problem has no constraints.");
	}

	if( method > 1 || method < 0) {
		pagmo_throw(value_error, "the death penalty method must be set to 0 for simple death or to 1 for Kuri death.");
	}

	set_bounds(m_original_problem->get_lb(),m_original_problem->get_ub());

}

/// Copy Constructor. Performs a deep copy
death_penalty::death_penalty(const death_penalty &prob):
	base((int)prob.get_dimension(),
		 prob.get_i_dimension(),
		 prob.get_f_dimension(),
		 prob.get_c_dimension(),
		 prob.get_ic_dimension(),
		 prob.get_c_tol()),
	m_original_problem(prob.m_original_problem->clone()),
	m_method(prob.m_method)
{
	set_bounds(m_original_problem->get_lb(),m_original_problem->get_ub());
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
		default:
			pagmo_throw(value_error, "Error: There are only 2 methods for the death penalty!");
			break;
		}
	}
}

/// Implementation of fitness vectors comparison.
/**
 * @brief compare_fitness_impl calls the compare_fitness method of the original problem.
 * @return true if v_f1 is dominating v_f2, false otherwise.
 */
bool death_penalty::compare_fitness_impl(const fitness_vector &v_f1, const fitness_vector &v_f2) const
{
	return m_original_problem->compare_fitness(v_f1,v_f2);
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
	}
	return m_original_problem->get_name() + " [death_penalty, method_" + method + "]";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::death_penalty);

