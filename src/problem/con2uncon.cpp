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
#include "base.h"
#include "con2uncon.h"

namespace pagmo { namespace problem {

/**
 * Constructor using initial constrained problem
 *
 * @param[in] problem base::problem to be con2unconed
 */
con2uncon::con2uncon(const base &problem, const method_type &method):
	base((int)problem.get_dimension(),
		 problem.get_i_dimension(),
		 problem.get_f_dimension(),
		 0,
		 0,
		 0.),
	m_original_problem(problem.clone()),
	m_method(method)
{
	set_bounds(m_original_problem->get_lb(),m_original_problem->get_ub());

	if((m_method < 0) || (m_method > 1)) {
		pagmo_throw(value_error,"the problem method must be either OPTIMALITY or FEASIBILITY.");
	}
}

/// Copy Constructor. Performs a deep copy
con2uncon::con2uncon(const con2uncon &prob):
	base((int)prob.get_dimension(),
		 prob.get_i_dimension(),
		 prob.get_f_dimension(),
		 0,
		 0,
		 0.),
	m_original_problem(prob.m_original_problem->clone()),
	m_method(prob.m_method)
{
	set_bounds(m_original_problem->get_lb(),m_original_problem->get_ub());
}

/// Clone method.
base_ptr con2uncon::clone() const
{
	return base_ptr(new con2uncon(*this));
}

/// Implementation of the objective functions.
/// (Wraps over the original implementation)
void con2uncon::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	switch(m_method) {
	case(OPTIMALITY): {
		m_original_problem->objfun(f,x);
		break;
	}
	case(FEASIBILITY): {
		// get the constraints dimension
		problem::base::c_size_type prob_c_dimension = m_original_problem->get_c_dimension();
		problem::base::c_size_type number_of_eq_constraints =
				m_original_problem->get_c_dimension() -
				m_original_problem->get_ic_dimension();

		constraint_vector c(m_original_problem->get_c_dimension(),0);
		m_original_problem->compute_constraints(c,x);

		const std::vector<double> &c_tol = m_original_problem->get_c_tol();

		double c_func = 0.;

		// update equality constraints
		for(problem::base::c_size_type j=0; j<number_of_eq_constraints; j++) {
			if(!m_original_problem->test_constraint(c,j)) {
				c_func += ( std::abs(c.at(j)) - c_tol.at(j) ) * ( std::abs(c.at(j)) - c_tol.at(j) );
			}
		}

		// update inequality constraints
		for(problem::base::c_size_type j=number_of_eq_constraints; j<prob_c_dimension; j++) {
			if(!m_original_problem->test_constraint(c,j)) {
				c_func += std::abs(c.at(j)) - c_tol.at(j);
			}
		}

		std::fill(f.begin(),f.end(), 0.);
		f[0] = c_func;
		break;
	}
	}
}

/// Implementation of fitness vectors comparison.
/**
 * @brief compare_fitness_impl calls the compare_fitness method of the original problem.
 * @return true if v_f1 is dominating v_f2, false otherwise.
 */
bool con2uncon::compare_fitness_impl(const fitness_vector &v_f1, const fitness_vector &v_f2) const
{
	return m_original_problem->compare_fitness(v_f1,v_f2);
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the type of constraint handling
 */
std::string con2uncon::human_readable_extra() const
{
	std::ostringstream oss;
	oss << m_original_problem->human_readable_extra() << std::endl;
	oss << "\n\tcon2unconed with method ";
	switch(m_method){
	case OPTIMALITY: {
		oss << "OPTIMALITY ";
		break;
	}
	case FEASIBILITY: {
		oss << "FEASIBILITY ";
		break;
	}
	}
	oss << std::endl;
	return oss.str();
}

std::string con2uncon::get_name() const
{
	std::string method;

	switch(m_method){
	case OPTIMALITY: {
		method = "OPTIMALITY ";
		break;
	}
	case FEASIBILITY: {
		method = "FEASIBILITY ";
		break;
	}
	}
	return m_original_problem->get_name() + " [con2uncon, method_" + method + "]";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::con2uncon);

