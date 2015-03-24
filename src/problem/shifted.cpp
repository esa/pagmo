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

#include "../exceptions.h"
#include "../types.h"
#include "../population.h"
#include "base.h"
#include "shifted.h"

namespace pagmo { namespace problem {

/**
 * Construct with a specified shift
 *
 * @param[in] p base::problem to be shifted
 * @param[in] translation The vector specifying the translation
 *
 * @see problem::base constructors.
 */

shifted::shifted(const base & p,
				 const decision_vector & translation):
		base_meta(
		 p,
		 p.get_dimension(),
		 p.get_i_dimension(),
		 p.get_f_dimension(),
		 p.get_c_dimension(),
		 p.get_ic_dimension(),
		 p.get_c_tol()),
		m_translation(translation)
{
	if (translation.size() != p.get_dimension()) {
		pagmo_throw(value_error,"The size of the shifting vector must be equal to the problem dimension");
	}
	configure_shifted_bounds(m_translation);
}


/**
 * Construct with a constant shift (across dimensions)
 *
 * @param[in] p pagmo::problem::base to be shifted
 * @param[in] t Specify the amount of translation, applied universally to all dimensions
 *
 * @see problem::base constructors.
 */

shifted::shifted(const base & p,
				 const double t):
		base_meta(
		 p,
		 p.get_dimension(),
		 p.get_i_dimension(),
		 p.get_f_dimension(),
		 p.get_c_dimension(),
		 p.get_ic_dimension(),
		 p.get_c_tol()),
		m_translation(decision_vector(p.get_dimension(), t))
{
	configure_shifted_bounds(m_translation);
}

/**
 * Construct with a random shift
 *
 * @param[in] p pagmo::problem::base to be shifted
 *
 * @see problem::base constructors.
 */

shifted::shifted(const base & p):
	base_meta(
		 p,
		 p.get_dimension(),
		 p.get_i_dimension(),
		 p.get_f_dimension(),
		 p.get_c_dimension(),
		 p.get_ic_dimension(),
		 p.get_c_tol()),
		m_translation(decision_vector(p.get_dimension(),0))
{
	for (size_t i=0; i< m_translation.size();++i) {
		m_translation[i] = (2*((double) rand() / (RAND_MAX))-1) * (p.get_ub()[i]-p.get_lb()[i]);
	}
	configure_shifted_bounds(m_translation);
}

/// Clone method.
base_ptr shifted::clone() const
{
	return base_ptr(new shifted(*this));
}

/// Set up new bounds for the shifted problem
void shifted::configure_shifted_bounds(const decision_vector & translation)
{
	decision_vector shifted_lb = m_original_problem->get_lb();
	decision_vector shifted_ub = m_original_problem->get_ub();

	for(problem::base::size_type i = 0; i < shifted_lb.size(); ++i) {
		shifted_lb[i] += translation[i];
		shifted_ub[i] += translation[i];
	}
	set_bounds(shifted_lb, shifted_ub);
}

/**
 * Returns the de-shifted version of the decision variables
 *
 * @param[in] x decision vector to be de-shifted
 */
decision_vector shifted::deshift(const decision_vector& x) const
{
	//TODO Need some assertion here?
	decision_vector x_translated(x.size(), 0);
	for(problem::base::size_type i = 0; i < x.size(); ++i) {
		x_translated[i] = x[i] - m_translation[i];
	}
	return x_translated;
}

/// Implementation of the objective function.
/// (Wraps over the original implementation with translated input x)
void shifted::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	decision_vector x_translated = deshift(x);
	m_original_problem->objfun(f, x_translated);
}

/// Implementation of the constraints computation.
/// (Wraps over the original implementation with translated input x)
void shifted::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	decision_vector x_translated = deshift(x);
	m_original_problem->compute_constraints(c, x_translated);
}

/**
 * Gets the shift vector which defines the problem
 *
 * @return the shift std::vector
 */

const decision_vector& shifted::get_shift_vector() const {
	return m_translation;
}

std::string shifted::get_name() const
{
	return m_original_problem->get_name() + " [Shifted]"; 
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the translation vector
 */
std::string shifted::human_readable_extra() const
{
	std::ostringstream oss;
	oss << m_original_problem->human_readable_extra() << std::endl;
	oss << "\n\tShift vector: " << m_translation << std::endl;
	return oss.str();
}
}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::shifted)
