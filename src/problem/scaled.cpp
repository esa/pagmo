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
#include "scaled.h"

namespace pagmo { namespace problem {

/**
 * Construct with a vector of new units
 *
 * @param[in] p base::problem to be shifted
 * @param[in] units The vector specifying the new fitness units
 *
 * @see problem::base constructors.
 */

scaled::scaled(const base & p, const decision_vector & units):
		base_meta(
		 p,
		 p.get_dimension(),
		 p.get_i_dimension(),
		 p.get_f_dimension(),
		 p.get_c_dimension(),
		 p.get_ic_dimension(),
		 p.get_c_tol()),
		m_units(units)
{
	if (units.size() != p.get_f_dimension()) {
		pagmo_throw(value_error,"The size of the provided units vector does not match the fitness dimension of the provided problem");
	}
}

/// Clone method.
base_ptr scaled::clone() const
{
	return base_ptr(new scaled(*this));
}

/**
 * Returns the de-scaled fitness 
 *
 * @param[in] f fitness vector to be de-scaled
 */
fitness_vector scaled::descale(const fitness_vector& f) const
{
	fitness_vector retval(f);
	for (fitness_vector::size_type i=0; i< f.size();++i) {
		retval[i] *= m_units[i];
	}
	return retval;
}

/**
 * Gets the units vector defining the problem
 *
 * @return the units fitness_vector
 */

const fitness_vector& scaled::get_units() const {
	return m_units;
}

/// Implementation of the objective function.
void scaled::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	m_original_problem->objfun(f, x);
	for (fitness_vector::size_type i=0; i< f.size();++i) {
		f[i] /= m_units[i];
	}
	
}

/// Implementation of the constraints computation.
void scaled::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	m_original_problem->compute_constraints(c, x);
}


std::string scaled::get_name() const
{
	return m_original_problem->get_name() + " [Scaled]"; 
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the translation vector
 */
std::string scaled::human_readable_extra() const
{
	std::ostringstream oss;
	oss << m_original_problem->human_readable_extra() << std::endl;
	oss << "\n\tScaled units vector: " << m_units << std::endl;
	return oss.str();
}
}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::scaled)
