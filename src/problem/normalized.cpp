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
#include "normalized.h"

namespace pagmo { namespace problem {

/**
 * Constructor
 *
 * @param[in] p base::problem to be shifted
 *
 * @see problem::base constructors.
 */

normalized::normalized(const base & p):
	base_meta(
		 p,
		 p.get_dimension(),
		 p.get_i_dimension(),
		 p.get_f_dimension(),
		 p.get_c_dimension(),
		 p.get_ic_dimension(),
		 p.get_c_tol()),
		 m_normalization_center(p.get_dimension(),0),
		 m_normalization_scale(p.get_dimension(),0)
{
	configure_new_bounds();
	std::vector<decision_vector> new_best_x = transform_x(p.get_best_x());
	set_best_x(new_best_x);
}

/// Clone method.
base_ptr normalized::clone() const
{
	return base_ptr(new normalized(*this));
}

/// Set up new bounds for the normalized problem
void normalized::configure_new_bounds()
{
	for(base::size_type i = 0; i < get_lb().size(); ++i){
		double center = (m_original_problem->get_ub()[i] + m_original_problem->get_lb()[i]) / 2;
		double spread = m_original_problem->get_ub()[i] - m_original_problem->get_lb()[i];
		m_normalization_center[i] = center;  // Center to origin
		m_normalization_scale[i] = spread/2; // Scale to [-1, 1] centered at origin
	}
	set_bounds(-1, 1);
}

/// Compute normalized vectors for the meta-problem from the original problem.
/*
 * @param[in] x vectors of the original problem
 * @param[out] vectors x normalized
 */
std::vector<decision_vector> normalized::transform_x(const std::vector<decision_vector> &x) const
{
	const base::size_type cnt = x.size();
	std::vector<decision_vector> new_x = x;

	for (base::size_type i = 0; i < cnt; ++i) {
		for (base::size_type j = 0; j < x[i].size(); ++j) {
			new_x[i][j] = (x[i][j] - m_normalization_center[j]) / m_normalization_scale[j];
		}
	}

	return new_x;
}

/// Returns the de-normalized version of the decision variables
decision_vector normalized::denormalize(const decision_vector& x) const
{
	decision_vector retval(x.size(), 0);
	for(base::size_type i = 0; i < x.size(); ++i){
		retval[i] = (x[i] * m_normalization_scale[i]) + m_normalization_center[i];
	}
	return retval;
}

/// Implementation of the objective function.
/// (Wraps over the original implementation with translated input x)
void normalized::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	decision_vector x_original = denormalize(x);
	m_original_problem->objfun(f, x_original);
}

/// Implementation of the constraints computation.
/// (Wraps over the original implementation with translated input x)
void normalized::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	decision_vector x_original = denormalize(x);
	m_original_problem->compute_constraints(c, x_original);
}


std::string normalized::get_name() const
{
	return m_original_problem->get_name() + " [Normalized]"; 
}

std::string normalized::human_readable_extra() const
{
	std::ostringstream oss;
	oss << m_original_problem->human_readable_extra() << std::endl;
	return oss.str();
}
}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::normalized)
