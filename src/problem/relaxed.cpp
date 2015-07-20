/*****************************************************************************
 *   Copyright (C) 2015 The PaGMO development team,                          *
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

#include "../exceptions.h"
#include "../types.h"
#include "../population.h"
#include "base.h"
#include "relaxed.h"

namespace pagmo { namespace problem {

/**
 * Construct the continuous relaxation of a (mixed-)integer problem
 *
 * @param[in] p base::problem to be relaxed
 *
 * @see problem::base constructors.
 */

relaxed::relaxed(const base & p):
		base_meta(
		 p,
		 p.get_dimension(),
		 0,
		 p.get_f_dimension(),
		 p.get_c_dimension(),
		 p.get_ic_dimension(),
		 p.get_c_tol())
{
}

/// Clone method.
base_ptr relaxed::clone() const
{
	return base_ptr(new relaxed(*this));
}

/// Implementation of the objective function.
/// (Wraps over the original implementation)
void relaxed::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	m_original_problem->objfun(f, x);
}

/// Implementation of the constraints computation.
/// (Wraps over the original implementation)
void relaxed::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	m_original_problem->compute_constraints(c, x);
}

std::string relaxed::get_name() const
{
	return m_original_problem->get_name() + " [Relaxed]";
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the translation vector
 */
std::string relaxed::human_readable_extra() const
{
	return m_original_problem->human_readable_extra();
}
}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::relaxed)
