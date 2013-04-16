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
#include "shifted.h"

namespace pagmo { namespace problem {

/**
 * Default constructor so that boost::serialization does not complain
 *
 */

shifted::shifted():
	base(30,0,2)
{
	//Or is it better to override the load_construct_data...?
}

/**
 * Will construct the shifted meta-problem.
 *
 * @param[translation]: Specify the amount of translation in each dimension
 *
 * @see problem::base constructors.
 */

shifted::shifted(const base_ptr & problem,
				 const decision_vector & translation):
	base((int)problem->get_dimension(), // But ambiguous without the cast?
		 problem->get_i_dimension(),
		 problem->get_f_dimension(),
		 problem->get_c_dimension(),
		 problem->get_ic_dimension(),
		 problem->get_c_tol()),
	m_original_problem(problem->clone()), //TODO: to clone or not to clone?
	m_translation(translation)
{
	//TODO Ok to do this here?
	configure_shifted_bounds(m_translation, problem->get_lb(), problem->get_ub());
}


/**
 * Will construct the shifted meta-problem.
 *
 * @param[t]: Specify the amount of translation, applied universally to all dimensions
 *
 * @see problem::base constructors.
 */

shifted::shifted(const base_ptr & problem,
				 const double t):
	base(problem->get_dimension(),0,problem->get_f_dimension()),
	m_original_problem(problem->clone()),
	m_translation(decision_vector(problem->get_dimension(), t))
{
	configure_shifted_bounds(m_translation, problem->get_lb(), problem->get_ub());
}

/// Clone method.
base_ptr shifted::clone() const
{
	return base_ptr(new shifted(*this));
}

/// Set up new bounds for the shifted problem
void shifted::configure_shifted_bounds(const decision_vector & translation,
							  const decision_vector & original_lb,
							  const decision_vector & original_ub)
{
	decision_vector shifted_lb = original_lb;
	decision_vector shifted_ub = original_ub;

	for(problem::base::size_type i = 0; i < shifted_lb.size(); ++i){
		shifted_lb[i] = original_lb[i] + translation[i];
		shifted_ub[i] = shifted_ub[i] + translation[i];
	}
	//Order of update matters here!
	//(As bound sanity check is invoked immediately after each update,
	// even if only lb or ub is changed...)
	for(problem::base::size_type i = 0; i < shifted_lb.size(); ++i){
		if(shifted_lb[i] > original_ub[i]){
			set_ub(i, shifted_ub[i]);
			set_lb(i, shifted_lb[i]);
		}
		else{
			set_lb(i, shifted_lb[i]);
			set_ub(i, shifted_ub[i]);
		}
	}
}

/// Returns the shifted version of the decision variables,
/// obtained from inversed translation
decision_vector shifted::get_shifted_vars(const decision_vector& x) const
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
	decision_vector x_translated = get_shifted_vars(x);
	m_original_problem->objfun(f, x_translated);
}

/// Implementation of the constraints computation.
/// (Wraps over the original implementation with translated input x)
void shifted::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	decision_vector x_translated = get_shifted_vars(x);
	m_original_problem->compute_constraints(c, x_translated);
}

std::string shifted::get_name() const
{
	return m_original_problem->get_name() + " [Shifted]"; 
}
}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::shifted);
