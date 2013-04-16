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
#include "rotated.h"

namespace pagmo { namespace problem {

/**
 * Default constructor so that boost::serialization does not complain
 *
 */

rotated::rotated():
	base(30,0,2)
{
	//Or is it better to override the load_construct_data...?
}

/**
 * Will construct the rotated meta-problem.
 *
 * @param[rotation]: Rotate the problem using this rotation matrix
 *
 * @see problem::base constructors.
 */

rotated::rotated(const base_ptr & problem,
				 const Eigen::MatrixXd & rotation):
	base((int)problem->get_dimension(), // But ambiguous without the cast?
		 problem->get_i_dimension(),
		 problem->get_f_dimension(),
		 problem->get_c_dimension(),
		 problem->get_ic_dimension(),
		 problem->get_c_tol()),
	m_original_problem(problem->clone()), //TODO: to clone or not to clone?
	m_Rotate(rotation)
{
	//m_InvRotate = m_Rotate.fullPivLu().inverse();
	m_InvRotate = m_Rotate.transpose();
	configure_shifted_bounds(rotation, problem->get_lb(), problem->get_ub());
}

/// Clone method.
base_ptr rotated::clone() const
{
	return base_ptr(new rotated(*this));
}

/// Update the new bounds due to the rotation transformation
// Slight twist here: Rotation causes the new bounds to be
// not seperable w.r.t the axes. Here, the approach taken
// is to construct a new constraint hyperbox that is a relaxation of
// the real non-separable constraints. All the points in the original
// search space will be searchable in the new constraint box. Some
// additional points (that may be invalid) will be introduced this way,
// so in the objfun_impl all these out-of-bounds point will be projected
// back to the valid bounds.
// (Note: Rot is used here just for othorgonality checking)
void rotated::configure_shifted_bounds(const Eigen::MatrixXd & Rot,
							  const decision_vector & original_lb,
							  const decision_vector & original_ub)
{	
	Eigen::MatrixXd check = m_InvRotate * Rot;
	if(!check.isIdentity()){
		// TODO: Is there a better way than this?
		std::cout<<"Warning: Rotation matrix is not orthonormal!"<<std::endl;
	}
	//Normalize to [-1, 1]
	m_normalize_translation.clear();
	m_normalize_scale.clear();
	for(base::size_type i = 0; i < original_lb.size(); i++){
		double mean_t = (original_ub[i] + original_lb[i]) / 2;
		double spread_t = original_ub[i] - original_lb[i]; // what if zero?!
		m_normalize_translation.push_back(mean_t); // Center to origin
		m_normalize_scale.push_back(spread_t/2); // Scale to [-1, 1] centered at origin
	}
	decision_vector l_updated_lb = normalize_to_center(original_lb);
	decision_vector l_updated_ub = normalize_to_center(original_ub);
	// Expand the box to cover the whole original search space
	for(base::size_type i = 0; i < original_lb.size(); i++){
		l_updated_lb[i] = sqrt(2.0) * l_updated_lb[i];
		l_updated_ub[i] = sqrt(2.0) * l_updated_ub[i];
	}
	set_bounds(l_updated_lb, l_updated_ub);
}

// Used to normalize the original upper and lower bounds
// to [-1, 1], at each dimension
decision_vector rotated::normalize_to_center(const decision_vector& x) const
{
	decision_vector normalized_x(x.size(), 0);
	for(base::size_type i = 0; i < x.size(); i++){
		normalized_x[i] = (x[i] - m_normalize_translation[i]) / m_normalize_scale[i];
	}
	return normalized_x;
}

// For a decision_vector in the normalized [-1, 1] space,
// perform the inverse operation of normalization to get it's
// location in the original space
decision_vector rotated::denormalize_to_original(const decision_vector& x_normed) const
{
	decision_vector denormalized_x(x_normed.size(), 0);
	for(base::size_type i = 0; i < x_normed.size(); i++){
		denormalized_x[i] = (x_normed[i] * m_normalize_scale[i]) + m_normalize_translation[i];
	}
	return denormalized_x;
}

// Clip the variables to the valid bounds
decision_vector rotated::projection_via_clipping(const decision_vector& x) const
{
	decision_vector x_clipped = x;
	decision_vector original_lb = m_original_problem->get_lb();
	decision_vector original_ub = m_original_problem->get_ub();
	for(base::size_type i = 0; i < x_clipped.size(); i++){
		x_clipped[i] = std::max(x_clipped[i], original_lb[i]);
		x_clipped[i] = std::min(x_clipped[i], original_ub[i]);
	}
	return x_clipped;
}

/// Returns the de-rotated version of the decision variables,
/// obtained from inversed rotation operation, which is ready to be fed
/// to the original obj. function
decision_vector rotated::get_inv_rotated_vars(const decision_vector& x_normed) const
{
	// This may be outside of the original domain, due to the 
	// relaxed variable bounds after rotation -- project it back if so.

	// 1. De-rotate the vector in the normalized space
	Eigen::VectorXd x_normed_vec = Eigen::VectorXd::Zero(x_normed.size());
	Eigen::VectorXd x_derotated_vec;
	for(base::size_type i = 0; i < x_normed.size(); i++){
		x_normed_vec(i) = x_normed[i];	
	}
	x_derotated_vec = m_InvRotate * x_normed_vec;

	// 2. De-normalize the de-rotated vector to the original bounds
	decision_vector x_derotated(x_normed.size(), 0);
	for(base::size_type i = 0; i < x_normed.size(); i++){
		x_derotated[i] = x_derotated_vec(i);
	}
	decision_vector x_wild = denormalize_to_original(x_derotated);

	// 3. The de-normalized vector may be out of bounds, if so project back in
	decision_vector x = projection_via_clipping(x_wild);

	return x;
}

/// Implementation of the objective function.
/// (Wraps over the original implementation with de-rotated input)
void rotated::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	decision_vector x_inv_rotated = get_inv_rotated_vars(x);
	m_original_problem->objfun(f, x_inv_rotated);
}

/// Implementation of the constraints computation.
/// (Wraps over the original implementation with de-rotated input)
void rotated::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	decision_vector x_inv_rotated = get_inv_rotated_vars(x);
	m_original_problem->compute_constraints(c, x_inv_rotated);
}

std::string rotated::get_name() const
{
	return m_original_problem->get_name() + " [Rotated]"; 
}
}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::rotated);

