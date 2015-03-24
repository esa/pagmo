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
#include "base.h"
#include "rotated.h"

namespace pagmo { namespace problem {

/**
 * Constructor using Eigen Matrix
 *
 * @param[in] p base::problem to be rotated
 * @param[in] rotation Eigen::MatrixXd expressing the problem rotation
 *
 * @see problem::base constructors.
 */

rotated::rotated(const base &p, const Eigen::MatrixXd &rotation ):
		base_meta(
		 p,
		 p.get_dimension(),
		 p.get_i_dimension(),
		 p.get_f_dimension(),
		 p.get_c_dimension(),
		 p.get_ic_dimension(),
		 p.get_c_tol()),
	m_Rotate(rotation), m_normalize_translation(), m_normalize_scale()
{
	m_InvRotate = m_Rotate.transpose();
	
	Eigen::MatrixXd check = m_InvRotate * m_Rotate;
	if(!check.isIdentity(1e-5)){
		pagmo_throw(value_error,"The input matrix seems not to be orthonormal (to a tolerance of 1e-5)");
	}
	if(p.get_i_dimension()>0){
		pagmo_throw(value_error,"Input problem has an integer dimension. Cannot rotate it.");
	}
	configure_new_bounds();
}

/**
 * Constructor using std::vector (for python exposition purposes)
 *
 * @param[in] p base::problem to be rotated
 * @param[in] rotation std::vector<std::vector<double> > expressing the problem rotation
 *
 * @see problem::base constructors.
 */
rotated::rotated(const base &p,
				 const std::vector<std::vector<double> > &rotation):
		base_meta(
		 p,
		 p.get_dimension(),
		 p.get_i_dimension(),
		 p.get_f_dimension(),
		 p.get_c_dimension(),
		 p.get_ic_dimension(),
		 p.get_c_tol()),
	m_Rotate(),m_normalize_translation(), m_normalize_scale()
{
	if(!(rotation.size()==get_dimension())){
			pagmo_throw(value_error,"The input matrix dimensions seem incorrect");
	}
	if(p.get_i_dimension()>0){
		pagmo_throw(value_error,"Input problem has an integer dimension. Cannot rotate it.");
	}
	m_Rotate.resize(rotation.size(),rotation.size());
	for (base::size_type i = 0; i < rotation.size(); ++i) {
		if(!(rotation.size()==rotation[i].size())){
			pagmo_throw(value_error,"The input matrix seems not to be square");
		}
		for (base::size_type j = 0; j < rotation[i].size(); ++j) {
			m_Rotate(i,j) = rotation[i][j];
		}
	}
	m_InvRotate = m_Rotate.transpose();
	
	Eigen::MatrixXd check = m_InvRotate * m_Rotate;
	if(!check.isIdentity(1e-5)){
		pagmo_throw(value_error,"The input matrix seems not to be orthonormal (to a tolerance of 1e-5)");
	}
	configure_new_bounds();
}

/**
 * Constructor with a random matrix
 *
 * @param[in] p base::problem to be rotated
 */

rotated::rotated(const base &p):
		base_meta(
		 p,
		 p.get_dimension(),
		 p.get_i_dimension(),
		 p.get_f_dimension(),
		 p.get_c_dimension(),
		 p.get_ic_dimension(),
		 p.get_c_tol()),
	m_normalize_translation(), m_normalize_scale()
{
	size_type dim = p.get_dimension();
	m_Rotate = Eigen::MatrixXd::Random(dim, dim).householderQr().householderQ();
	m_InvRotate = m_Rotate.transpose();

	Eigen::MatrixXd check = m_InvRotate * m_Rotate;
	if(!check.isIdentity(1e-5)){
		pagmo_throw(value_error,"The input matrix seems not to be orthonormal (to a tolerance of 1e-5)");
	}
	if(p.get_i_dimension()>0){
		pagmo_throw(value_error,"Input problem has an integer dimension. Cannot rotate it.");
	}
	configure_new_bounds();
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
// so in the objfun_impl all these out-of-bounds point will have to be projected
// back to valid bounds.

void rotated::configure_new_bounds()
{	

	for(base::size_type i = 0; i < get_lb().size(); i++){
		double mean_t = (m_original_problem->get_ub()[i] + m_original_problem->get_lb()[i]) / 2;
		double spread_t = m_original_problem->get_ub()[i] - m_original_problem->get_lb()[i]; // careful, if zero needs to be acounted for later
		m_normalize_translation.push_back(mean_t); // Center to origin
		m_normalize_scale.push_back(spread_t/2); // Scale to [-1, 1] centered at origin
	}
	// Expand the box to cover the whole original search space. We may here call directly
	// the set_bounds(const double &, const double &) as all dimensions are now equal
	set_bounds(-sqrt(2), sqrt(2));
}

// Used to normalize the original upper and lower bounds
// to [-1, 1], at each dimension
decision_vector rotated::normalize_to_center(const decision_vector& x) const
{
	decision_vector normalized_x(x.size(), 0);
	for(base::size_type i = 0; i < x.size(); i++) {
		if (m_normalize_scale[i] == 0) { //If the bounds witdth is zero
			normalized_x[i] = 0;
		} 
		else {
			normalized_x[i] = (x[i] - m_normalize_translation[i]) / m_normalize_scale[i];
		}
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
	for(base::size_type i = 0; i < x_clipped.size(); i++){
		x_clipped[i] = std::max(x_clipped[i], m_original_problem->get_lb()[i]);
		x_clipped[i] = std::min(x_clipped[i], m_original_problem->get_ub()[i]);
	}
	return x_clipped;
}

/// Returns the original version of the decision variables ready to be fed
/// to the original problem
decision_vector rotated::derotate(const decision_vector& x_normed) const
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
	m_original_problem->objfun(f, derotate(x));
}

/// Implementation of the constraints computation.
/// (Wraps over the original implementation with de-rotated input)
void rotated::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	m_original_problem->compute_constraints(c, derotate(x));
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the string representation of the rotation matrix 
 */
std::string rotated::human_readable_extra() const
{
	std::ostringstream oss;
	oss << m_original_problem->human_readable_extra() << std::endl;
	oss << "\n\tRotation matrix: " << std::endl;
	if (m_Rotate.cols() > 5) {
		oss << m_Rotate.block(0,0,5,5) << std::endl;
		oss << "..." << std::endl;
	}
	else {
		oss << m_Rotate << std::endl;
	}
	return oss.str();
}

std::string rotated::get_name() const
{
	return m_original_problem->get_name() + " [Rotated]"; 
}

/**
 * Gets the rotation matrix
 *
 * @return an orthonormal Eigen::MatrixXd defining part of the transformation applied to the original problem
 */
const Eigen::MatrixXd& rotated::get_rotation_matrix() const {
	return m_Rotate;
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::rotated)

