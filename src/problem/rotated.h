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

#ifndef PAGMO_PROBLEM_ROTATED_H
#define PAGMO_PROBLEM_ROTATED_H

#include <string>

#include "../serialization.h"
#include "../types.h"
#include "ackley.h"
#include "base_meta.h"
#include "../Eigen/Dense"

namespace pagmo{ namespace problem {

/// Shifted meta-problem
/**
 * Implements a meta-problem class that wraps some other problems,
 * resulting in a rotated version of the underlying problem.
 *
 * @author Yung-Siang Liau (liauys@gmail.com)
 */

class __PAGMO_VISIBLE rotated : public base_meta
{
	public:
		//constructors
		rotated(const base &, const Eigen::MatrixXd &);
		rotated(const base &, const std::vector<std::vector<double> > &);
		rotated(const base & = ackley(1));
		
		base_ptr clone() const;
		std::string get_name() const;
		
		decision_vector derotate(const decision_vector &) const;
		const Eigen::MatrixXd& get_rotation_matrix() const;

	protected:
		std::string human_readable_extra() const;
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		void compute_constraints_impl(constraint_vector &, const decision_vector &) const;

	private:
		void configure_new_bounds();

		decision_vector normalize_to_center(const decision_vector& x) const;
		decision_vector denormalize_to_original(const decision_vector& x) const;
		decision_vector projection_via_clipping(const decision_vector& x) const;
	
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base_meta>(*this);
			ar & m_Rotate;
			ar & m_InvRotate;
			ar & m_normalize_translation;
			ar & m_normalize_scale;
		}
		Eigen::MatrixXd m_Rotate;
		Eigen::MatrixXd m_InvRotate;
		decision_vector m_normalize_translation;
		decision_vector m_normalize_scale;


};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::rotated)

#endif // PAGMO_PROBLEM_ROTATED_H
