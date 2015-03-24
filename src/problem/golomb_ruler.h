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

#ifndef PAGMO_PROBLEM_GOLOMB_RULER_H
#define PAGMO_PROBLEM_GOLOMB_RULER_H

#include <cstddef>
#include <string>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "base.h"

namespace pagmo { namespace problem {

/// Golomb ruler problem.
/**
 * \image html golomb_ruler.png "Optimal and perfect Golomb ruler of order 4."
 * \image latex golomb_ruler.png "Optimal and perfect Golomb ruler of order 4." width=2cm
 *
 * A Golomb ruler is a set of marks at integer positions along an imaginary ruler such that no two pairs of marks
 * are the same distance apart. The number of marks on the ruler is its order, and the largest distance between two of its marks is its length.
 * A Golomb ruler is optimal if no shorter Golomb ruler of the same order exists.
 *
 * This problem is setup to look for optimal Golomb rulers of a given order. The problem has dimension (order - 1), with the decision vector
 * representing distances between successive marks. The objective is to minimise the length of the ruler, with the equality constraint that
 * for each possible pair of marks the distance must be unique (i.e., the number of duplicate distances must be null in order to satisfy the constraint).
 * Note that when this constraint is not satisfied, the ruler is not a Golomb ruler.
 *
 * @see http://en.wikipedia.org/wiki/Golomb_ruler
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE golomb_ruler: public base
{
	public:
		golomb_ruler(int = 5,int = 10);
		base_ptr clone() const;
		std::string get_name() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		void compute_constraints_impl(constraint_vector &, const decision_vector &) const;
		bool equality_operator_extra(const base &) const;
	private:
		void compute_marks_and_dist(const decision_vector &) const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & const_cast<std::size_t &>(m_max_length);
			ar & m_tmp_x;
			ar & m_tmp_marks;
			ar & m_tmp_dist;
		}
		const std::size_t	m_max_length;
		mutable decision_vector	m_tmp_x;
		mutable decision_vector	m_tmp_marks;
		mutable decision_vector	m_tmp_dist;
};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::problem::golomb_ruler)

#endif
