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

#include <algorithm>
#include <boost/integer_traits.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <climits>
#include <cstddef>
#include <exception>
#include <iterator>
#include <typeinfo>

#include <string>

#include "../exceptions.h"
#include "base.h"
#include "golomb_ruler.h"

namespace pagmo { namespace problem {

static inline int check_golomb_order(int n)
{
	if (n < 2) {
		pagmo_throw(value_error,"Golomb ruler problem must have at least order 2");
	}
	return n;
}

/// Constructor from order and maximum distance between consecutive marks.
/**
 * Construct a Golomb ruler problem of order n and whose maximum distance between consecutive marks is m.
 *
 * @param[in] n order of the Golomb ruler.
 * @param[in] m upper limit for the distance between consecutive marks.
 */
golomb_ruler::golomb_ruler(int n, int m):base(check_golomb_order(n) - 1,n - 1,1,1,0),m_max_length(boost::numeric_cast<std::size_t>(m))
{
	if (!m_max_length || m_max_length > static_cast<std::size_t>(INT_MAX)) {
		pagmo_throw(value_error,"maximum distance between consecutive marks must be in the ]0,32767] range");
	}
	if (get_dimension() == boost::integer_traits<size_type>::const_max) {
		pagmo_throw(std::overflow_error,"size overflow in Golomb ruler problem");
	}
	if (double(m_max_length) * get_dimension() > static_cast<double>(INT_MAX)) {
		pagmo_throw(std::overflow_error,"fitness value overflow in Golomb ruler problem");
	}
	// Lower bound is already set to 0.
	set_ub(m_max_length);
}

/// Clone method.
base_ptr golomb_ruler::clone() const
{
	return base_ptr(new golomb_ruler(*this));
}

/// Additional requirements for equality.
/**
 * @return true if upper limits for length are equal, false otherwise.
 */
bool golomb_ruler::equality_operator_extra(const base &other) const
{
	pagmo_assert(typeid(*this) == typeid(other));
	return (m_max_length == dynamic_cast<golomb_ruler const &>(other).m_max_length);
}

/// Implementation of the objective function.
/**
 * Will return the distance of the ruler.
 *
 * @param[out] f pagmo::fitness_vector that stores the output fitness.
 * @param[in] x pagmo::decision_vector whose fitness will be calculated.
 */
void golomb_ruler::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == 1 && x.size() == get_dimension());
	compute_marks_and_dist(x);
	// Fitness is the maximum distance.
	f[0] = *std::max_element(m_tmp_dist.begin(),m_tmp_dist.end());
}

/// Implementation of constraint calculation.
/**
 * Implements an equality constraint on the number of equal distances between the marks pairs. If this number is 0,
 * the constraint is satisfied and the ruler is a Golomb ruler.
 *
 * @param[out] c pagmo::constraint_vector that stores the output constraint.
 * @param[in] x pagmo::decision_vector whose feasibility will be calculated.
 */
void golomb_ruler::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	pagmo_assert(c.size() == 1 && x.size() == get_dimension());
	compute_marks_and_dist(x);
	// Sort the vector of distances.
	std::sort(m_tmp_dist.begin(),m_tmp_dist.end());
	// Now compute how many duplicate distances are there.
	c[0] = boost::numeric_cast<double>(m_tmp_dist.size()) - std::distance(m_tmp_dist.begin(),std::unique(m_tmp_dist.begin(),m_tmp_dist.end()));
}

// Compute marks and distances of x and store them internally.
void golomb_ruler::compute_marks_and_dist(const decision_vector &x) const
{
	// We already computed distances and marks of this decision vector, do not do anything.
	if (x == m_tmp_x) {
		return;
	}
	m_tmp_x = x;
	const size_type size = m_tmp_x.size(), marks_size =  size + 1;
	m_tmp_marks.resize(marks_size);
	m_tmp_marks[0] = 0;
	// Write marks into temporary vector.
	for (size_type i = 0; i < size; ++i) {
		m_tmp_marks[i + 1] = m_tmp_marks[i] + m_tmp_x[i];
	}
	// Clear distances vector.
	m_tmp_dist.clear();
	for (size_type i = 0; i < marks_size - 1; ++i) {
		for (size_type j = i + 1; j < marks_size; ++j) {
			const double tmp = m_tmp_marks[j] - m_tmp_marks[i];
			pagmo_assert(tmp >= 0);
			m_tmp_dist.push_back(tmp);
		}
	}
}

std::string golomb_ruler::get_name() const
{
	return "Golomb ruler";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::golomb_ruler)
