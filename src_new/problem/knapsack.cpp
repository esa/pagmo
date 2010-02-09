/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
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

#include <boost/numeric/conversion/cast.hpp>
#include <vector>

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "knapsack.h"

namespace pagmo { namespace problem {

/// Constructor from vectors and maximum weight.
/**
 * Initialise the values and weights of the items from vectors, and maximum weight to max_weight. Will fail if max_weight is negative,
 * if vector sizes are not equal or null, or if any weight/value is negative.
 */
knapsack::knapsack(const std::vector<double> &values, const std::vector<double> &weights, const double &max_weight):
	base(boost::numeric_cast<int>(values.size()),boost::numeric_cast<int>(values.size()),1,1,1),
	m_values(values),m_weights(weights),m_max_weight(max_weight)
{
	verify_init();
}

/// Clone method.
base_ptr knapsack::clone() const
{
	return base_ptr(new knapsack(*this));
}

/// Implementation of the objective function.
void knapsack::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == 1);
	pagmo_assert(x.size() == get_dimension() && x.size() == m_values.size());
	f[0] = 0;
	for (size_type i = 0; i < get_dimension(); ++i) {
		f[0] += m_values[i] * x[i];
	}
}

/// Re-implement default fitness comparison from problem::base.
/**
 * We need to maximise the weight, while the default fitness comparison minimises.
 */
bool knapsack::compare_f_impl(const fitness_vector &v_f1, const fitness_vector &v_f2) const
{
	pagmo_assert(v_f1.size() == 1 && v_f2.size() == v_f1.size());
	return v_f1[0] > v_f2[0];
}

/// Re-implement constraint computation,
void knapsack::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	pagmo_assert(c.size() == 1 && x.size() == m_values.size());
	c[0] = 0;
	for (size_type i = 0; i < get_dimension(); ++i) {
		c[0] += m_weights[i] * x[i];
	}
	c[0] -= m_max_weight;
}

/// Additional requirements for equality.
/**
 * Return true if items values and weights are the same.
 */
bool knapsack::equality_operator_extra(const base &other) const
{
	pagmo_assert(typeid(*this) == typeid(other));
	return (m_values == dynamic_cast<knapsack const &>(other).m_values && m_weights == dynamic_cast<knapsack const &>(other).m_weights);
}

// Verify that sane values have been input during construction.
void knapsack::verify_init() const
{
	if (m_values.size() != m_weights.size() || m_max_weight < 0) {
		pagmo_throw(value_error,"invalid value(s) in construction of the knapsack problem");
	}
	for (std::vector<double>::size_type i = 0; i < m_values.size(); ++i) {
		if (m_values[i] < 0 || m_weights[i] < 0) {
			pagmo_throw(value_error,"invalid value(s) in construction of the knapsack problem");
		}
	}
}

}
}