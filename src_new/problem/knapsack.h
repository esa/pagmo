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

#ifndef PAGMO_PROBLEM_KNAPSACK_H
#define PAGMO_PROBLEM_KNAPSACK_H

#include <boost/numeric/conversion/cast.hpp>
#include <cstddef>
#include <vector>

#include "../config.h"
#include "../types.h"
#include "base.h"

namespace pagmo { namespace problem {

/// N-dimensional knapsack problem.
/**
 *
 * \image html knapsack.png "Knapsack problem."
 * \image latex knapsack.png "Knapsack problem." width=3cm
 *
 * Classical formulation of the 0-1 knapsack problem: given N items, each one with a weight and a monetary value, determine
 * which items to include in the knapsack so that the total weight is less than a given limit and the total value is as large as possible.
 *
 * Mathematically, the problem is formulated as follows:
 * \f[
 * 	\begin{array}{ll}
 * 	\textnormal{maximise:} & \sum_{i=1}^Np_ix_i, \\
 * 	\textnormal{subject to:} & \sum_{i=1}^Nw_ix_i \leq W, x_i \in \left\{ 0,1 \right\},
 * 	\end{array}
 * \f]
 * where \f$ p_i \f$ is the value of the item and \f$ w_i \f$ its weight.
 *
 * In PaGMO's terminology, this problem has global and integer dimensions equal to N, fitness dimension equal to 1, global and inequality constraints
 * dimensions equal to 1.
 *
 * @see http://en.wikipedia.org/wiki/Knapsack_problem
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE knapsack: public base
{
	public:
		knapsack(const std::vector<double> &, const std::vector<double> &, const double &);
		/// Constructor from raw arrays and maximum weight.
		/**
		 * Initialise the values and weights of the items from raw arrays, and maximum weight to max_weight. Will fail if max_weight is negative,
		 * if N is 0 or any weight/value is negative.
		 */
		template <std::size_t N>
		knapsack(const double (&values)[N], const double (&weights)[N], const double &max_weight):base(boost::numeric_cast<int>(N),boost::numeric_cast<int>(N),1,1,1),
			m_values(values,values + N),m_weights(weights,weights + N),m_max_weight(max_weight)
		{
			verify_init();
		}
		base_ptr clone() const;
	protected:
		void compute_constraints_impl(constraint_vector &, const decision_vector &) const;
		bool compare_fitness_impl(const fitness_vector &, const fitness_vector &) const;
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		bool equality_operator_extra(const base &) const;
	private:
		void verify_init() const;
	private:
		const std::vector<double>	m_values;
		const std::vector<double>	m_weights;
		const double			m_max_weight;
};

}
}

#endif
