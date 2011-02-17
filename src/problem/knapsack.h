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
#include <string>
#include <vector>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "base_aco.h"

namespace pagmo { namespace problem {

/// N-dimensional knapsack problem.
/**
 *
 * \image html knapsack.png "Knapsack problem."
 * \image latex knapsack.png "Knapsack problem." width=3cm
 *
 * This is a constrained integer single-objective problem.
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
 * NOTE: this problem calls the virtual function base_aco::set_heuristic_information_matrix() in its constructor, hence it is advisable *not* to
 * use this class as a base for another class. See http://www.artima.com/cppsource/nevercall.html for a discusssion.
 *
 * @see http://en.wikipedia.org/wiki/Knapsack_problem
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE knapsack: public base_aco
{
	public:
		knapsack();
		knapsack(const std::vector<double> &, const std::vector<double> &, const double &);
		/// Constructor from raw arrays and maximum weight.
		/**
		 * Initialise the values and weights of the items from raw arrays, and maximum weight to max_weight. Will fail if max_weight is negative,
		 * if N is 0 or any weight/value is negative.
		 *
		 * @param[in] values raw array of values.
		 * @param[in] weights raw array of weights.
		 * @param[in] max_weight maximum weight.
		 */
		template <std::size_t N>
		knapsack(const double (&values)[N], const double (&weights)[N], const double &max_weight):base_aco(boost::numeric_cast<int>(N),1,1),
			m_values(values,values + N),m_weights(weights,weights + N),m_max_weight(max_weight)
		{
			verify_init();
			set_heuristic_information_matrix();
		}
		base_ptr clone() const;
		std::string get_name() const;
		bool check_partial_feasibility(const decision_vector &x) const;
	protected:
		void compute_constraints_impl(constraint_vector &, const decision_vector &) const;
		bool compare_fitness_impl(const fitness_vector &, const fitness_vector &) const;
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		bool equality_operator_extra(const base &) const;
		std::string human_readable_extra() const;
	private:
		void verify_init() const;
		void set_heuristic_information_matrix();
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & const_cast<std::vector<double> &>(m_values);
			ar & const_cast<std::vector<double> &>(m_weights);
			ar & const_cast<double &>(m_max_weight);
		}
		const std::vector<double>	m_values;
		const std::vector<double>	m_weights;
		const double			m_max_weight;
};

}
}

BOOST_CLASS_EXPORT_KEY(pagmo::problem::knapsack);

#endif
