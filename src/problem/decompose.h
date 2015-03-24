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

#ifndef PAGMO_PROBLEM_DECOMPOSE_H
#define PAGMO_PROBLEM_DECOMPOSE_H

#include <string>

#include "../serialization.h"
#include "../types.h"
#include "base_meta.h"
#include "base.h"
#include "zdt.h"

namespace pagmo{ namespace problem {


/// Decompose meta-problem
/**
 * Implements a meta-problem class resulting in a decomposed version
 * of the multi-objective input problem, i.e. a single-objective problem
 * having as fitness function some kind of combination of the original fitness functions
 *
 * Being
 * \f$ F(X) = (F_1(X), \ldots, F_n(X)) \f$
 * the original multi-objective fitness function and
 * \f$ w = (w_1, \ldots, w_n) \f$
 * \f$ z = (z_1, \ldots, z_n) \f$
 * respectively a weight vector and a reference point,
 * the decomposed problem has as single-objective fitness function one of the following
 * according to which decomposition methods is choosed:
 *
 * WEIGHTED: \f$ F_d(X) = \sum_{i=1}^n w_i F_i(X) \f$
 *
 * TCHEBYCHEFF \f$ F_d(X) = max_{1 \leq i \leq m} w_i \vert F_i(X) - z_i \vert   \f$
 *
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 * @see "Q. Zhang -- MOEA/D: A Multiobjective Evolutionary Algorithm Based on Decomposition"
 */

class __PAGMO_VISIBLE decompose : public base_meta
{
	public:
		/// Mechanism used to perform the problem decomposition
		enum method_type {
		 WEIGHTED=0, ///< The fitness function is the weighted sum of the multiple original fitnesses
		 TCHEBYCHEFF=1, ///< The Tchebycheff method is used to perform the decomposition
		 BI=2 ///< The Boundary Intersection method is used to perform the decomposition
		};

		decompose(const base & = zdt(1,2),
				  method_type = WEIGHTED,
				  const std::vector<double> & = std::vector<double>(),
				  const std::vector<double> & = std::vector<double>(),
				  const bool = false);
		base_ptr clone() const;
		std::string get_name() const;
		const std::vector<double>& get_weights() const;
		void compute_decomposed_fitness(fitness_vector &, const fitness_vector &) const;
		void compute_decomposed_fitness(fitness_vector &, const fitness_vector &, const fitness_vector &) const;
		void compute_original_fitness(fitness_vector &, const decision_vector &) const;
		fitness_vector get_ideal_point() const;
		void set_ideal_point(const fitness_vector &f);


	protected:
		std::string human_readable_extra() const;
		void objfun_impl(fitness_vector &, const decision_vector &) const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base_meta>(*this);
			ar & m_method;
			ar & m_weights;
			ar & m_z;
			ar & const_cast<bool&>(m_adapt_ideal);
		}
		method_type m_method;
		fitness_vector m_weights;
		mutable fitness_vector m_z;
		const bool m_adapt_ideal;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::decompose)

#endif // PAGMO_PROBLEM_DECOMPOSE_H
