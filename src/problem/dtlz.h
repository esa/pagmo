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

#ifndef PAGMO_PROBLEM_DTLZ_H
#define PAGMO_PROBLEM_DTLZ_H

#include <string>

#include "../serialization.h"
#include "../types.h"
#include "base_dtlz.h"

namespace pagmo{ namespace problem {

/// DTLZ problem test suite
/**
 * All problems in this test suite are box-constrained continuous n-dimensional multi-objective problems, scalable in fitness dimension.
 * The dimension of the decision space is \f$ k + fdim - 1 \f$, whereas fdim is the number of objectives and k a paramter.
 * Properties of the decision space and the Pareto-front of each problems are as follow:
 * 
 * DTLZ1:
 *
 * The optimal pareto front lies on a linear hyperplane \f$ \sum_{m=1}^M f_m = 0.5 \f$ .
 *
 * DTLZ2:
 *
 * The search space is continous, unimodal and the problem is not deceptive. 
 *
 * DTLZ3:
 *
 * The search space is continous, unimodal and the problem is not deceptive.
 * It is supposed to be harder to converge towards the optimal pareto front than DTLZ2
 * 
 * DTLZ4:
 * 
 * The search space contains a dense area of solutions next to the f_M/f_1 plane.
 * 
 * DTLZ5:
 * 
 * This problem will test an MOEA's ability to converge to a cruve and will also allow an easier way to visually demonstrate
 * (just by plotting f_M with any other objective function) the performance of an MOEA. Since there is a natural bias for
 * solutions close to this Pareto-optimal curve, this problem may be easy for an algorithmn to solve. Because of its simplicity
 * its recommended to use a higher number of objectives \f$ M \in [5, 10]\f$.
 * 
 * DTLZ6:
 *
 * A more difficult version of the DTLZ5 problem: the non-linear distance function g makes it harder to convergence
 * against the pareto optimal curve.
 * 
 * DTLZ7:
 * 
 * This problem has disconnected Pareto-optimal regions in the search space.
 * 
 * @see K. Deb, L. Thiele, M. Laumanns, E. Zitzler, Scalable test problems for evoulationary multiobjective optimization
 * @author Marcus Maertens (mmarcusx@gmail.com)
 * 
 */

class __PAGMO_VISIBLE dtlz : public base_dtlz
{
	public:
		dtlz(size_type id = 1, size_type k = 5, fitness_vector::size_type fdim = 3, const size_t alpha = 100);
		base_ptr clone() const;
		std::string get_name() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
	private:
                void f1_objfun_impl(fitness_vector &, const decision_vector &) const;
                void f23_objfun_impl(fitness_vector &, const decision_vector &) const;
                void f4_objfun_impl(fitness_vector &, const decision_vector &) const;
                void f56_objfun_impl(fitness_vector &, const decision_vector &) const;
                void f7_objfun_impl(fitness_vector &, const decision_vector &) const;
		double g13_func(const decision_vector &) const;
		double g245_func(const decision_vector &) const;
		double g6_func(const decision_vector &) const;
		double g7_func(const decision_vector &) const;
		double h7_func(const fitness_vector &, const double) const;	// used only by DTLZ7

		double g_func(const decision_vector &) const;	// general g-function for wrapping switch
                
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base_dtlz>(*this);
                        ar & const_cast<unsigned int&>(m_problem_number);
			ar & const_cast<double &>(m_alpha);
                }

                const unsigned int m_problem_number;
		const double m_alpha; 		// used only for DTLZ4
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::dtlz)

#endif // PAGMO_PROBLEM_DTLZ_H
