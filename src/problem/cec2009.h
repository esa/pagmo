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

#ifndef PAGMO_PROBLEM_CEC2009_H
#define PAGMO_PROBLEM_CEC2009_H

#include <string>

#include "../serialization.h"
#include "../types.h"
#include "base.h"

namespace pagmo{ namespace problem {

/// The CEC 2009 problems: Competition on "Performance Assessment of Constrained / Bound
///  Constrained Multi-Objective Optimization Algorithms"
/**
 *
 * This class instantiates any of the problems from CEC2009's competition
 * on multi-objective optimization algorithms, commonly referred to by the literature
 * as UF1-UF10 (unconstrained) and CF1-CF10 (constrained).
 *
 * Note: The three problems constructed by some transformation on DTLZ2, DTLZ3
 * and WFG1 problems as described in the technical report are not included in
 * this implementation.
 *
 * @see http://www3.ntu.edu.sg/home/EPNSugan/index_files/CEC09-MOEA/CEC09-MOEA.htm
 *
 * @author Yung-Siang Liau (liauys@gmail.com)
 */

class __PAGMO_VISIBLE cec2009 : public base
{
	public:
		cec2009(unsigned int = 1, problem::base::size_type = 30, bool = false);
		base_ptr clone() const;
		std::string get_name() const;

	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		void compute_constraints_impl(constraint_vector &, const decision_vector &) const;

	private:
		void UF1(const double *x, double *f, const unsigned int nx) const;
		void UF2(const double *x, double *f, const unsigned int nx) const;
		void UF3(const double *x, double *f, const unsigned int nx) const;
		void UF4(const double *x, double *f, const unsigned int nx) const;
		void UF5(const double *x, double *f, const unsigned int nx) const;
		void UF6(const double *x, double *f, const unsigned int nx) const;
		void UF7(const double *x, double *f, const unsigned int nx) const;
		void UF8(const double *x, double *f, const unsigned int nx) const;
		void UF9(const double *x, double *f, const unsigned int nx) const;
		void UF10(const double *x, double *f, const unsigned int nx) const;
		
		void CF1(const double *x, double *f, double *c, const unsigned int nx) const;
		void CF2(const double *x, double *f, double *c, const unsigned int nx) const;
		void CF3(const double *x, double *f, double *c, const unsigned int nx) const;
		void CF4(const double *x, double *f, double *c, const unsigned int nx) const;
		void CF5(const double *x, double *f, double *c, const unsigned int nx) const;
		void CF6(const double *x, double *f, double *c, const unsigned int nx) const;
		void CF7(const double *x, double *f, double *c, const unsigned int nx) const;
		void CF8(const double *x, double *f, double *c, const unsigned int nx) const;
		void CF9(const double *x, double *f, double *c, const unsigned int nx) const;
		void CF10(const double *x, double *f, double *c, const unsigned int nx) const;

		static fitness_vector::size_type cec2009_fitness_dimension(int);
		static constraint_vector::size_type cec2009_ic_dimension(int);

		void configure_bounds();

		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & const_cast<unsigned int&>(m_problem_number);
			ar & const_cast<bool&>(m_is_constrained);
		}

		const unsigned int m_problem_number;
		const bool m_is_constrained;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::cec2009)

#endif
