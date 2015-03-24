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

#ifndef PAGMO_PROBLEM_ZDT_H
#define PAGMO_PROBLEM_ZDT_H

#include <string>

#include "../serialization.h"
#include "../types.h"
#include "base_unc_mo.h"

namespace pagmo{ namespace problem {

/// ZDT problem test suite
/**
 * ZDT1:
 *
 * This is a box-constrained continuous n-dimensional (n>1) multi-objecive problem.
 * \f[
 *	g\left(x\right) = 1 + 9 \left(\sum_{i=2}^{n} x_i \right) / \left( n-1 \right)
 * \f]
 * \f[
 * 	F_1 \left(x\right) = x_1
 * \f]
 * \f[
 *      F_2 \left(x\right) = g(x) \left[ 1 - \sqrt{x_1 / g(x)} \right]  x \in \left[ 0,1 \right].
 *
 * \f]
 *
 * ZDT2:
 *
 * This is a box-constrained continuous n-dimension multi-objecive problem.
 * \f[
 *      g\left(x\right) = 1 + 9 \left(\sum_{i=2}^{n} x_i \right) / \left( n-1 \right)
 * \f]
 * \f[
 *      F_1 \left(x\right) = x_1
 * \f]
 * \f[
 *      F_2 \left(x\right) = g(x) \left[ 1 - \left(x_1 / g(x)\right)^2 \right]  x \in \left[ 0,1 \right].
 *
 * \f]
 *
 * ZDT3:
 *
 * This is a box-constrained continuous n-dimension multi-objecive problem.
 * \f[
 *      g\left(x\right) = 1 + 9 \left(\sum_{i=2}^{n} x_i \right) / \left( n-1 \right)
 * \f]
 * \f[
 *      F_1 \left(x\right) = x_1
 * \f]
 * \f[
 *      F_2 \left(x\right) = g(x) \left[ 1 - \sqrt{x_1 / g(x)} - x_1/g(x) \sin(10 \pi x_1) \right]  x \in \left[ 0,1 \right].
 *
 * \f]
 *
 * ZDT4:
 *
 * This is a box-constrained continuous n-dimension multi-objecive problem.
 * \f[
 *      g\left(x\right) = 91 + \sum_{i=2}^{n} \left[x_i^2 - 10 \cos \left(4\pi x_i \right) \right]
 * \f]
 * \f[
 *      F_1 \left(x\right) = x_1
 * \f]
 * \f[
 *      F_2 \left(x\right) = g(x) \left[ 1 - \sqrt{x_1 / g(x)} \right]  x_1 \in [0,1], x_i \in \left[ -5,5 \right] i=2, \cdots, 10.
 *
 * \f]
 *
 * ZDT5
 *
 * This is a box-constrained integer n-dimension multi-objecive problem.
 * \f[
 *       F_1\left(x\right) = 1 + u \left(x_{1} \right)
 * \f]
 * \f[
 *       g\left(x\right) = \sum_{i=2}^{11} v \left(u \left(x_{i} \right) \right)
 * \f]
 * \f[
 *       v\left(u\left(x_{i}\right)\right) =  2 + u \left(x_{i} \right)    if u \left(x_{i} \right) < 5
 *       v\left(u\left(x_{i}\right)\right) =  1                            if u \left(x_{i} \right) = 5
 * \f]
 * \f[
 *       F_2 = g \left(x \right) * 1/F_1 \left(x \right)
 * \f]
 *
 * ZDT6
 *
 * This is a box-constrained continuous 30-dimension multi-objecive problem.
 * \f[
 *      g\left(x\right) = 1 + 9 \left[\left(\sum_{i=2}^{n} x_i \right) / \left( n-1 \right)\right]^{0.25}
 * \f]
 * \f[
 *      F_1 \left(x\right) = 1 - \exp(-4 x_1) \sin^6(6 \pi \ x_1)
 * \f]
 * \f[
 *      F_2 \left(x\right) = g(x) \left[ 1 - (f_1(x) / g(x))^2  \right]  x \in \left[ 0,1 \right].
 *
 * \f]
 *
 * @see = http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.30.5848&rep=rep1&type=pdf
 * @see http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.18.4257&rep=rep1&type=pdf
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 * @author Dario Izzo (dario.izzo@googlemail.com)
 * @author Jainit Purohit (mjainit@gmail.com)
 * @author Marcus Märtens (mmarcusx@gmail.com)
 *
 */

class __PAGMO_VISIBLE zdt : public base_unc_mo
{
	public:
		zdt(size_type = 1, size_type = 30);
		base_ptr clone() const;
		std::string get_name() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		double convergence_metric(const decision_vector &) const;
	private:
				void g01_objfun_impl(fitness_vector &, const decision_vector &) const;
				double g0123_convergence_metric(const decision_vector &) const;
				void g02_objfun_impl(fitness_vector &, const decision_vector &) const;
				void g03_objfun_impl(fitness_vector &, const decision_vector &) const;
				void g04_objfun_impl(fitness_vector &, const decision_vector &) const;
				double g04_convergence_metric(const decision_vector &) const;
				void g05_objfun_impl(fitness_vector &, const decision_vector &) const;
				double g05_convergence_metric(const decision_vector &) const;
				void g06_objfun_impl(fitness_vector &, const decision_vector &) const;
				double g06_convergence_metric(const decision_vector &) const;

				friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base_unc_mo>(*this);
						ar & const_cast<unsigned int&>(m_problem_number);
				}

				const unsigned int m_problem_number;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::zdt)

#endif // PAGMO_PROBLEM_ZDT_H
