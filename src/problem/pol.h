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

#ifndef PAGMO_PROBLEM_POL_H
#define PAGMO_PROBLEM_POL_H

#include <string>

#include "../serialization.h"
#include "../types.h"
#include "base.h"

namespace pagmo{ namespace problem {

/// Poloni's study
/**
 *
 * This is a box-constrained continuous three-dimension multi-objecive problem.
 * \f[
 * 	F_1 \left(x\right) = 1 + (A_1 - B_1)^2 + (A_2 - B_2)^2
 * \f]
 * \f[
 *      F_2 \left(x\right) = (x_1 + 3)^2 + (x_2+1)^2  x \in \left[ -\pi,\pi \right].
 * \f]
 * \f[
 *      A_1 = 0.5 \sin(1) - 2 \cos(1) + \sin(2) -1.5 \cos(2)
 * \f]
 * \f[
 *      A_2 = 1.5 \sin(1) - \cos(1) + 2 \sin(2) - 0.5 \cos(2)
 * \f]
 * \f[
 *      B_1 = 0.5 \sin(x_1) - 2 \cos(x_1) + \sin(x_2) - 1.5 \cos(x_2)
 * \f]
 * \f[
 *      B_2 = 1.5 \sin(x_1) - \cos(x_1) + 2 \sin(x_2) - 0.5 \cos(x_2)
 * \f]
 *
 * @see http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.18.4257&rep=rep1&type=pdf
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 */

class __PAGMO_VISIBLE pol : public base
{
	public:
		pol();
		base_ptr clone() const;
		std::string get_name() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
		}
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::pol)

#endif // PAGMO_PROBLEM_POL_H
