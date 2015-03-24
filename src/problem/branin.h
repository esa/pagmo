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

#ifndef PAGMO_PROBLEM_BRANIN_H
#define PAGMO_PROBLEM_BRANIN_H

#include <string>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "base.h"

namespace pagmo { namespace problem {

/// Branin rcos test function.
/**
 * \image html branin.png "Branin's rcos function."
 * \image latex branin.png "Branin's rcos function." width=5cm
 *
 * The Branin rcos function is a global optimization test function defined as:
 * \f[
 * 	f_{br}\left(x_1,x_2\right) = a \left( x_2 - bx_1^2 +cx_1 - d \right)^2 + e\left( 1 - f \right) \cos x_1 + e,
 * \f]
 * with
 * \f[
 * 	\begin{array}{rcl}
 * 		a & = & 1, \\
 * 		b & = & \frac{5.1}{4\pi^2}, \\
 * 		c & = & \frac{5}{\pi}, \\
 * 		d & = & 6, \\
 * 		e & = & 10, \\
 * 		f & = & \frac{1}{8\pi},
 * 	\end{array}
 * \f]
 * and
 * \f[
 * 	\begin{array}{rcl}
 * 		x_1 & \in & \left[ -5,10 \right ],\\
 * 		x_2 & \in & \left[ 0,15 \right ].
 * 	\end{array}
 * \f]
 * The function has three global minima, located at \f$ \left( -\pi, 12.275 \right) \f$, \f$ \left( \pi, 2.275 \right) \f$ and
 * \f$ \left( 9.42478, 2.475 \right) \f$, where the value of the function is \f$ 0.397887 \f$.
 *
 * @see Branin, F. K.: A widely convergent method for finding multiple solutions of simultaneous nonlinear equations. IBM J. Res. Develop., pp. 504-522, Sept., 1972. 
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE branin: public base
{
	public:
		branin();
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

} }

BOOST_CLASS_EXPORT_KEY(pagmo::problem::branin)

#endif
