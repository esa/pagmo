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

#ifndef PAGMO_PROBLEM_HIMMELBLAU_H
#define PAGMO_PROBLEM_HIMMELBLAU_H

#include <string>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "base.h"

namespace pagmo { namespace problem {

/// Himmelblau function minimisation problem.
/**
 * \image html himmelblau.png "Plot of Himmelblau's function."
 * \image latex himmelblau.png "Plot of Himmelblau's function." width=3cm
 *
 * The Himmelblau function is a multi-modal function, often used to test the performance of optimization algorithms. The function is defined by:
 * \f[
 * 	f(x, y) = \left(x^2+y-11\right)^2 + \left(x+y^2-7\right)^2.
 * \f]
 * The function has one local maximum at \f$ x = -0.270844 \f$ and \f$ y = -0.923038 \f$, where \f$ f(x,y) = 181.616 \f$, and four identical local minimums:
 * \f$ f(3.0, 2.0) = 0.0 \f$, \f$ f(-2.805118, 3.131312) = 0.0 \f$, \f$ f(-3.779310, -3.283186) = 0.0 \f$, \f$ f(3.584428, -1.848126) = 0.0 \f$.
 *
 * The bounds of the search space are in the \f$ [-6,6] \f$ range.
 *
 * @see http://en.wikipedia.org/wiki/Himmelblau%27s_function
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE himmelblau: public base
{
	public:
		himmelblau();
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

}}

BOOST_CLASS_EXPORT_KEY(pagmo::problem::himmelblau)

#endif
