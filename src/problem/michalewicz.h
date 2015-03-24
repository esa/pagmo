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

#ifndef PAGMO_PROBLEM_MICHALEWICZ_H
#define PAGMO_PROBLEM_MICHALEWICZ_H

#include <string>

#include "../serialization.h"
#include "../types.h"
#include "base.h"

namespace pagmo{ namespace problem {

/// The Michalewicz problem.
/**
 * \image html michalewicz.jpg "Two-dimensional Michalewicz function."
 * \image latex michalewicz.jpg "Two-dimensional Michalewicz function." width=5cm
 *
 * This is a box-constrained continuous single-objecive problem.
 * The objective function is the  n-dimensional Michalewicz function:
 * \f[
 * 	F \left(x_1,\ldots,x_n\right) = \sum_{i=1}^n sin(x_i) \left[sin \left(\frac{i x_i^2}{\pi} \right) \right]^2m, \quad x_i \in \left[ 0,\pi \right].
 * \f]
 * Several global minima, one local minimum (ex: \f$ m=10, n=5 f(x) = -4.687658\f$).
 *
 * @see http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2376.htm
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 */

class __PAGMO_VISIBLE michalewicz : public base
{
	public:
		michalewicz(int n = 1,int m = 10);
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
			ar & m_m;
		}
		int m_m;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::michalewicz)

#endif // PAGMO_PROBLEM_MICHALEWICZ_H
