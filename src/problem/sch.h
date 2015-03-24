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

#ifndef PAGMO_PROBLEM_SCH_H
#define PAGMO_PROBLEM_SCH_H

#include <string>

#include "../serialization.h"
#include "../types.h"
#include "base.h"

namespace pagmo{ namespace problem {

/// Shaffer's study problem
/**
 *
 * This is a box-constrained continuous one dimensional multi-objective problem.
 * \f[
 * 	F_1 \left( x\right) = x^2
 * \f]
 * \f[
 * 	F_2 \left( x\right) = (x-2)^2  x \in \left[ -10^3,10^3 \right].
 * \f]
 * The optimal solutions (minimization problem) are \f$ x \in \left[0,2\right]\f$
 *
 * @see http://portal.acm.org/citation.cfm?id=657079
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 */

class __PAGMO_VISIBLE sch : public base
{
	public:
		sch();
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

BOOST_CLASS_EXPORT_KEY(pagmo::problem::sch)

#endif // PAGMO_PROBLEM_SCH_H
