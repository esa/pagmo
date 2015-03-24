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

#ifndef PAGMO_PROBLEM_GTOC_1_H
#define PAGMO_PROBLEM_GTOC_1_H

#include <string>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "../AstroToolbox/mga.h"
#include "base.h"

namespace pagmo{ namespace problem {

/// GTOC_1 MGA problem
/**
 * This problem draws inspiration from the first edition of the Global Trajectory Optimisation Competition (GTOC1)
 * It is an 8 dimensional MGA problem with a rather long fly-by sequence including mainly Earth and Venus. The final target
 * is the asteroid TW229. The objective of the mission is to maximise the change in sami-major axis of the
 * asteroid orbit following an anaelastic impact of the spacecraft with the asteroid.
 *
 * The problem is also part of the Global Trajectory Optimization database (GTOP)
 *
 * gtoc_1 is a box constrained single objective, continuous optimization problem of dimension eight.
 *
 * @see http://www.esa.int/gsp/ACT/inf/op/globopt/evevejsa.htm
 * @author Dario Izzo (dario.izzo@esa.int)
 */
class __PAGMO_VISIBLE gtoc_1: public base
{
	public:
		gtoc_1();
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
			ar & problem;
			ar & Delta_V;
			ar & rp;
			ar & t;
		}
		mgaproblem problem;
		//Variables used in the call to MGA
		mutable std::vector<double> Delta_V;
		mutable std::vector<double> rp;
		mutable std::vector<double> t;
};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::problem::gtoc_1)

#endif // PAGMO_PROBLEM_GTOC_1_H
