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

#ifndef PAGMO_PROBLEM_CASSINI_2_H
#define PAGMO_PROBLEM_CASSINI_2_H

#include <string>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "../AstroToolbox/mga_dsm.h"
#include "base.h"

namespace pagmo{ namespace problem {

/// Cassini MGA-DSM Problem
/**
 *
 * The problem is transcribed as an MGA-DSM problem allowing one chemical manouvre per trajectory leg.
 * The objective function is defined
 *
 * cassini_2 is a box constrained single objective, continuous optimization problem of dimension 22.
 *
 * @see http://www.esa.int/gsp/ACT/inf/op/globopt/TandEM.htm
 * @author Dario Izzo (dario.izzo@esa.int)
 */
class __PAGMO_VISIBLE cassini_2: public base
{
	public:
		cassini_2();
		base_ptr clone() const;
		std::string pretty(const std::vector<double> &x) const;
		std::string get_name() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		void set_sparsity(int &, std::vector<int> &, std::vector<int> &) const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			//[DS] "sequence" array is static constant so it doesn't need serialization
			ar & problem;
		}
		static const int sequence[6];
		mgadsmproblem problem;

};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::problem::cassini_2)

#endif // PAGMO_PROBLEM_CASSINI_2_H
