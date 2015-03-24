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

#ifndef LUKSAN_VLCEK_2_H
#define LUKSAN_VLCEK_2_H

#include <string>
#include <vector>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "base.h"

namespace pagmo { namespace problem {

/// Test problem from the Luksan and Vlcek book.
/**
 * Implementation of the Example 5.2 in "Sparse and Parially Separable
 * Test Problems for Unconstrained and Equality Constrained
 * Optimization" by Luksan and Vlcek. Code adapted from the Ipopt scalable problem examples
 *
 * @author Dario Izzo (dario.izzo@esa.int)
 */

class __PAGMO_VISIBLE luksan_vlcek_2: public base
{
	public:
		luksan_vlcek_2(int = 16, const double & = 0, const double & = 0);
		base_ptr clone() const;
		std::string get_name() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		void compute_constraints_impl(constraint_vector &, const decision_vector &) const;
		void set_sparsity(int &, std::vector<int> &, std::vector<int> &) const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & m_clb;
			ar & m_cub;
		}
		std::vector<double>	m_clb;
		std::vector<double>	m_cub;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::luksan_vlcek_2)

#endif // LUKSAN_VLCEK_2_H
