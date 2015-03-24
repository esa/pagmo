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

#ifndef PAGMO_PROBLEM_SCALED_H
#define PAGMO_PROBLEM_SCALED_H

#include <string>

#include "../serialization.h"
#include "ackley.h"
#include "../types.h"
#include "base_meta.h"

namespace pagmo{ namespace problem {

/// Scaled meta-problem
/**
 * Scales the fitness of the input problem.
 *
 * @author Dario Izzo (dario.izzo@gmail.com)
 */

class __PAGMO_VISIBLE scaled : public base_meta
{
	public:
		//constructors
		scaled(const base & = ackley(1), const fitness_vector &units = std::vector<double>(1));
		
		base_ptr clone() const;
		std::string get_name() const;
		
		fitness_vector descale(const fitness_vector &) const;
		const fitness_vector& get_units() const;
		
	protected:
		std::string human_readable_extra() const;
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		void compute_constraints_impl(constraint_vector &, const decision_vector &) const;
	private:

		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base_meta>(*this);
			ar & const_cast<decision_vector &>(m_units);
		}
		const decision_vector m_units;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::scaled)

#endif // PAGMO_PROBLEM_SCALED_H
