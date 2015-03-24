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

#ifndef PAGMO_PROBLEM_SHIFTED_H
#define PAGMO_PROBLEM_SHIFTED_H

#include <string>

#include "../serialization.h"
#include "ackley.h"
#include "../types.h"
#include "base_meta.h"

namespace pagmo{ namespace problem {

/// Shifted meta-problem
/**
 * Implements a meta-problem class that wraps some other problems,
 * resulting in a shifted (translated) version of the underlying problem.
 *
 * @author Yung-Siang Liau (liauys@gmail.com)
 */

class __PAGMO_VISIBLE shifted : public base_meta
{
	public:
		//constructors
		shifted(const base &, const decision_vector &);
		shifted(const base &, const double);
		shifted(const base & = ackley(1));
		
		base_ptr clone() const;
		std::string get_name() const;
		
		decision_vector deshift(const decision_vector &) const;
		const decision_vector& get_shift_vector() const;
		
	protected:
		std::string human_readable_extra() const;
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		void compute_constraints_impl(constraint_vector &, const decision_vector &) const;
	private:
		void configure_shifted_bounds(const decision_vector &);

		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base_meta>(*this);
			ar & m_translation;
		}
		decision_vector m_translation;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::shifted)

#endif // PAGMO_PROBLEM_SHIFTED_H
