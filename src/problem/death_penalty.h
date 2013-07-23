/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
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

#ifndef PAGMO_PROBLEM_DEATH_PENALTY_H
#define PAGMO_PROBLEM_DEATH_PENALTY_H

#include <string>

#include "../serialization.h"
#include "../types.h"
#include "cec2006.h"
#include "base.h"

namespace pagmo{ namespace problem {

/// Constrainted death penalty meta-problem
/**
 * Implements a meta-problem class that wraps some other constrained problems,
 * resulting in death penalty constraints handling.
 *
 * Two implementations of the death penalty are availlable. The first one
 * is the most common simple death penalty. The second one is the death
 * penalty defined by Angel Kuri Morales et al. Simple death penalty penalizes the fitness function with a high value, Kuri method penalizes the
 * fitness function according to the rate of satisfied constraints.
 *
 * Note: This constraints handling technique can only be used for <b>MINIMIZATION</b> problems.
 *
 * @see Coello Coello, C. A. (2002). Theoretical and numerical constraint-handling techniques used with evolutionary algorithms: a survey of the state of the art. Computer methods in applied mechanics and engineering, 191(11), 1245-1287.
 * @see Kuri Morales, A. and Quezada, C.C. A Universal eclectic genetic algorithm for constrained optimization, Proceedings 6th European Congress on Intelligent Techniques & Soft Computing, EUFIT'98, 518-522, 1998. 
 *
 * @author Jeremie Labroquere (jeremie.labroquere@gmail.com)
 */

class __PAGMO_VISIBLE death_penalty : public base
{
	public:
		/// Type of death penalty.
		/**
		* Definition of two types of death penalty simple and kuri.
		* Simple death penalty penalizes the fitness function with a high value, Kuri method penalizes the
		* fitness function according to the rate of satisfied constraints.
		*/
		//death penalty type simple or kuri
		enum method_type {SIMPLE = 0, KURI = 1};

		//constructors
		death_penalty(const base & = cec2006(4), const method_type = SIMPLE);

		//copy constructor
		death_penalty(const death_penalty &);
		base_ptr clone() const;
		std::string get_name() const;

	protected:
		std::string human_readable_extra() const;
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		bool compare_fitness_impl(const fitness_vector &v_f1, const fitness_vector &v_f2) const;

	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & m_original_problem;
			ar & const_cast<method_type &>(m_method);
		}
		base_ptr m_original_problem;

		const method_type m_method;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::death_penalty);

#endif // PAGMO_PROBLEM_DEATH_PENALTY_H
