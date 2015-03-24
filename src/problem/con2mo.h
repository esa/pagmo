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

#ifndef PAGMO_PROBLEM_CON2MO_H
#define PAGMO_PROBLEM_CON2MO_H

#include <string>

#include "../serialization.h"
#include "../types.h"
#include "cec2006.h"
#include "base_meta.h"

namespace pagmo{ namespace problem {

/// Constrainted to Multi-Objective meta-problem
/**
 * Implements a meta-problem class that wraps some other constrained problems,
 * resulting in multi-objective problem.
 *
 * Three implementations of the constrained to multi-objective are available. For a problem with m constraints,
 * m+nobj objective functions, the first objectives functions are the original objective functions.
 * The first implementation is the constrained to multi-objective defined by Coello Coello. The
 * objectives defined from constraints includes number of violated constraints and objective functions.
 * The second implementation is the COMOGA multi-objective problem: a nobj+1 problem with the last
 * objective the sum of the violations of the constraints.
 * The third implementation is the same as the second one but splitting the sum of violations between equality
 * and inequality constraints, resulting in a total of nobj+2 objectives problem.
 *
 * Note: This constraints handling technique can only be used for <b>MINIMIZATION</b> problems.
 *
 * @see Coello Coello, C. A. (2002). Theoretical and numerical constraint-handling techniques used with evolutionary algorithms: a survey of the state of the art. Computer methods in applied mechanics and engineering, 191(11), 1245-1287.
 *
 * @see Coello, C. A. C. (2000). Treating constraints as objectives for single-objective evolutionary optimization. Engineering Optimization+ A35, 32(3), 275-308.
 *
 * @author Jeremie Labroquere (jeremie.labroquere@gmail.com)
 */

class __PAGMO_VISIBLE con2mo : public base_meta
{
	public:
		/// Mechanism used to deal with constraints in the objectives
		/**
		* Definition of three types of constrained to multi-objective.
		* OBJ_CSTRS is the approach suggest by Coello: the single objective constrained problem is transformed into
		* a multiobjective unconstrained problem with m+1 objectives functions (where m is the number of constraints)
		* OBJ_CSTRSVIO: the single objective constrained problem is transformed into a nobj+1 objectives unconstrained 
		* problem with the original fitness functions as first objectives and the aggregation of the constraints violation 
		* as last objective
		* OBJ_EQVIO_INEQVIO: the single objective constrained problem is transformed into a nobj+2 objectives unconstrained 
		* problem with the original fitness functions as first objectives, the aggregation of the inequality constraints violations 
		* as second last objective and the sum of violation of the equality constraints violations as last objective.
		*/
		enum method_type {
			 OBJ_CSTRS = 0, ///< Each constraint violation is transformed into one objective
			 OBJ_CSTRSVIO = 1, ///< The total constraint violation is addd as one objective
			 OBJ_EQVIO_INEQVIO = 2 ///< The total constraint violation is addd as two objectives (equalities + inequalities)
			};

		//constructors
		con2mo(const base & = cec2006(4), const method_type = OBJ_CSTRS);

		base_ptr clone() const;
		std::string get_name() const;

	protected:
		std::string human_readable_extra() const;
		void objfun_impl(fitness_vector &, const decision_vector &) const;

	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base_meta>(*this);
			ar & const_cast<method_type &>(m_method);
		}
		const method_type m_method;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::con2mo)

#endif // PAGMO_PROBLEM_CON2MO_H
