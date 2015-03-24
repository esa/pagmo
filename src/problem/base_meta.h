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

#ifndef PAGMO_PROBLEM_BASE_META_H
#define PAGMO_PROBLEM_BASE_META_H

#include <string>

#include "../serialization.h"
#include "ackley.h"
#include "../types.h"
#include "base.h"

namespace pagmo{ namespace problem {

/// Meta=problems base class
/**
 * All meta-problems can inherit directly from this class. The virtual functions
 * to compare fitness and constraints vectors are thus passed to the meta-problem.
 * The copy constructor is also implemented as to provide a deep copy of the object.
 *
 * NOTE: The original problem can also have the virtual function compare_fitness_impl implemented
 * the meta-problem will make use of it. The custom implementation, in this case, is expected to work
 * for generic dimensions of the fitness vector as metaproblems may transform this dimension at will.
 *
 * @author Dario Izzo (dario,izzo@gmail.com)
 */

class __PAGMO_VISIBLE base_meta : public base
{
	public:
		/// Constructor
		base_meta(const base &p = ackley(1), int n=1, int ni=0, int nf=1, int nc=0, int nic=0, const std::vector<double>&c_tol = std::vector<double>()):
			 base(n,ni,nf,nc,nic,c_tol), m_original_problem(p.clone()) {
			 	//Setting the bounds according to the original problem
				set_bounds(m_original_problem->get_lb(),m_original_problem->get_ub());
			 }
		/// Copy constructor
		base_meta(const base_meta &p):base(p), m_original_problem(p.m_original_problem->clone()) {}
	protected:
		bool compare_fitness_impl(const fitness_vector &f1, const fitness_vector &f2) const 
			{return m_original_problem->compare_fitness_impl(f1,f2);}
		//NOTE: It is not possible to use the same trick for the other two virtual compares as they also depend from
		//the class parameter m_ic_dimension which in m_original_problem may be different than in the meta problem
		bool compare_constraints_impl(const constraint_vector &c1, const constraint_vector &c2) const
			{return m_original_problem->compare_constraints_impl(c1,c2);}
		bool compare_fc_impl(const fitness_vector &f1, const constraint_vector &c1, const fitness_vector &f2, const constraint_vector &c2) const
			{return m_original_problem->compare_fc_impl(f1,c1,f2,c2);}
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & m_original_problem;
		}
	protected:
		/// Smart pointer to the original problem instance
		base_ptr m_original_problem;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::base_meta)

#endif // PAGMO_PROBLEM_BASE_META_H
