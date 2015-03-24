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

#ifndef PAGMO_PROBLEM_BASE_UNC_MO_H
#define PAGMO_PROBLEM_BASE_UNC_MO_H

#include "base.h"
#include "../serialization.h"

namespace pagmo{ namespace problem {

/// Base Class for Unconstrained Multi-objective problems
/**
 * This class adds to pagmo::problem::base the virtual method p_distance
 * which returns the convergence metric of an individual or a population if
 * the user re-implemented the convergence_metric virtual method.
 * Only unconstrained multi-objective problems can derive from this class.
 * 
 * In problems where it is possible, the user needs to derive from
 * pagmo::problem::base_unc_mo and reimplement the virtual method
 * double base_unc_mo::convergence_metric(const decision_vector &x)
 * (typically returning 0 if x lies on the Pareto Front)
 *
 * @author Dario Izzo (dario.izzo@gmail.com)
 */
class __PAGMO_VISIBLE base_unc_mo : public base
{
	public:
		base_unc_mo(base::size_type, base::size_type, base::f_size_type);
		
		// p_distance methods
		double p_distance(const decision_vector &) const;
		double p_distance(const pagmo::population &) const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
		}
	protected:
		virtual double convergence_metric(const decision_vector &) const;

};

}} //namespaces

BOOST_SERIALIZATION_ASSUME_ABSTRACT(pagmo::problem::base_unc_mo)

#endif //PAGMO_PROBLEM_BASE_UNC_MO_H
