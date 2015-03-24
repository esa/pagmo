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

#ifndef PAGMO_PROBLEM_BASE_DTLZ_H
#define PAGMO_PROBLEM_BASE_DTLZ_H

#include "base_unc_mo.h"
#include "../serialization.h"

namespace pagmo{ namespace problem {

/// Base DTLZ Multi-objective optimization problem.
/**
 * The DTLZ test suite was introduced by Deb et. al in order to create benchmark problems
 * with various features that are scalable to any number of objectives and dimensions. The first
 * seven DTLZ problems are constructed with a bottom up approach: First, the pareto-optimal front is described by a surface. 
 * Then, parallel layers of this surface are used to build the remainder of a decision space. 
 * For this, a distance function g is introduced which is minimized for pareto-optimal solutions. 
 * By evaluation of the g-function we thus get information how close the a solution is to the Pareto-optimal front. 
 * This is used to define a simple convergence metric (p-distance) that is average value of the g-function over all individuals.
 * 
 * The first seven DTLZ problems are already implemented in PaGMO but you can create your own DTLZ-style problems
 * by subclassing this base class. In order to do this, you have to implement a g-function. 
 * The shape of the front itself has to be defined in the implementation
 * of the objective function. The p-distance and a generic 3d-plot will then automatically be 
 * available for your new problem.
 * 
 * @see K. Deb, L. Thiele, M. Laumanns, E. Zitzler, Scalable test problems for evoulationary multiobjective optimization
 * 
 * @author Marcus Maertens (mmarcusx@gmail.com)
 */
class __PAGMO_VISIBLE base_dtlz : public base_unc_mo
{
	public:
		base_dtlz(int, int);
		/// Clone method.
		virtual base_ptr clone() const = 0;
	protected:
		/// Distance function
		/**
		 * This pure virtual function is re-implemented in the derived classes
		 * and is used to compute the distance of a point from the Pareto front
		 **/
		virtual double g_func(const decision_vector &) const = 0;
	private:
		double convergence_metric(const decision_vector &) const;
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base_unc_mo>(*this);
		}
};

}} //namespaces

BOOST_SERIALIZATION_ASSUME_ABSTRACT(pagmo::problem::base_dtlz)

#endif
