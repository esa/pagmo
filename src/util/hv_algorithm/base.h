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

#ifndef PAGMO_UTIL_HV_ALGORITHM_BASE_H
#define PAGMO_UTIL_HV_ALGORITHM_BASE_H

#include <iostream>
#include <string>
#include <typeinfo>

#include "../../config.h"
#include "../../exceptions.h"
#include "../../serialization.h"
#include "../../types.h"

namespace pagmo { namespace util {
/// Hypervolume algorithm namespace.
/**
 * This namespace contains all the algorithms implemented for the purpose of calculating the hypervolume indicator.
 */
namespace hv_algorithm {

/// Base hypervolume algorithm class.
class base;

/// Base hypervolume algorithm class
typedef boost::shared_ptr<base> base_ptr;

/// Base hypervolume class.
/**
 * This class represents the abstract hypervolume algorithm used for computing
 * the hypervolume indicator (also known as lebesgue measure, or S-metric).
 *
 * Every hypervolume algorithm that extends this base class implement the following methods:
 * - compute() which takes a list of points and a reference point.
 * - verify_before_compute() which raises an exception for special cases of misue of this method (e.g. method working only for some number of dimensions).
 *
 * @author Krzysztof Nowak (kn@kiryx.net)
 */
class __PAGMO_VISIBLE base
{
	public:
		/// Compute method
		/**
		 * This method computes the hypervolume.
		 * It accepts a list of points as an input, and the distinguished "reference point".
		 * Hypervolume is then computed as a joint hypervolume of hypercubes, generated pairwise with the reference point and each point from the set.
		 *
		 * @param[in] points - vector of fitness_vectors for which the hypervolume is computed
		 * @param[in] r_point - distinguished "reference point".
		 */
		virtual double compute(const std::vector<fitness_vector> & points, const fitness_vector & r_point) = 0;

		/// Verification of input
		/**
		 * This method serves as a verification method.
		 * Not every algorithm is suited of every type of problem.
		 *
		 * @param[in] points - vector of fitness_vectors for which the hypervolume is computed
		 * @param[in] r_point - distringuished "reference point".
		 */
		virtual void verify_before_compute(const std::vector<fitness_vector> & points, const fitness_vector & r_point) = 0;

		/// Clone method.
		/**
		 * @return pagmo::util::hv_algorithm::base_ptr to a copy of this.
		 */
		virtual base_ptr clone() const = 0;
		virtual ~base();
	protected:
		void assert_maximal_reference_point(const std::vector<fitness_vector> &points, const fitness_vector &r_point);
};

} } }


#endif
