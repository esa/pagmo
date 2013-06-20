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

#include "hypervolumes.h"

namespace pagmo { namespace util {

/// Computes hypervolume indicator for given pareto set using the LebMeasure algorithm.
/**
 * @see L. While, "A new analysis of the LebMeasure Algorithm for Calculating Hypervolume"
 *
 * @param[in] r reference point of the hypervolume
 * @param[in] point_set points describing the hypervolume
 *
 * @throws value_error if sizes of reference point and pareto set do not match, when the pareto set is empty and when the dimension is lesser than 3
 *
 * @return hypervolume of the pareto set.
 */
double hypervolumes::lebmeasure(const std::vector<fitness_vector> &points, const fitness_vector &r)
{
	if ( points.size() == 0 ) {
		pagmo_throw(value_error, "Point set cannot be empty.");
	}
	if ( r.size() < 3 ) {
		pagmo_throw(value_error, "Hypervolume of dimension lesser than 3 is not allowed for this method");
	}
	for (std::vector<fitness_vector>::size_type idx = 0 ; idx < points.size() ; ++idx) {
		if ( points[idx].size() != r.size() ) {
			pagmo_throw(value_error, "Point set dimensions and reference point dimension must be equal.");
		}
	}

	fitness_vector::size_type f_dim = points[0].size();

	//prepare the original points (lebmeasure requires to track the the spawn dimension corresponding to each of the points)
	lebmeasure_points point_set;
	for (population::size_type idx = 0 ; idx < points.size() ; ++idx) {
		point_set.push_back(std::make_pair(points[idx], f_dim));
	}
	return lebmeasure::compute_hypervolume(point_set, r);
}

}}
