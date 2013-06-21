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

#include "optimal2d.h"

namespace pagmo { namespace util {

/// Computes hypervolume indicator for given pareto set using the optimal 2D algorithm.
/**
 * This method should be used both as a solution to 2D cases, and as a general termination method for algorithms that reduce n-dimensional problem to 2D.
 *
 * Computational complexity: n*log(n)
 *
 * @param[in] r reference point of the hypervolume
 * @param[in] vector of points containing the pareto set
 *
 * @return hypervolume of the pareto set.
 */
double optimal2d::compute_hypervolume(std::vector<fitness_vector> points, const fitness_vector &r)
{
	if (r.size() != 2) {
		pagmo_throw(value_error, "Optimal 2D hypervolume method works only for 2-dimensional cases.");
	}
	sort(points.begin(), points.end(), compare_fitness);
	double hypervolume = 0.0;
	for( std::vector<fitness_vector>::size_type idx = 0; idx < points.size() - 1 ; ++idx) {
		double area = (points[idx][0] - points[idx+1][0]) * (points[idx][1] - r[1]);
		hypervolume += fabs(area);
	}
	fitness_vector &last = points.back();
	hypervolume += fabs((r[0] - last[0]) * (r[1] - last[1]));

	return hypervolume;
}

inline bool compare_fitness(const fitness_vector &a, const fitness_vector &b) {
	return a[0] < b[0];
}

}}
