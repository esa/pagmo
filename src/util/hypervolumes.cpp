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

/// Computes Hypervolume Indicator by LebMeasure
/**
 * Computes hypervolume indicator for given pareto set using the LebMeasure algorithm.
 *
 * @see L. While, "A new analysis of the LebMeasure Algorithm for Calculating Hypervolume"
 *
 * @param[in] r reference point of the hypervolume
 * @param[in] point_set points describing the hypervolume
 *
 * @return hypervolume of the pareto set.
 */
double hypervolumes::hypervolume_lebmeasure(const std::vector<fitness_vector> &points, const fitness_vector &r)
{
	//std::vector<std::vector<population::size_type> > pareto_fronts = compute_pareto_fronts();
	fitness_vector::size_type f_dim = points[0].size();

	std::deque<std::pair<fitness_vector, fitness_vector::size_type> > point_set;
	for (population::size_type idx = 0 ; idx < points.size() ; ++idx) {
		point_set.push_back(std::make_pair(points[idx], f_dim));
	}
	return lebmeasure::hypervolume_lebmeasure(point_set, r);
}

}}
