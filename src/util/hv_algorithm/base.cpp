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

#include "base.h"


namespace pagmo { namespace util { namespace hv_algorithm {

// Destructor
base::~base() {}

/// Assert that reference point dominates every other point from the set.
/**
 * This is a method that can be referenced from verify_before_compute method.
 * The method checks whether provided reference point is "no worse" and in at least one objective "better" for each of the points from the set.
 *
 * @param[in] points - vector of fitness_vectors for which the hypervolume is computed
 * @param[in] r_point - distringuished "reference point".
*/
void base::assert_maximal_reference_point(const std::vector<fitness_vector> &points, const fitness_vector &r_point) {
	for(std::vector<fitness_vector>::size_type idx = 0 ; idx < points.size() ; ++idx) {
		bool outside_bounds = false;
		bool all_equal = true;

		for(fitness_vector::size_type f_idx = 0 ; f_idx < points[idx].size() ; ++f_idx) {
			outside_bounds |= (r_point[f_idx] < points[idx][f_idx]);
			all_equal &= (r_point[f_idx] == points[idx][f_idx]);
		}
		if (all_equal || outside_bounds) {
			std::stringstream ss;
			std::string str_p("("), str_r("(");
			for(fitness_vector::size_type f_idx = 0 ; f_idx < points[idx].size() ; ++f_idx) {
					str_p += boost::lexical_cast<std::string>(points[idx][f_idx]);
					str_r += boost::lexical_cast<std::string>(r_point[f_idx]);
				if (f_idx < points[idx].size() - 1) {
					str_p += ", ";
					str_r += ", ";
				} else {
					str_p += ")";
					str_r += ")";
				}
			}
			ss << "Reference point is invalid: another point seems to be outside the reference point boundary, or be equal to it:" << std::endl;
			ss << " P[" << idx << "]\t= " << str_p << std::endl;
			ss << " R\t= " << str_r << std::endl;
			pagmo_throw(value_error, ss.str());
		}
	}
}

} } }
