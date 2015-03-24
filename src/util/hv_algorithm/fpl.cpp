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


#include "fpl.h"

namespace pagmo { namespace util { namespace hv_algorithm {

/// Compute hypervolume
/**
 * @param[in] points vector of points containing the D-dimensional points for which we compute the hypervolume
 * @param[in] r_point reference point for the points
 *
 * @return hypervolume.
 */
double fpl::compute(std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	// Prepare the initial data to suit the signature of the function 'fpli_hv'
	unsigned int fdim = points[0].size();
	std::vector<double> data;
	data.reserve(points.size() * fdim);
	for (unsigned int p_idx = 0 ; p_idx < points.size() ; ++p_idx) {
		for (unsigned int d_idx = 0 ; d_idx < fdim ; ++d_idx) {
			data.push_back(points[p_idx][d_idx]);
		}
	}

	double hv = fpli_hv(&data[0], fdim, points.size(), &r_point[0]);
	return hv;
}

/// Verify before compute
/**
 * Verifies whether given algorithm suits the requested data.
 *
 * @param[in] points vector of points containing the 4-dimensional points for which we compute the hypervolume
 * @param[in] r_point reference point for the vector of points
 *
 * @throws value_error when trying to compute the hypervolume for the non-maximal reference point
 */
void fpl::verify_before_compute(const std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	base::assert_minimisation(points, r_point);
}

/// Clone method.
base_ptr fpl::clone() const
{
	return base_ptr(new fpl(*this));
}

/// Algorithm name
std::string fpl::get_name() const
{
	return "FPL algorithm";
}

} } }

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::util::hv_algorithm::fpl)
