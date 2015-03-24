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


#include "hv4d.h"

namespace pagmo { namespace util { namespace hv_algorithm {

/// Compute hypervolume
/**
 * @param[in] points vector of points containing the D-dimensional points for which we compute the hypervolume
 * @param[in] r_point reference point for the points
 *
 * @return hypervolume.
 */
double hv4d::compute(std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	// Prepare the initial data to suit the original code
	double* data = new double[points.size() * 4];
	double refpoint[4];
	for (unsigned int d_idx = 0 ; d_idx < 4 ; ++d_idx) {
		refpoint[d_idx] = r_point[d_idx];
	}
	unsigned int data_idx = 0;
	for (unsigned int p_idx = 0 ; p_idx < points.size() ; ++p_idx) {
		for (unsigned int d_idx = 0 ; d_idx < 4 ; ++d_idx) {
			data[data_idx++] = points[p_idx][d_idx];
		}
	}

	double hv = guerreiro_hv4d(data, points.size(), refpoint);
	delete[] data;
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
void hv4d::verify_before_compute(const std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	if (r_point.size() != 4) {
		pagmo_throw(value_error, "Algorithm HV4D works only for 4-dimensional cases");
	}
	base::assert_minimisation(points, r_point);
}

/// Clone method.
base_ptr hv4d::clone() const
{
	return base_ptr(new hv4d(*this));
}

/// Algorithm name
std::string hv4d::get_name() const
{
	return "Four-dimensional hypervolume by Andreia P. Guerreiro";
}

} } }

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::util::hv_algorithm::hv4d)
