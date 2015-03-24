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


#include "hoy.h"
#include "base.h"
#include <algorithm>
#include <bitset>
#include <cmath>
#include <limits>
#include <boost/numeric/conversion/cast.hpp>

namespace pagmo { namespace util { namespace hv_algorithm {

/// Constructor
hoy::hoy() { }

/// Compute hypervolume
/**
 * @param[in] points vector of points containing the D-dimensional points for which we compute the hypervolume
 * @param[in] r_point reference point for the points
 *
 * @return hypervolume.
 */
double hoy::compute(std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	m_dimension = r_point.size();
	m_total_size = points.size();
	m_sqrt_size = sqrt(boost::numeric_cast<double>(m_total_size));
	m_volume = 0.0;

	sort(points.begin(), points.end(), fitness_vector_cmp(m_dimension - 1, '<'));

	m_region_low = new double[m_dimension - 1];
	m_region_up = new double[m_dimension - 1];
	m_boundaries = new double[m_total_size];
	m_no_boundaries = new double[m_total_size];
	m_piles = new int[m_total_size];
	m_trellis = new double[m_dimension - 1];

	// initialize the D-1 dimensional region vectors and D-dimensional reference point
	for (int i = 0 ; i < m_dimension - 1 ; ++i) {
		m_region_up[i] = r_point[i];
		m_region_low[i] = std::numeric_limits<double>::max();
	}

	double** initial_points = new double*[m_total_size];
	for (int n = 0 ; n < m_total_size ; ++n) {
		initial_points[n] = new double[m_dimension];

		initial_points[n][m_dimension - 1] = points[n][m_dimension - 1];
		for (int i = 0 ; i < m_dimension - 1 ; ++i) {
			initial_points[n][i] = points[n][i];
			if(initial_points[n][i] < m_region_low[i]) {
				m_region_low[i] = (double)initial_points[n][i];
			}
		}
	}

	// call stream initially
	stream(m_region_low, m_region_up, initial_points, m_total_size, 0, r_point[m_dimension - 1], 0);

	// free the memory for child node points
	for (unsigned int n = 0; n < m_child_points.size() ; ++n) {
		delete[] m_child_points[n];
	}
	m_child_points.clear();

	// free the memory of the initial points
	for (int n = 0; n < m_total_size ; ++n) {
		delete[] initial_points[n];
	}
	delete[] initial_points;

	// free the member variables
	delete[] m_region_low;
	delete[] m_region_up;
	delete[] m_boundaries;
	delete[] m_no_boundaries;
	delete[] m_piles;
	delete[] m_trellis;

	return m_volume;
}

bool hoy::covers(const double cub[], const double reg_low[]) const
{
	for (int i = 0; i < m_dimension - 1; ++i) {
		if (cub[i] > reg_low[i]) {
			return false;
		}
	}
	return true;
}

bool hoy::part_covers(const double cub[], const double reg_up[]) const
{
	for (int i = 0; i < m_dimension - 1; ++i) {
		if (cub[i] >= reg_up[i]) {
			return false;
		}
	}
	return true;
}

int hoy::contains_boundary(const double cub[], const double reg_low[], const int split) const
{
	// condition only checked for split > 0
	if (reg_low[split] >= cub[split]){
		// boundary in m_dimension split not contained in region, thus
		// boundary is no candidate for the splitting line
		return -1;
	}
	else {
		for (int j = 0 ; j < split; ++j) { // check boundaries
			if (reg_low[j] < cub[j]) {
				// boundary contained in region
				return 1;
			}
		}
	}
	// no boundary contained in region
	return 0;
}

double hoy::get_measure(const double reg_low[], const double reg_up[]) const
{
	double vol = 1.0;
	for (int i = 0 ; i < m_dimension - 1 ; ++i) {
		vol *= (reg_up[i] - reg_low[i]);
	}
	return vol;
}

int hoy::is_pile(const double cub[], const double reg_low[]) const
{

	int pile = m_dimension;
	// check all dimensions of the node
	for (int k = 0 ; k < m_dimension - 1 ; ++k) {
		// k-boundary of the node's region contained in the cuboid?
		if (cub[k] > reg_low[k]) {
			if (pile != m_dimension) {
				// second dimension occured that is not completely covered
				// ==> cuboid is no pile
				return -1;
			}
			pile = k;
		}
	}
	// if pile == this.dimension then
	// cuboid completely covers region
	// case is not possible since covering cuboids have been removed before
		
	// region in only one dimenison not completly covered
	// ==> cuboid is a pile
	return pile;
}

double hoy::compute_trellis(const double reg_low[], const double reg_up[], const double trellis[]) const
{

	double total = 1.0; // total area region
	double empty = 1.0; // inner "empty" box for trellis
	for (int i = 0 ; i < m_dimension - 1 ; ++i) {
		// no abs value is required as reg_up[i] > reg_low[i] and trellis[i] >= reg_low[i]
		total *= (reg_up[i] - reg_low[i]);
		empty *= (trellis[i] - reg_low[i]);
	}
	return total - empty;
}

double hoy::get_median(double* bounds, unsigned int n) const
{
	if (n < 3) {
		return bounds[n - 1];
	}
	unsigned int n2 = n/2;
	std::partial_sort(bounds, bounds + n2, bounds + n);
	return bounds[n2];
}

/// Recursive calculation of the hypervolume.
void hoy::stream(double m_region_low[], double m_region_up[], double** points, const unsigned int n_points, int split, double cover, unsigned int rec_level) const
{
	double cover_old = cover;
	unsigned int cover_index = 0;
		
	// Identify first covering cuboid.
	// Compute the D-1 dimensional sweeping hyperplane area
	double plane_area = get_measure(m_region_low, m_region_up);
	while (cover == cover_old && cover_index < n_points) {
		if ( covers(points[cover_index], m_region_low) ) {
			// new cover value
			cover = points[cover_index][m_dimension - 1];
			m_volume += plane_area * (cover_old - cover);
		}
		else ++cover_index;
	}

	/* cover_index shall be the index of the first point in points which
	 * is ignored in the remaining process
	 *
	 * It may occur that that some points in front of cover_index have the same
	 * d-th coordinate as the point at cover_index. This points must be discarded
	 * and therefore the following for-loop checks for this points and reduces
	 * cover_index if necessary.
	 */
	for (int c = cover_index ; c > 0; --c) {
		if (points[c - 1][m_dimension - 1] == cover) {
			--cover_index;
		}
	}
	
	// Abort if points is empty
	if (cover_index == 0) {
		return;
	}
	// Note: in the remainder points is only considered to index cover_index
	
	// Make sure we have a trellis (each box must be an i-pile).
	bool all_piles = true;
	
	for (unsigned int i = 0; i < cover_index; ++i) {
		m_piles[i] = is_pile(points[i], m_region_low);
		if (m_piles[i] == -1) {
			all_piles = false;
			break;
		}
	}
	
	/*
	 * m_trellis[i] contains the values of the minimal i-coordinate of
	 * the i-piles.
	 * If there is no i-pile the default value is the upper bpund of the region.
	 * The 1-dimensional KMP of the i-piles is: reg[1][i] - trellis[i]
	 *
	 */
	if (all_piles) { // Proceed with the sweeping.
		// Initialize trellis with region's upper bound
		for (int c = 0 ; c < m_dimension - 1; ++c) {
			m_trellis[c] = m_region_up[c];
		}
		
		double current = 0.0;
		double next = 0.0;
		unsigned int i = 0;
		do { // while(next != cover)
			current = points[i][m_dimension-1];
			do { // while(next == current)
				if (points[i][m_piles[i]] < m_trellis[m_piles[i]]) {
					m_trellis[m_piles[i]] = points[i][m_piles[i]];
				}
				++i; // index of next point
				if (i < cover_index) {
					next = points[i][m_dimension - 1];
				}
				else {
					next = cover;
				}
			} while(next == current);		
			m_volume += compute_trellis(m_region_low, m_region_up, m_trellis) * (next - current);
		} while(next != cover);
	}
	// If the trellis was NOT detected, partite the node.
	else{
		// Find the splitting boundary.
		double bound = 0.0;
		bool not_found = true;
		unsigned int n_bounds = 0;
		unsigned int n_no_bounds = 0;

		do {
			for (unsigned int i = 0; i < cover_index; ++i) {
				int contained = contains_boundary(points[i], m_region_low, split);
				if (contained == 1) {
					m_boundaries[n_bounds++] = points[i][split];
				} else if (contained == 0) {
					m_no_boundaries[n_no_bounds++] = points[i][split];
				}
			}
			
			//if (boundaries.size() > 0) {
			if (n_bounds > 0) {
				bound = get_median(m_boundaries, n_bounds);
				not_found = false;
			}
			else if (n_no_bounds > m_sqrt_size) {
				bound = get_median(m_no_boundaries, n_no_bounds);
				not_found = false;
			}
			else {
				++split;
			}
		} while (not_found);

		// if new frame for child_points is required
		if (rec_level >= m_child_points.size()) {
			m_child_points.push_back(new double*[m_total_size]);
		}

		unsigned int n_cp = 0;
	
		double d_last = m_region_up[split];
		// Left child
		m_region_up[split] = bound;
		for (unsigned int i = 0; i < cover_index ; ++i) {
			if (part_covers(points[i], m_region_up)) {
				m_child_points[rec_level][n_cp++] = points[i];
			}
		}
		if (n_cp > 0) {
			stream(m_region_low, m_region_up, m_child_points[rec_level], n_cp, split, cover, rec_level + 1);
		}

		// Right child
		n_cp = 0;
		m_region_up[split] = d_last;
		d_last = m_region_low[split];
		m_region_low[split] = bound;
		for (unsigned int i = 0 ; i < cover_index ; ++i) {
			if (part_covers(points[i], m_region_up)) {
				m_child_points[rec_level][n_cp++] = points[i];
			}
		}
		if (n_cp > 0) {
			stream(m_region_low, m_region_up, m_child_points[rec_level], n_cp, split, cover, rec_level + 1);
		}
		m_region_low[split] = d_last;
	}
}

/// Verify before compute method
/**
 * Verifies whether given algorithm suits the requested data.
 *
 * @param[in] points vector of points containing the D-dimensional points for which we compute the hypervolume
 * @param[in] r_point reference point for the vector of points
 *
 * @throws value_error when trying to compute the hypervolume for the non-maximal reference point
 */
void hoy::verify_before_compute(const std::vector<fitness_vector> &points, const fitness_vector &r_point) const
{
	base::assert_minimisation(points, r_point);
}

/// Clone method.
base_ptr hoy::clone() const
{
	return base_ptr(new hoy(*this));
}

/// Algorithm name
std::string hoy::get_name() const
{
	return "HOY algorithm";
}

} } }

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::util::hv_algorithm::hoy)
