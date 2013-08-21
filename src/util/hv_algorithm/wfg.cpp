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


#include "wfg.h"
#include "base.h"
#include <algorithm>
#include <boost/bind.hpp>

#define DEBUG_ON 0
#define SPACER_ON 1
#define SPACER_CHAR "-"
#define SPACER_VAL rec_level
#define DEBUG(x) if(DEBUG_ON) {std::cout << #x << ": " << (x) << std::endl;}
#define DEBUGV(x) if(DEBUG_ON) {std::cout << #x << ":"; for(unsigned int _idx = 0; _idx < (x).size() ; ++ _idx){ std::cout << " " << (x)[_idx];} std::cout << std::endl;}
#define DEBUGA(x,s) if(DEBUG_ON) {std::cout << #x << ":"; for(unsigned int _idx = 0; _idx < (s) ; ++ _idx){ std::cout << " " << (x)[_idx];} std::cout << std::endl;}
#define DEBUGVP(x) if(DEBUG_ON) {std::cout << #x << ":"; for(unsigned int _idx = 0; _idx < (x).size() ; ++ _idx){ std::cout << " (" << (x)[_idx].first << "," << (x)[_idx].second << ")";} std::cout << std::endl;}
#define DEBUGVV(x) if(DEBUG_ON) {std::cout << #x << ":" << std::endl; for(unsigned int _idx = 0; _idx < (x).size() ; ++_idx) {std::cout << " " << #x << "[" << _idx << "]:"; for(unsigned int _idx2 = 0; _idx2 < ((x)[_idx]).size() ; ++ _idx2) {std::cout << " " << (x)[_idx][_idx2]; } std::cout << std::endl; } }
#define PRINT(__x) if(DEBUG_ON) {std::cout << __x << std::endl;}
#define _SPACER() if(SPACER_ON){for(unsigned int _ii= 0 ;  _ii < SPACER_VAL ; ++_ii) {std::cout << SPACER_CHAR;}}
#define SDEBUG(x) if(DEBUG_ON) {_SPACER(); std::cout << #x << ": " << (x) << std::endl;}
#define SDEBUGV(x) if(DEBUG_ON) { _SPACER();std::cout << #x << ":"; for(unsigned int _idx = 0; _idx < (x).size() ; ++ _idx){ std::cout << " " << (x)[_idx];} std::cout << std::endl;}
#define SDEBUGVV(x) if(DEBUG_ON) { _SPACER(); std::cout << #x << ":" << std::endl; for(unsigned int _idx = 0; _idx < (x).size() ; ++_idx) { _SPACER(); std::cout << " " << #x << "[" << _idx << "]:"; for(unsigned int _idx2 = 0; _idx2 < ((x)[_idx]).size() ; ++ _idx2) { std::cout << " " << (x)[_idx][_idx2]; } std::cout << std::endl; } }
#define SPRINT(__x) if(DEBUG_ON) { _SPACER();std::cout << __x << std::endl;}

namespace pagmo { namespace util { namespace hv_algorithm {

bool wfg::cmp_func(double* a, double* b) {
	for(int i = m_current_slice - 1; i >=0 ; --i){
		if (a[i] > b[i]) {
			return true;
		} else if(a[i] < b[i]) {
			return false;
		}
	}
	return false;
}

// Constructor
wfg::wfg(const unsigned int stop_dimension) : m_current_slice(0), m_stop_dimension(stop_dimension) {
	if (stop_dimension < 2 ) {
		pagmo_throw(value_error, "Stop dimension for WFG must be greater than or equal to 2");
	}
}

/// Compute hypervolume 
/**
 * @param[in] points vector of points containing the D-dimensional points for which we compute the hypervolume
 * @param[in] r_point reference point for the points
 *
 * @return hypervolume.
 */
double wfg::compute(std::vector<fitness_vector> &points, const fitness_vector &r_point) {
	m_max_points = points.size();
	m_max_dim = r_point.size();

	m_refpoint = new double[m_max_dim];
	for(unsigned int d_idx = 0 ; d_idx < m_max_dim ; ++d_idx) {
		m_refpoint[d_idx] = r_point[d_idx];
	}

	// Reserve the space beforehand for each level or recursion.
	// WFG with slicing feature will not go recursively deeper than the dimension size.
	m_frames = new double**[m_max_dim];
	m_frames_size = new unsigned int[m_max_dim];

	// Copy the initial set into the frame at index 0.
	double** fr = new double*[m_max_points];
	for(unsigned int p_idx = 0 ; p_idx < m_max_points ; ++p_idx) {
		fr[p_idx] = new double[m_max_dim];
		for(unsigned int d_idx = 0 ; d_idx < m_max_dim ; ++d_idx) {
			fr[p_idx][d_idx] = points[p_idx][d_idx];
		}
	}
	m_frames[0] = fr;
	m_frames_size[0] = m_max_points;
	m_n_frames = 1;

	// variable holding the current "depth" of dimension slicing. We progress by slicing dimensions from the end.
	m_current_slice = m_max_dim;

	double hv = compute_hv(m_frames[0], m_frames_size[0], 1);

	// Free the memory.
	delete[] m_refpoint;

	for(unsigned int fr_idx = 0 ; fr_idx < m_n_frames ; ++fr_idx) {
		for(unsigned int p_idx = 0; p_idx < m_max_points ; ++p_idx) {
			delete[] m_frames[fr_idx][p_idx];
		}
		delete[] m_frames[fr_idx];
	}
	delete[] m_frames;
	delete[] m_frames_size;

	return hv;
}

void wfg::limitset(double** points, const unsigned int n_points, const unsigned int p_idx, const unsigned int rec_level) {
	int no_points = 0;

	double* p = points[p_idx];
	double** frame = m_frames[rec_level];

	for(unsigned int idx = p_idx + 1; idx < n_points; ++idx) {

		for(fitness_vector::size_type f_idx = 0; f_idx < m_current_slice; ++f_idx) {
			frame[no_points][f_idx] = fmax(points[idx][f_idx], p[f_idx]);
		}

		int cmp_results[no_points];
		double* s = frame[no_points];

		bool keep_s = true;

		// check whether any point is dominating the point 's'
		for(int q_idx = 0; q_idx < no_points; ++q_idx) {
			cmp_results[q_idx] = base::dom_cmp(s, frame[q_idx], m_current_slice);
			if (cmp_results[q_idx] == base::DOM_CMP_B_DOMINATES_A) {
				keep_s = false;
				break;
			}
		}
		// if neither is, remove points dominated by 's' (we store that during the first loop)
		if( keep_s ) {
			int prev = 0;
			int next = 0;
			while(next < no_points) {
				if( cmp_results[next] != base::DOM_CMP_A_DOMINATES_B && cmp_results[next] != base::DOM_CMP_A_B_EQUAL) {
					if(prev < next) {
						for(unsigned int d_idx = 0; d_idx < m_current_slice ; ++d_idx) {
							frame[prev][d_idx] = frame[next][d_idx];
						}
					}
					++prev;
				} 
				++next;
			}
			// Append 's' at the end, if prev==next it's not necessary as it's already there.
			if(prev < next) {
				for(unsigned int d_idx = 0; d_idx < m_current_slice ; ++d_idx) {
					frame[prev][d_idx] = s[d_idx];
				}
			}
			no_points = prev + 1;
		}
	}

	m_frames_size[rec_level] = no_points;
}

double wfg::compute_hv(double** points, const unsigned int n_points, const unsigned int rec_level) {

	// If already sliced to dimension at which we use another algorithm.
	if (m_current_slice == m_stop_dimension) {

		if (m_stop_dimension == 2) {
			// Use a very efficient version of native2d.
			return native2d().compute(points, n_points, m_refpoint);
		} else {
			// Let hypervolume object pick the best method otherwise.

			// Copy arrays into vectors
			std::vector<fitness_vector> points_cpy;
			points_cpy.reserve(n_points);
			for(unsigned int i = 0 ; i < n_points ; ++i) {
				points_cpy.push_back(fitness_vector(points[i], points[i] + m_current_slice));
			}
			fitness_vector r_cpy(m_refpoint, m_refpoint + m_current_slice);

			hypervolume hv = hypervolume(points_cpy, false);
			hv.set_copy_points(false);
			return hv.compute(r_cpy);
		}
	} else {
		// Otherwise, sort the points in preparation for the next recursive step
		// Bind the cmp_func method so it can be used as a comparator function, while still being private and having access to member variables
		std::sort(points, points + n_points, boost::bind(&wfg::cmp_func, this, _1, _2));
	}

	double H = 0.0;
	--m_current_slice;
	for(unsigned int i = 0 ; i < n_points ; ++i) {
		H += fabs((points[i][m_current_slice] - m_refpoint[m_current_slice]) * exclusive_hv(points, n_points, i, rec_level));
	}
	++m_current_slice;
	return H;
}

double wfg::exclusive_hv(double** points, const unsigned int n_points, const unsigned int p_idx, const unsigned int rec_level) {
	if(rec_level >= m_n_frames) {
		double** fr = new double*[m_max_points];
		for(unsigned int i = 0 ; i < m_max_points ; ++i) {
			fr[i]=new double[m_current_slice];
		}
		m_frames[m_n_frames] = fr;
		m_frames_size[m_n_frames] = 0;
		++m_n_frames;
	}

	limitset(points, n_points, p_idx, rec_level);

	double H = base::volume_between(points[p_idx], m_refpoint, m_current_slice);

	if (m_frames_size[rec_level] == 1) {
		H -= base::volume_between(m_frames[rec_level][0], m_refpoint, m_current_slice);
	} else if (m_frames_size[rec_level] > 1) {
		H -= compute_hv(m_frames[rec_level], m_frames_size[rec_level], rec_level + 1);
	}

	return H;
}

// verify_before_compute
/**
 * Verifies whether given algorithm suits the requested data.
 *
 * @param[in] points vector of points containing the D-dimensional points for which we compute the hypervolume
 * @param[in] r_point reference point for the vector of points
 *
 * @throws value_error when trying to compute the hypervolume for the non-maximal reference point
 */
void wfg::verify_before_compute(const std::vector<fitness_vector> &points, const fitness_vector &r_point) {
	base::assert_maximal_reference_point(points, r_point);
}

/// Clone method.
base_ptr wfg::clone() const
{
	return base_ptr(new wfg(*this));
}

/// Algorithm name
std::string wfg::get_name() const {
	return "WFG algorithm";
}

} } }

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::util::hv_algorithm::wfg);
