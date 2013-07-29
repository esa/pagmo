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


#include "beume3d.h"

namespace pagmo { namespace util { namespace hv_algorithm {


/// Compute hypervolume 
/**
 * This method should be used both as a solution to 3D cases, and as a general termination method for algorithms that reduce n-dimensional problem to 3D.
 *
 * Computational complexity: n*log(n)
 *
 * @param[in] points vector of points containing the d dimensional points for which we compute the hypervolume
 * @param[in] r_point reference point for the points
 *
 * @return hypervolume.
 */
double beume3d::compute(const std::vector<fitness_vector> &points, const fitness_vector &r_point)
{
	// copy the initial set
	std::vector<fitness_vector> points_cpy(points.begin(), points.end());
	sort(points_cpy.begin(), points_cpy.end(), fitness_vector_cmp(2,'<'));
	double V = 0.0; // hypervolume
	double A = 0.0; // varying area of the sweeping plane
	std::multiset<fitness_vector, fitness_vector_cmp> T(fitness_vector_cmp(0, '>'));

	// sentinel points (r_point[0], -INF, r_point[2]) and (-INF, r_point[1], r_point[2])
	const double INF = std::numeric_limits<double>::max();
	fitness_vector sA(r_point.begin(), r_point.end()); sA[1] = -INF;
	fitness_vector sB(r_point.begin(), r_point.end()); sB[0] = -INF;

	T.insert(sA);
	T.insert(sB);
	double z3 = points_cpy[0][2];
	T.insert(points_cpy[0]);
	A = fabs((points_cpy[0][0] - r_point[0]) * (points_cpy[0][1] - r_point[1]));

	std::multiset<fitness_vector>::iterator p;
	std::multiset<fitness_vector>::iterator q;
	for(std::vector<fitness_vector>::size_type idx = 1 ; idx < points_cpy.size() ; ++idx) {
		p = T.insert(points_cpy[idx]);
		q = (p);
		++q; //setup q to be a successor of p
		if ( (*q)[1] <= (*p)[1] ) { // current point is dominated
			T.erase(p); // disregard the point from further calculation
		} else {
			V += A * fabs(z3 - (*p)[2]);
			z3 = (*p)[2];
			std::multiset<fitness_vector>::reverse_iterator rev_it(q);
			++rev_it;

			std::multiset<fitness_vector>::reverse_iterator erase_begin (rev_it);
			std::multiset<fitness_vector>::reverse_iterator rev_it_pred;
			while((*rev_it)[1] >= (*p)[1] ) {
				rev_it_pred = rev_it;
				++rev_it_pred;
				A -= fabs(((*rev_it)[0] - (*rev_it_pred)[0])*((*rev_it)[1] - (*q)[1]));
				++rev_it;
			}
			A += fabs(((*p)[0] - (*(rev_it))[0])*((*p)[1] - (*q)[1]));
			T.erase(rev_it.base(),erase_begin.base());
		}
	}
	V += A * fabs(z3 - r_point[2]);

	return V;
}


// verify_before_compute
/**
 * Verifies whether given algorithm suits the requested data.
 *
 * @param[in] points vector of points containing the d dimensional points for which we compute the hypervolume
 * @param[in] r_point reference point for the vector of points
 *
 * @throws value_error when trying to compute the hypervolume for the dimension other than 3 or non-maximal reference point
 */
void beume3d::verify_before_compute(const std::vector<fitness_vector> &points, const fitness_vector &r_point) {
	if (r_point.size() != 3) {
		pagmo_throw(value_error, "Algorithm Beume3D works only for 3-dimensional cases");
	}

	base::assert_maximal_reference_point(points, r_point);
}

/// Clone method.
base_ptr beume3d::clone() const
{
	return base_ptr(new beume3d(*this));
}

/// Algorithm name
std::string beume3d::get_name() const {
	return "Beume3D algorithm";
}


} } }

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::util::hv_algorithm::beume3d);
