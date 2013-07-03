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


#include "optimal3d.h"

namespace pagmo { namespace util { namespace hv_algorithm {


// Computes hypervolume indicator for given pareto set using the optimal 3D algorithm.
/**
 * This method should be used both as a solution to 3D cases, and as a general termination method for algorithms that reduce n-dimensional problem to 3D.
 *
 * Computational complexity: n*log(n)
 *
 * @param[in] points vector of points containing the d dimensional points for which we compute the hypervolume
 * @param[in] r reference point for the points
 *
 * @return hypervolume.
 */
double optimal3d::compute(const std::vector<fitness_vector> & points, const fitness_vector & r)
{
	// copy the initial set
	std::vector<fitness_vector> points_cpy(points.begin(), points.end());
	sort(points_cpy.begin(), points_cpy.end(), compare_fitness);
	double V = 0.0; // hypervolume
	double A = 0.0; // varying area of the sweeping plane
	std::multiset<fitness_vector, ltcmp> T;

	// sentinel points (r[0], -INF, r[2]) and (-INF, r[1], r[2])
	const double INF = std::numeric_limits<double>::max();
	fitness_vector sA(r.begin(), r.end()); sA[1] = -INF;
	fitness_vector sB(r.begin(), r.end()); sB[0] = -INF;

	T.insert(sA);
	T.insert(sB);
	double z3 = points_cpy[0][2];
	T.insert(points_cpy[0]);
	A = fabs((points_cpy[0][0] - r[0]) * (points_cpy[0][1] - r[1]));

	std::multiset<fitness_vector>::iterator p;
	std::multiset<fitness_vector>::iterator q;
	for(std::vector<fitness_vector>::size_type idx = 1 ; idx < points_cpy.size() ; ++idx) {
		p = T.insert(points_cpy[idx]);
		q = (p);
		++q; //set up q to be a successor of p
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
	V += A * fabs(z3 - r[2]);

	return V;

}

// Verifies whether given method suits the problem
/**
 * This method should be used both as a solution to 3D cases, and as a general termination method for algorithms that reduce n-dimensional problem to 3D.
 *
 * Computational complexity: n*log(n)
 *
 * @param[in] points vector of points containing the d dimensional points for which we compute the hypervolume
 * @param[in] r reference point for the points
 *
 * @throws value_error when trying to compute the hypervolume for the dimension other than 3
 */
void optimal3d::verify_before_compute(const std::vector<fitness_vector> & points, const fitness_vector & r) {
	if (r.size() != 3) {
		pagmo_throw(value_error, "optimal3d method method works only for 3-dimensional cases.");
	}
}

bool compare_fitness(const fitness_vector &a, const fitness_vector &b) {
	return a[2] < b[2];
}

bool ltcmp::operator()(const fitness_vector & a, const fitness_vector & b){
	return a[0] > b[0];
}

} } }
