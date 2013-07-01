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

#include "hypervolume.h"

namespace pagmo { namespace util {

// hypervolume constructor
/**
 * Constructs a hypervolume object, where points are elicited from the referenced population object.
 *
 * @param[in] pop reference to population object from which pareto front is computed
 */
hypervolume::hypervolume(const population & pop) : m_pop(&pop) {
	this->m_f_dim = this->m_pop->problem().get_f_dimension();
	std::vector<std::vector<population::size_type> > pareto_fronts = this->m_pop->compute_pareto_fronts();
	this->m_points.resize(pareto_fronts[0].size());
	std::deque<std::pair<fitness_vector, fitness_vector::size_type> > point_set;
	for (population::size_type idx = 0 ; idx < pareto_fronts[0].size() ; ++idx) {
		this->m_points[idx] = fitness_vector(this->m_pop->get_individual(pareto_fronts[0][idx]).cur_f);
	}

	verify_after_construct();
}

// hypervolume constructor
/**
 * Constructs a hypervolume object from a provided set of points.
 *
 * @params[in] 
 */
hypervolume::hypervolume(const std::vector<fitness_vector> & points) : m_pop(NULL), m_points(points) {
	this->m_f_dim = m_points[0].size();

	verify_after_construct();
}

// verify after construct
/**
 * Verifies whether basic requirements are met for the initial set of points
 *
 * @throws value_error if point size is empty or when the dimensions among the points differ
 */
void hypervolume::verify_after_construct() {
	if ( this->m_points.size() == 0 ) {
		pagmo_throw(value_error, "Point set cannot be empty.");
	}
	fitness_vector::size_type reference_size = this->m_points[0].size();
	for (std::vector<fitness_vector>::size_type idx = 1 ; idx < this->m_points.size() ; ++idx) {
		if ( this->m_points[idx].size() != reference_size ) {
			pagmo_throw(value_error, "All point set dimensions must be equal.");
		}
	}
}

// verify before compute
/**
 * Verifies whether reference point and the hypervolume method meet certain criteria
 *
 * @params[in] reference_point fitness vector describing the reference point
 *
 * @throws value_error if reference_point's and point set dimension do not agree
 */
void hypervolume::verify_before_compute(const fitness_vector & reference_point, hv_algorithm::base &hv_algorithm) {
	if ( this->m_points[0].size() != reference_point.size() ) {
		pagmo_throw(value_error, "Point set dimensions and reference_point dimension must be equal.");
	}
	hv_algorithm.verify_before_compute(this->m_points, reference_point);
}

// compute hypervolume
/*
 * Computes hypervolume provided a reference_point and an algorithm object
 *
 * @params[in] reference_point fitness vector describing the reference point
 * @params[in] hv_algorithm algorithm object used for computing the hypervolume
 *
 * @return hypervolume
 */
double hypervolume::compute(const fitness_vector & reference_point, hv_algorithm::base & hv_algorithm) {
	this->verify_before_compute(reference_point, hv_algorithm);
	return hv_algorithm.compute(this->m_points, reference_point);
}

}}
