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

#ifndef PAGMO_UTIL_HV_ALGORITHM_WFG_H
#define PAGMO_UTIL_HV_ALGORITHM_WFG_H

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iterator>

#include "base.h"
#include "../hypervolume.h"

namespace pagmo { namespace util { namespace hv_algorithm {

/// WFG hypervolume algorithm
/**
 * This is the class containing the implementation of the WFG algorithm for the computation of hypervolume indicator.
 *
 * @see "While, Lyndon, Lucas Bradstreet, and Luigi Barone. "A fast way of calculating exact hypervolumes." Evolutionary Computation, IEEE Transactions on 16.1 (2012): 86-95."
 * @see "Lyndon While and Lucas Bradstreet. Applying the WFG Algorithm To Calculate Incremental Hypervolumes. 2012 IEEE Congress on Evolutionary Computation. CEC 2012, pages 489-496. IEEE, June 2012."
 *
 * @author Krzysztof Nowak (kn@kiryx.net)
 */
class __PAGMO_VISIBLE wfg : public base
{
public:
	wfg(const unsigned int stop_dimension = 2);
	double compute(std::vector<fitness_vector> &, const fitness_vector &) const;
	std::vector<double> contributions(std::vector<fitness_vector> &, const fitness_vector &) const;
	void verify_before_compute(const std::vector<fitness_vector> &, const fitness_vector &) const;
	base_ptr clone() const;
	std::string get_name() const;

private:
	void limitset(const unsigned int, const unsigned int, const unsigned int) const;
	double exclusive_hv(const unsigned int, const unsigned int) const;
	double compute_hv(const unsigned int) const;

	bool cmp_points(double* a, double* b) const;

	void allocate_wfg_members(std::vector<fitness_vector> &, const fitness_vector &) const;
	void free_wfg_members() const;

	/**
	 * 'compute' and 'extreme_contributor' method variables section.
	 *
	 * Variables below (especially the pointers m_frames, m_frames_size and m_refpoint) are initialized
	 * at the beginning of the 'compute' and 'extreme_contributor' methods, and freed afterwards.
	 * The state of the variables is irrelevant outside the scope of the these methods.
	 */

	// Current slice depth
	mutable unsigned int m_current_slice;

	// Array of point sets for each recursive level.
	mutable double*** m_frames;

	// Maintains the number of points at given recursion level.
	mutable unsigned int* m_frames_size;

	// Keeps track of currently allocated number of frames.
	mutable unsigned int m_n_frames;

	// Copy of the reference point
	mutable double* m_refpoint;

	// Size of the original front
	mutable unsigned int m_max_points;

	// Size of the dimension
	mutable unsigned int m_max_dim;
	/**
	 * End of 'compute' method variables section.
	 */

	// Dimension at which WFG stops the slicing
	const unsigned int m_stop_dimension;

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<unsigned int &>(m_stop_dimension);
	}
};

} } }

BOOST_CLASS_EXPORT_KEY(pagmo::util::hv_algorithm::wfg)

#endif
