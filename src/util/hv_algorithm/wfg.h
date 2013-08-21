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
class __PAGMO_VISIBLE wfg : public base {
	public:
		wfg(const unsigned int stop_dimension = 2);
		double compute(std::vector<fitness_vector> &, const fitness_vector &);
		void verify_before_compute(const std::vector<fitness_vector> &, const fitness_vector &);
		base_ptr clone() const;
		std::string get_name() const;


	private:

		void limitset(double**, const unsigned int, const unsigned int, const unsigned int);
		double exclusive_hv(double**, const unsigned int, const unsigned int, const unsigned int);
		double compute_hv(double**, const unsigned int, const unsigned int);

		bool cmp_func(double* a, double* b);

		mutable unsigned int m_current_slice;
		const unsigned int m_stop_dimension;

		// Array of point sets for each recursive level.
		double*** m_frames;

		// Maintains the number of points at given recursion level.
		unsigned int* m_frames_size;

		// Keeps track of currently allocated number of frames.
		unsigned int m_n_frames;

		double* m_refpoint;

		unsigned int m_max_points;
		unsigned int m_max_dim;

		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & m_current_slice;
			ar & const_cast<unsigned int &>(m_stop_dimension);
		}
};

} } }

BOOST_CLASS_EXPORT_KEY(pagmo::util::hv_algorithm::wfg);

#endif
