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

#ifndef PAGMO_UTIL_HV_ALGORITHM_HV3D_H
#define PAGMO_UTIL_HV_ALGORITHM_HV3D_H

#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <cmath>
#include <set>
#include <algorithm>

#include "base.h"

namespace pagmo { namespace util { namespace hv_algorithm {

/// hv3d hypervolume algorithm class
/**
 * This class contains the implementation of efficient algorithms for the hypervolume computation in 3-dimensions.
 *
 * 'compute' method relies on the efficient algorithm as it was presented by Nicola Beume et al.
 * 'least[greatest]_contributor' methods rely on the HyCon3D algorithm by Emmerich and Fonseca.
 *
 * @see "On the Complexity of Computing the Hypervolume Indicator", Nicola Beume, Carlos M. Fonseca, Manuel Lopez-Ibanez, Luis Paquete, Jan Vahrenhold. IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTATION, VOL. 13, NO. 5, OCTOBER 2009
 * @see "Computing hypervolume contribution in low dimensions: asymptotically optimal algorithm and complexity results", Michael T. M. Emmerich, Carlos M. Fonseca
 *
 * @author Krzysztof Nowak (kn@kiryx.net)
 */
class __PAGMO_VISIBLE hv3d : public base
{
public:
	hv3d(bool initial_sorting = true);
	double compute(std::vector<fitness_vector> &, const fitness_vector &) const;
	std::vector<double> contributions(std::vector<fitness_vector> &, const fitness_vector &) const;

	void verify_before_compute(const std::vector<fitness_vector> &, const fitness_vector &) const;
	base_ptr clone() const;
	std::string get_name() const;

private:
	// flag stating whether the points should be sorted in the first step of the algorithm
	const bool m_initial_sorting;

	struct box3d
	{
		box3d(double _lx, double _ly, double _lz, double _ux, double _uy, double _uz)
			: lx(_lx), ly(_ly), lz(_lz), ux(_ux), uy(_uy), uz(_uz) { }
		double lx;
		double ly;
		double lz;
		double ux;
		double uy;
		double uz;
	};

	struct hycon3d_tree_cmp
	{
		bool operator()(const std::pair<fitness_vector, int> &, const std::pair<fitness_vector, int> &);
	};

	static bool hycon3d_sort_cmp(const std::pair<fitness_vector, unsigned int> &, const std::pair<fitness_vector, unsigned int> &);
	static double box_volume(const box3d &b);

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<bool &>(m_initial_sorting);
	}
};

} } }

BOOST_CLASS_EXPORT_KEY(pagmo::util::hv_algorithm::hv3d)

#endif
