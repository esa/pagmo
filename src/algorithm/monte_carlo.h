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

#ifndef PAGMO_ALGORITHM_MONTE_CARLO_H
#define PAGMO_ALGORITHM_MONTE_CARLO_H

#include <cstddef>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"

namespace pagmo { namespace algorithm {

/// Monte Carlo algorithm.
/**
 * This algorithm will simply evaluate random values within the problem bounds
 * for the number of times specified and return the best.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE monte_carlo: public base
{
	public:
		monte_carlo(int = 1);
		base_ptr clone() const;
		void evolve(population &) const;
		std::string get_name() const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & const_cast<std::size_t &>(m_max_eval);
		}  
		std::string human_readable_extra() const;
		const std::size_t m_max_eval;
};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::monte_carlo)

#endif
