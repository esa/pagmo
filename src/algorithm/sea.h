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

#ifndef PAGMO_ALGORITHM_SEA_H
#define PAGMO_ALGORITHM_SEA_H

#include "../config.h"
#include "../problem/base.h"
#include "../serialization.h"
#include "base.h"

namespace pagmo { namespace algorithm {

/// (N+1)-EA Simple Evolutionary Algorithm
/**
 * Used in many works dealing with
 *
 * At each generation the best individual is selected and offspring is generated performing a global
 * mutation on it (for each dimension of the problem the gene is forced to a random gene inside the
 * contrains).
 *
 * This particular implementation of EA is able to solve integer box-contrained single-objective problems.
 *
 * @see Oliveto, Pietro S., Jun He, and Xin Yao. "Time complexity of evolutionary algorithms for
 * combinatorial optimization: A decade of results." International Journal of Automation and Computing
 * 4.3 (2007): 281-293.
 *
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 *
 */

class __PAGMO_VISIBLE sea: public base
{
public:
	sea(int gen  = 1);
	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;
protected:
	std::string human_readable_extra() const;
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<int &>(m_gen);
	}  
	//Number of generations
	const int m_gen;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::sea)

#endif // PAGMO_ALGORITHM_SEA_H
