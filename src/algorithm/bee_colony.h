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

#ifndef PAGMO_ALGORITHM_BEE_COLONY_H
#define PAGMO_ALGORITHM_BEE_COLONY_H

#include <string>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"


namespace pagmo { namespace algorithm {

/// The Artificial Bee Colony Solver (ABC)
/**
 * \image html bee.jpg "Bee"
 * \image latex bee.jpg  "Bee" width=3cm
 * The Artificial Bee Colony (ABC) algorithm is a meta-heuristic algorithm inspired by the
 * behaviour of bees.
 *
 * At each call of the evolve method a number of function evaluations equal
 * to 2 * gen * pop.size() is performed.
 *
 * NOTE: when called on mixed-integer problems ABC treats the integer part as fixed and optimizes
 * the continuous part.
 *
 * @see http://mf.erciyes.edu.tr/abc/pub/ABC.C
 * @see http://www.scholarpedia.org/article/Artificial_bee_colony_algorithm
 *
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 */

class __PAGMO_VISIBLE bee_colony: public base
{
public:
	bee_colony(int gen = 1, int limit = 20);
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
		ar & const_cast<int &>(m_iter);
		ar & const_cast<int &>(m_limit);   
	}  
	// Number of generations.
	const int m_iter;
	const int m_limit;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::bee_colony)

#endif // BEECOLONY_H
