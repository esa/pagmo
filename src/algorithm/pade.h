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

#ifndef PAGMO_ALGORITHM_PADE_H
#define PAGMO_ALGORITHM_PADE_H

#include "../config.h"
#include "../serialization.h"
#include "base.h"
#include "jde.h"



namespace pagmo { namespace algorithm {

/// Parallel Decomposition (PaDe)
/**
 *
 * This class implement a multi-objective optimization algorithm based on parallel decomposition.
 * For each element of the population a different single objective problem is generated having as fitness function a random
 * convex combination of the different objectives. Those single-objective problems are thus solved in parallel.
 * At the end of the evolution the population is set as the best individual for each single-objective problem.
 *
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 **/

class __PAGMO_VISIBLE pade: public base
{
public:
    pade(int gen=100, const pagmo::algorithm::base & = pagmo::algorithm::jde());
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
        ar & m_solver;
	}
	//Number of generations
	const int m_gen;
    base_ptr m_solver;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::pade);

#endif // PAGMO_ALGORITHM_PADE_H
