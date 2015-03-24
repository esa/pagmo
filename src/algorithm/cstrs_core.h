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

#ifndef PAGMO_ALGORITHM_CSTRS_CORE_H
#define PAGMO_ALGORITHM_CSTRS_CORE_H

#include <string>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"
#include "jde.h"
#include "sga.h"

namespace pagmo { namespace algorithm {

/// CORE constraints handling meta-algorithm
/**
 *
 * CORE is a meta-algorithm that allow to solve constrained optimization problems.
 * The key idea of this constraint handling technique is to repair infeasible individuals
 * of the population at fixed intervals.
 *
 * The original algorithm is using a random evolutionary search algorithm. This implementation
 * is extended to be used with any algorithm.
 *
 * This meta-algorithm is based on the population population::repair() method.
 *
 * @see Belur, S.V. CORE: Constrained optimization by random evolution, Late Breaking Papers
 * at the Genetic Programming Conference, Stanford University, 280-286, 1997.
 *
 * @author Jeremie Labroquere (jeremie.labroquere@gmail.com)
 */

class __PAGMO_VISIBLE cstrs_core: public base
{
public:
    cstrs_core(const base & = jde(1), const base & = jde(1),
               int = 1,
			   int = 10,
			   double = 1.,
			   double = 1e-15, double = 1e-15);
	cstrs_core(const cstrs_core &);
	base_ptr clone() const;

public:
	void evolve(population &) const;
	std::string get_name() const;
	base_ptr get_algorithm() const;
	void set_algorithm(const base &);
    base_ptr get_repair_algorithm() const;
    void set_repair_algorithm(const base &);

protected:
	std::string human_readable_extra() const;

private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & m_original_algo;
        ar & m_repair_algo;
        ar & const_cast<int &>(m_gen);
		ar & const_cast<int &>(m_repair_frequency);
		ar & const_cast<double &>(m_repair_ratio);
		ar & const_cast<double &>(m_ftol);
		ar & const_cast<double &>(m_xtol);
	}
	base_ptr m_original_algo;
    base_ptr m_repair_algo;
	//Number of generations
	const int m_gen;
    // repair constants
	const int m_repair_frequency;
	const double m_repair_ratio;

	// tolerance
	const double m_ftol;
	const double m_xtol;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::cstrs_core)

#endif // PAGMO_ALGORITHM_CSTRS_CORE_H
