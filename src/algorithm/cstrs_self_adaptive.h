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

#ifndef PAGMO_ALGORITHM_CSTRS_SELF_ADAPTIVE_H
#define PAGMO_ALGORITHM_CSTRS_SELF_ADAPTIVE_H

#include <string>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"
#include "jde.h"

namespace pagmo { namespace algorithm {

/// Self-Adaptive Fitness constraints handling meta-algorithm
/**
 *
 * Seld-Adaptive Fitness constraints handling is a meta-algorithm that allows
 * to solve constrained optimization problems. The key idea of this constraint
 * handling technique is to represent the constraint violation by a single
 * infeasibility measure, and to adapt dynamically the penalization of infeasible solutions.
 *
 * This meta-algorithm is based on the problem self-adaptive.
 *
 * Note: This constraints handling technique can only be used for <b>MINIMIZATION</b> problems
 *
 * @see Farmani, R., & Wright, J. A. (2003). Self-adaptive fitness formulation for constrained optimization.
 * Evolutionary Computation, IEEE Transactions on, 7(5), 445-455 for the paper introducing the method.
 *
 * @author Jeremie Labroquere (jeremie.labroquere@gmail.com)
 */
		
class __PAGMO_VISIBLE cstrs_self_adaptive: public base
{
public:
	cstrs_self_adaptive(const base & = jde(), int gen = 1,
						double = 1e-15, double = 1e-15);
    cstrs_self_adaptive(const cstrs_self_adaptive &);
	base_ptr clone() const;

public:
	void evolve(population &) const;
	std::string get_name() const;
	base_ptr get_algorithm() const;
	void set_algorithm(const base &);

protected:
	std::string human_readable_extra() const;

private:

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & m_original_algo;
		ar & const_cast<int &>(m_gen);
		ar & const_cast<double &>(m_ftol);
		ar & const_cast<double &>(m_xtol);
	}
	base_ptr m_original_algo;
	//Number of generations
	const int m_gen;

	// tolerance
	const double m_ftol;
	const double m_xtol;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::cstrs_self_adaptive)

#endif // PAGMO_ALGORITHM_cstrs_self_adaptive_H
