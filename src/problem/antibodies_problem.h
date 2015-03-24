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

#ifndef PAGMO_PROBLEM_ANTIBODIES_PROBLEM_H
#define PAGMO_PROBLEM_ANTIBODIES_PROBLEM_H

#include <string>
#include <boost/functional/hash.hpp>
#include <boost/serialization/map.hpp>

#include "../serialization.h"
#include "../types.h"
#include "cec2006.h"
#include "base.h"
#include "../algorithm/cstrs_immune_system.h"

namespace pagmo{ namespace problem {

/// Antibodies meta-problem
/**
 * Implements a meta-problem class that simulates the immune system. The objective
 * function computes the distance to the given antigenes either by using the
 * HAMMING or EUCLIDEAN distance. The original problem is needed to set up the
 * boundaries necessary for the HAMMING distance.
 *
 * @see Hajela, P., & Lee, J. (1996). Constrained genetic search via schema adaptation: an immune
 * network solution. Structural optimization, 12(1), 11-15.
 * @see Coello, C. A. C., Cortés, N. C., San, C., & Zacatenco, P. (2001). Use of emulations of
 * the immune system to handle constraints in evolutionary algorithms.
 *
 * @author Jeremie Labroquere (jeremie.labroquere@gmail.com)
 */
class __PAGMO_VISIBLE antibodies_problem : public base
{
//public:
//	// hamming distance, euclidean distance
//	enum method_type {HAMMING = 0, EUCLIDEAN = 1};

public:
	//constructors
	antibodies_problem(const base & = cec2006(4), const algorithm::cstrs_immune_system::distance_method_type = algorithm::cstrs_immune_system::HAMMING);

	//copy constructor
	antibodies_problem(const antibodies_problem &);
	base_ptr clone() const;
	std::string get_name() const;

	void set_antigens(const std::vector<decision_vector> &);

protected:
	std::string human_readable_extra() const;
	void objfun_impl(fitness_vector &, const decision_vector &) const;
	bool compare_fitness_impl(const fitness_vector &, const fitness_vector &) const;

private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & m_original_problem;
		ar & m_pop_antigens;
		ar & const_cast<algorithm::cstrs_immune_system::distance_method_type &>(m_method);
		ar & m_bit_encoding;
		ar & m_max_encoding_integer;
	}

	base_ptr m_original_problem;
	std::vector<decision_vector> m_pop_antigens;

	const algorithm::cstrs_immune_system::distance_method_type m_method;

	// encoding size
	int m_bit_encoding;
	int m_max_encoding_integer;

	// function to compute the distance
	double compute_distance(const decision_vector &x) const;

	// function for the hamming distance
	std::vector<int> double_to_binary(const double &number, const double &lb, const double &ub) const;

};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::antibodies_problem)

#endif // PAGMO_PROBLEM_antibodies_problem_H
