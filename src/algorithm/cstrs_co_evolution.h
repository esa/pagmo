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

#ifndef PAGMO_ALGORITHM_CSTRS_CO_EVOLUTION_H
#define PAGMO_ALGORITHM_CSTRS_CO_EVOLUTION_H

#include <string>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"
#include "jde.h"
#include "sga.h"

namespace pagmo { namespace algorithm {

/// Co-Evolution constraints handling meta-algorithm
/**
 *
 * Co-Evolution constraints handling is a meta-algorithm that allow
 * to solve constrained optimization problems. The key idea of this constraint
 * handling technique is to use two different populations that evolves the one after the
 * other. The first one has an objective function penalized with constraints. The penalty
 * coefficients are encoded in the population 2 which evolves depending on the average
 * population 1 fitness.
 *
 * This meta-algorithm is based on the problems cstrs_co_evolution and cstrs_co_evolution_2.
 *
 * Note: This constraints handling technique can only be used for <b>MINIMIZATION</b> problems.
 *
 * @see Coello Coello, C. A. (2000). Use of a self-adaptive penalty approach for engineering optimization problems.
 * Computers in Industry, 41(2), 113-127.
 * @see Qie He and Ling Wang. 2007. An effective co-evolutionary particle swarm optimization for constrained engineering design problems.
 * Eng. Appl. Artif. Intell. 20, 1 (February 2007), 89-99.
 * DOI=10.1016/j.engappai.2006.03.003 http://dx.doi.org/10.1016/j.engappai.2006.03.003
 *
 * @author Jeremie Labroquere (jeremie.labroquere@gmail.com)
 */
		
class __PAGMO_VISIBLE cstrs_co_evolution: public base
{
public:
	/// Type of co-evolution.
	/**
	* Definition of three types of co-evolution: SIMPLE, SPLIT_NEQ_EQ and SPLIT_CONSTRAINTS.
	* The SIMPLE, is co-evolution defined by COELLO. The SPLIT_NEQ_EQ, splits equalities and
	* inequalities constraints (4 penalty coefficients). The SPLIT_CONSTRAINTS split the
	* number of coefficients upon the number of penalty coefficients (2 * c_dimension).
	*/
	// co-evolution simple, split_neq_eq, split_constraints
	enum method_type {SIMPLE = 0, SPLIT_NEQ_EQ = 1, SPLIT_CONSTRAINTS = 2};

	cstrs_co_evolution(const base & = jde(), const base & = sga(1), int pop_penalties_size = 30, int gen = 1,
					   method_type method = SIMPLE, double pen_lower_bound = 0.,
					   double pen_upper_bound = 100000.,
					   double = 1e-15, double = 1e-15);
	cstrs_co_evolution(const cstrs_co_evolution &);
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
		ar & m_original_algo_penalties;
		ar & const_cast<int &>(m_gen);
		ar & m_method;
		ar & m_pop_penalties_size;
		ar & m_pen_lower_bound;
		ar & m_pen_upper_bound;
		ar & const_cast<double &>(m_ftol);
		ar & const_cast<double &>(m_xtol);
	}
	base_ptr m_original_algo;
	base_ptr m_original_algo_penalties;
	//Number of generations
	const int m_gen;
	// population encoding penalties size
	int m_pop_penalties_size;
	// problem associated to population penalties variables
	method_type m_method;
	double m_pen_lower_bound;
	double m_pen_upper_bound;

	// tolerance
	const double m_ftol;
	const double m_xtol;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::cstrs_co_evolution)

#endif // PAGMO_ALGORITHM_CSTRS_CO_EVOLUTION_H
