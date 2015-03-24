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

#ifndef PAGMO_PROBLEM_CSTRS_CO_EVOLUTION_H
#define PAGMO_PROBLEM_CSTRS_CO_EVOLUTION_H

#include <string>
#include <boost/functional/hash.hpp>
#include <boost/serialization/map.hpp>

#include "../serialization.h"
#include "../types.h"
#include "cec2006.h"
#include "base.h"
#include "base_meta.h"
#include "../algorithm/cstrs_co_evolution.h"

///Doxygen will ignore whatever is in //! @cond As this problem is only to be used by the equally named algorithm
//! @cond

namespace pagmo{ namespace problem {

/// Constrainted co evolution meta-problem
/**
 * Implements a meta-problem class that wraps some other constrained problems,
 * resulting in self adaptive constraints handling.
 *
 * The key idea of this constraint handling technique is to .
 *
 * @see R., & Wright, J. A. (2003). Self-adaptive fitness formulation for constrained optimization.
 * Evolutionary Computation, IEEE Transactions on, 7(5), 445-455 for the paper introducing the method.
 *
 * @author Jeremie Labroquere (jeremie.labroquere@gmail.com)
 */

class __PAGMO_VISIBLE cstrs_co_evolution : public base_meta
{
public:
	//constructors
	cstrs_co_evolution(const base & = cec2006(4), const algorithm::cstrs_co_evolution::method_type = algorithm::cstrs_co_evolution::SIMPLE);
	cstrs_co_evolution(const base &, const population&, const algorithm::cstrs_co_evolution::method_type = algorithm::cstrs_co_evolution::SIMPLE);

	base_ptr clone() const;
	std::string get_name() const;

	void set_penalty_coeff(const std::vector<double> &);
	int get_penalty_coeff_size();

protected:
	std::string human_readable_extra() const;
	void objfun_impl(fitness_vector &, const decision_vector &) const;

private:
	void compute_penalty(std::vector<double> &, std::vector<int> &, const constraint_vector &) const;

private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & m_penalty_coeff;
		ar & const_cast<algorithm::cstrs_co_evolution::method_type &>(m_method);
		ar & m_map_fitness;
		ar & m_map_constraint;
	}

	std::vector<double> m_penalty_coeff;

	const algorithm::cstrs_co_evolution::method_type m_method;

	//caches for fitness and constraints values to not recompute them
	std::map<std::size_t, fitness_vector> m_map_fitness;
	std::map<std::size_t, constraint_vector> m_map_constraint;

	// no need to serialize the hasher (and impossible)
	boost::hash< std::vector<double> > m_decision_vector_hash;
};


//The cstrs_co_evolution_penalty cannot inherit from base_meta as the problem dimesnion changes 
// and the constructor of base_meta uses set_bounds using the original problem dimension ....
class __PAGMO_VISIBLE cstrs_co_evolution_penalty : public base
{
public:

	//constructors
	cstrs_co_evolution_penalty(const base & = cec2006(4), int dimension = 2, int size = 30);

	//copy constructor
	cstrs_co_evolution_penalty(const cstrs_co_evolution_penalty &);
	base_ptr clone() const;
	std::string get_name() const;

	void update_penalty_coeff(population::size_type &, const decision_vector &, const population  &);

protected:
	std::string human_readable_extra() const;
	void objfun_impl(fitness_vector &, const decision_vector &) const;
	bool compare_fitness_impl(const fitness_vector &, const fitness_vector &) const;

private:
	void compute_penalty(double &, int &, const decision_vector &) const;

private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & m_original_problem;
		ar & m_pop_2_x_vector;
		ar & m_feasible_count_vector;
		ar & m_feasible_fitness_sum_vector;
		ar & m_max_feasible_fitness;
		ar & m_total_sum_viol;
		ar & m_total_num_viol;
	}
	base_ptr m_original_problem;

	std::vector<decision_vector> m_pop_2_x_vector;

	// penalty coefficients
	std::vector<int> m_feasible_count_vector;
	std::vector<double> m_feasible_fitness_sum_vector;
	double m_max_feasible_fitness;

	// penalty for infeasible
	std::vector<double> m_total_sum_viol;
	std::vector<int> m_total_num_viol;
};
}} //namespaces

//! @endcond

BOOST_CLASS_EXPORT_KEY(pagmo::problem::cstrs_co_evolution)
BOOST_CLASS_EXPORT_KEY(pagmo::problem::cstrs_co_evolution_penalty)

#endif // PAGMO_PROBLEM_cstrs_co_evolution_H
