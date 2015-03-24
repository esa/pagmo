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

#ifndef PAGMO_PROBLEM_CSTRS_SELF_ADAPTIVE_H
#define PAGMO_PROBLEM_CSTRS_SELF_ADAPTIVE_H

#include <string>

#include "../serialization.h"
#include "../types.h"
#include "cec2006.h"
#include "base_meta.h"

#include <boost/functional/hash.hpp>
#include <boost/serialization/map.hpp>

///Doxygen will ignore whatever is in //! @cond As this problem is only to be used by the equally named algorithm
//! @cond

namespace pagmo{ namespace problem {

/// Constrainted self adaptive meta-problem
/**
 * Implements a meta-problem class that wraps some other constrained problems,
 * resulting in self adaptive constraints handling.
 *
 * This meta problem is not intended to be used alone.
 * It has been set to be used by the self-adaptive meta-algorithm.
 *
 * The key idea of this constraint handling technique is to represent the
 * constraint violation by a single infeasibility measure, and to adapt
 * dynamically the penalization of infeasible solutions. As the penalization
 * process depends on a given population, this implementation is valid for
 * the provided population only.
 *
 * @see Farmani R., & Wright, J. A. (2003). Self-adaptive fitness formulation for constrained optimization.
 * Evolutionary Computation, IEEE Transactions on, 7(5), 445-455 for the paper introducing the method.
 *
 * @author Jeremie Labroquere (jeremie.labroquere@gmail.com)
 */

class __PAGMO_VISIBLE cstrs_self_adaptive : public base_meta
{
public:
	//constructors
    cstrs_self_adaptive(const base & = cec2006(4));
    cstrs_self_adaptive(const base &, const population &);

	base_ptr clone() const;
	std::string get_name() const;

	//update penalty and fitnesses with the population
	void update_penalty_coeff(const population &pop);

protected:
	std::string human_readable_extra() const;
	void objfun_impl(fitness_vector &, const decision_vector &) const;

private:
	void update_c_scaling(const population &pop);
	double compute_solution_infeasibility(const constraint_vector &c) const;

private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base_meta>(*this);
		ar & const_cast<bool &>(m_apply_penalty_1);
		ar & m_scaling_factor;
		ar & m_c_scaling;
		ar & m_f_hat_down;
		ar & m_f_hat_up;
		ar & m_f_hat_round;
		ar & m_i_hat_down;
		ar & m_i_hat_up;
		ar & m_i_hat_round;
		ar & m_map_fitness;
		ar & m_map_constraint;
	}

	bool m_apply_penalty_1;
	double m_scaling_factor;

	constraint_vector m_c_scaling;

	fitness_vector m_f_hat_down;
	fitness_vector m_f_hat_up;
	fitness_vector m_f_hat_round;

	double m_i_hat_down;
	double m_i_hat_up;
	double m_i_hat_round;

	std::map<std::size_t, fitness_vector> m_map_fitness;
	std::map<std::size_t, constraint_vector> m_map_constraint;

	// no need to serialize the hasher (and impossible)
	boost::hash< std::vector<double> > m_decision_vector_hash;
};
}} //namespaces

//! @endcond

BOOST_CLASS_EXPORT_KEY(pagmo::problem::cstrs_self_adaptive)

#endif // PAGMO_PROBLEM_cstrs_self_adaptive_H
