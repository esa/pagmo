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

#include <cmath>
#include <algorithm>

#include "../exceptions.h"
#include "../types.h"
#include "../population.h"
#include "cstrs_self_adaptive.h"

///Doxygen will ignore whatever is in //! @cond As this problem is only to be used by the equally named algorithm
//! @cond

namespace pagmo { namespace problem {

/**
 * Constructor using initial constrained problem
 *
 * Note: This problem is not inteded to be used by itself. Instead use the
 * self-adaptive algorithm if you want to solve constrained problems.
 *
 * @param[in] problem base::problem to be modified to use a self-adaptive
 * as constraints handling technique.
 * @param[in] pop population to be used to set up the penalty coefficients.
 *
 */
cstrs_self_adaptive::cstrs_self_adaptive(const base &problem, const population &pop):
	base_meta(
		 problem,
		 problem.get_dimension(),
		 problem.get_i_dimension(),
		 problem.get_f_dimension(),
		 0,
		 0,
		 std::vector<double>()),
	m_apply_penalty_1(false),
	m_scaling_factor(0.0),
	m_c_scaling(problem.get_c_dimension(),0.0),
	m_f_hat_down(problem.get_f_dimension(),0.0),
	m_f_hat_up(problem.get_f_dimension(),0.0),
	m_f_hat_round(problem.get_f_dimension(),0.0),
	m_i_hat_down(0.0),
	m_i_hat_up(0.0),
	m_i_hat_round(0.0),
	m_map_fitness(),
	m_map_constraint(),
	m_decision_vector_hash()
{
	if(m_original_problem->get_c_dimension() <= 0){
		pagmo_throw(value_error,"The original problem has no constraints.");
	}

	// check that the dimension of the problem is 1
	if (m_original_problem->get_f_dimension() != 1) {
		pagmo_throw(value_error,"The original fitness dimension of the problem must be one, multi objective problems can't be handled with self adaptive meta problem.");
	}

	if(problem != pop.problem()) {
		pagmo_throw(value_error,"The problem linked to the population is not the same as the problem given in argument.");
	}

	update_penalty_coeff(pop);
}

cstrs_self_adaptive::cstrs_self_adaptive(const base &problem):
	base_meta(
		 problem,
		 problem.get_dimension(),
		 problem.get_i_dimension(),
		 problem.get_f_dimension(),
		 0,
		 0,
		 std::vector<double>()),
	m_apply_penalty_1(false),
	m_scaling_factor(0.0),
	m_c_scaling(problem.get_c_dimension(),0.0),
	m_f_hat_down(problem.get_f_dimension(),0.0),
	m_f_hat_up(problem.get_f_dimension(),0.0),
	m_f_hat_round(problem.get_f_dimension(),0.0),
	m_i_hat_down(0.0),
	m_i_hat_up(0.0),
	m_i_hat_round(0.0),
	m_map_fitness(),
	m_map_constraint(),
	m_decision_vector_hash()
{
	population pop(*m_original_problem,0);

	if(m_original_problem->get_c_dimension() <= 0){
		pagmo_throw(value_error,"The original problem has no constraints.");
	}

	// check that the dimension of the problem is 1
	if (m_original_problem->get_f_dimension() != 1) {
		pagmo_throw(value_error,"The original fitness dimension of the problem must be one, multi objective problems can't be handled with self adaptive meta problem.");
	}
	update_penalty_coeff(pop);
}

/// Clone method.
base_ptr cstrs_self_adaptive::clone() const
{
	return base_ptr(new cstrs_self_adaptive(*this));
}

/// Implementation of the objective function.
/// (Wraps over the original implementation)
/**
 *  Returns the penalized objective value.
 */
void cstrs_self_adaptive::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	std::map<std::size_t, fitness_vector>::const_iterator it_f;
	std::map<std::size_t, constraint_vector>::const_iterator it_c;

	double solution_infeasibility;

	it_f = m_map_fitness.find(m_decision_vector_hash(x));
	if(it_f != m_map_fitness.end()) {
		// we suppose that the constraints is also available
		it_c = m_map_constraint.find(m_decision_vector_hash(x));

		f = it_f->second;

		solution_infeasibility = compute_solution_infeasibility(it_c->second);
	} else {
		// we compute the function f
		m_original_problem->objfun(f, x);

		// we compute the constraints
		constraint_vector c(m_original_problem->get_c_dimension(), 0.);
		m_original_problem->compute_constraints(c,x);

		solution_infeasibility = compute_solution_infeasibility(c);
	}

	if(solution_infeasibility > 0.) {
		// apply penalty 1
		if(m_apply_penalty_1) {
			double inf_tilde = 0.;
			inf_tilde = (solution_infeasibility - m_i_hat_down) /
					(m_i_hat_up - m_i_hat_down);

			f[0] += inf_tilde * (m_f_hat_down[0] - m_f_hat_up[0]);
		}

		// apply penalty 2
		f[0] += m_scaling_factor * std::fabs(f[0]) * ( (std::exp(2. * solution_infeasibility) - 1.) / (std::exp(2.) - 1.) );
	}
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the type of constraint handling
 */
std::string cstrs_self_adaptive::human_readable_extra() const
{
	std::ostringstream oss;
	oss << m_original_problem->human_readable_extra() << std::endl;
	oss << "\n\tConstraints handled with self-adaptive method ";
	oss << std::endl;
	return oss.str();
}

std::string cstrs_self_adaptive::get_name() const
{
	return m_original_problem->get_name() + " [cstrs_self_adaptive]";
}

/// Updates the penalty coefficients needed to compute the penalized
/// objective function.
/**
 * Updates the penalty coefficients with the given population.
 * @param[in] population pop.
 */
void cstrs_self_adaptive::update_penalty_coeff(const population &pop)
{
	if(*m_original_problem != pop.problem()) {
		pagmo_throw(value_error,"The problem linked to the population is not the same as the problem given in argument.");
	}

	// Let's store some useful variables.
	const population::size_type pop_size = pop.size();

	// Get out if there is nothing to do.
	if (pop_size == 0) {
		return;
	}

	m_map_fitness.clear();
	m_map_constraint.clear();
	// store f and c in maps depending on the the x hash
	for(population::size_type i=0; i<pop_size; i++) {
		const population::individual_type &current_individual = pop.get_individual(i);
		// m_map_fitness.insert(std::pair<std::size_t, fitness_vector>(m_decision_vector_hash(current_individual.cur_x),current_individual.cur_f));
		m_map_fitness[m_decision_vector_hash(current_individual.cur_x)]=current_individual.cur_f;
		m_map_constraint[m_decision_vector_hash(current_individual.cur_x)]=current_individual.cur_c;
	}

	std::vector<population::size_type> feasible_idx(0);
	std::vector<population::size_type> infeasible_idx(0);

	// store indexes of feasible and non feasible individuals
	for(population::size_type i=0; i<pop_size; i++) {
		const population::individual_type &current_individual = pop.get_individual(i);

		if(m_original_problem->feasibility_c(current_individual.cur_c)) {
			feasible_idx.push_back(i);
		} else {
			infeasible_idx.push_back(i);
		}
	}

	// if the population is only feasible, then nothing is done
	if(infeasible_idx.size() == 0) {
		update_c_scaling(pop);		
		m_apply_penalty_1 = false;
		m_scaling_factor = 0.;
		return;
	}
	m_apply_penalty_1 = false;
	m_scaling_factor = 0.;

	// updates the c_scaling, needed for solution infeasibility computation
	update_c_scaling(pop);

	// evaluate solutions infeasibility
	//compute_pop_solution_infeasibility(solution_infeasibility, pop);

	std::vector<double> solution_infeasibility(pop_size);
	std::fill(solution_infeasibility.begin(),solution_infeasibility.end(),0.);

	// evaluate solutions infeasibility
	solution_infeasibility.resize(pop_size);
	std::fill(solution_infeasibility.begin(),solution_infeasibility.end(),0.);

	for(population::size_type i=0; i<pop_size; i++) {
		const population::individual_type &current_individual = pop.get_individual(i);

		// compute the infeasibility of the constraint
		solution_infeasibility[i] = compute_solution_infeasibility(current_individual.cur_c);
	}

	// search position of x_hat_down, x_hat_up and x_hat_round
	population::size_type hat_down_idx = -1;
	population::size_type hat_up_idx = -1;
	population::size_type hat_round_idx = -1;

	// first case, the population contains at least one feasible solution
	if(feasible_idx.size() > 0) {
		// initialize hat_down_idx
		hat_down_idx = feasible_idx.at(0);

		// x_hat_down = feasible individual with lowest objective value in p
		for(population::size_type i=0; i<feasible_idx.size(); i++) {
			const population::size_type current_idx = feasible_idx.at(i);
			const population::individual_type &current_individual = pop.get_individual(current_idx);

			if(m_original_problem->compare_fitness(current_individual.cur_f, pop.get_individual(hat_down_idx).cur_f)) {
				hat_down_idx = current_idx;
			}
		}

		// hat down is now available
		fitness_vector f_hat_down = pop.get_individual(hat_down_idx).cur_f;

		// x_hat_up value depends if the population contains infeasible individual with objective
		// function better than f_hat_down
		bool pop_contains_infeasible_f_better_x_hat_down = false;
		for(population::size_type i=0; i<infeasible_idx.size(); i++) {
			const population::size_type current_idx = infeasible_idx.at(i);
			const population::individual_type &current_individual = pop.get_individual(current_idx);

			if(m_original_problem->compare_fitness(current_individual.cur_f, f_hat_down)) {
				pop_contains_infeasible_f_better_x_hat_down = true;

				// initialize hat_up_idx
				hat_up_idx = current_idx;

				break;
			}
		}

		if(pop_contains_infeasible_f_better_x_hat_down) {
			// hat_up_idx is already initizalized

			// gets the individual with maximum infeasibility and objfun lower than f_hat_down
			for(population::size_type i=0; i<infeasible_idx.size(); i++) {
				const population::size_type current_idx = infeasible_idx.at(i);
				const population::individual_type &current_individual = pop.get_individual(current_idx);

				if(m_original_problem->compare_fitness(current_individual.cur_f, f_hat_down) &&
					(solution_infeasibility.at(current_idx) >= solution_infeasibility.at(hat_up_idx)) ) {

					if(solution_infeasibility.at(current_idx) == solution_infeasibility.at(hat_up_idx)) {
						if(m_original_problem->compare_fitness(current_individual.cur_f, pop.get_individual(hat_up_idx).cur_f)) {
							hat_up_idx = current_idx;
						}
					} else {
						hat_up_idx = current_idx;
					}
				}
			}

			// apply penalty 1
			m_apply_penalty_1 = true;

		} else {
			// all the infeasible soutions have an objective function value greater than f_hat_down
			// the worst is the one that has the maximum infeasibility
			// initialize hat_up_idx
			hat_up_idx = infeasible_idx.at(0);

			for(population::size_type i=0; i<infeasible_idx.size(); i++) {
				const population::size_type current_idx = infeasible_idx.at(i);
				const population::individual_type &current_individual = pop.get_individual(current_idx);

				if(solution_infeasibility.at(current_idx) >= solution_infeasibility.at(hat_up_idx)) {
					if(solution_infeasibility.at(current_idx) == solution_infeasibility.at(hat_up_idx)) {
						if(m_original_problem->compare_fitness(pop.get_individual(hat_up_idx).cur_f, current_individual.cur_f)) {
							hat_up_idx = current_idx;
						}
					} else {
						hat_up_idx = current_idx;
					}
				}
			}

			// do not apply penalty 1
			m_apply_penalty_1 = false;
		}

	} else { // case where there is no feasible solution in the population
		// best is the individual with the lowest infeasibility
		hat_down_idx = 0;
		hat_up_idx = 0;

		for(population::size_type i=0; i<pop_size; i++) {
			const population::individual_type &current_individual = pop.get_individual(i);

			if(solution_infeasibility.at(i) <= solution_infeasibility.at(hat_down_idx)) {
				if(solution_infeasibility.at(i) == solution_infeasibility.at(hat_down_idx)) {
					if(m_original_problem->compare_fitness(current_individual.cur_f, pop.get_individual(hat_down_idx).cur_f)) {
						hat_down_idx = i;
					}
				} else {
					hat_down_idx = i;
				}
			}
		}

		// worst individual
		for(population::size_type i=0; i<pop_size; i++) {
			const population::individual_type &current_individual = pop.get_individual(i);

			if(solution_infeasibility.at(i) >= solution_infeasibility.at(hat_up_idx)) {
				if(solution_infeasibility.at(i) == solution_infeasibility.at(hat_up_idx)) {
					if(m_original_problem->compare_fitness(pop.get_individual(hat_up_idx).cur_f, current_individual.cur_f)) {
						hat_up_idx = i;
					}
				} else {
					hat_up_idx = i;
				}
			}
		}

		// apply penalty 1 to the population
		m_apply_penalty_1 = true;
	}

	// stores the hat round idx, i.e. the solution with highest objective
	// function value in the population
	hat_round_idx = 0;
	for(population::size_type i=0; i<pop_size; i++) {
		const population::individual_type &current_individual = pop.get_individual(i);

		if(m_original_problem->compare_fitness(pop.get_individual(hat_round_idx).cur_f, current_individual.cur_f)) {
			hat_round_idx = i;
		}
	}

	// get the objective function values of the three individuals
	m_f_hat_round = pop.get_individual(hat_round_idx).cur_f;
	m_f_hat_down =  pop.get_individual(hat_down_idx).cur_f;
	m_f_hat_up = pop.get_individual(hat_up_idx).cur_f;

	// get the solution infeasibility values of the three individuals
	m_i_hat_round = solution_infeasibility.at(hat_round_idx);
	m_i_hat_down = solution_infeasibility.at(hat_down_idx);
	m_i_hat_up = solution_infeasibility.at(hat_up_idx);

	// computes the scaling factor
	m_scaling_factor = 0.;
	// evaluates scaling factor
	if(m_original_problem->compare_fitness(m_f_hat_down, m_f_hat_up)) {
		m_scaling_factor = (m_f_hat_round[0] - m_f_hat_up[0]) / m_f_hat_up[0];
	} else {
		m_scaling_factor = (m_f_hat_round[0] - m_f_hat_down[0]) / m_f_hat_down[0];
	}
	if(m_f_hat_up[0] == m_f_hat_round[0]) {
		m_scaling_factor = 0.;
	}

}

/// Updates the constraints scaling.
/**
 * Updates the constraints scaling vector with the given population.
 * @param[in] population pop.
 */
void cstrs_self_adaptive::update_c_scaling(const population &pop)
{
	if(*m_original_problem != pop.problem()) {
		pagmo_throw(value_error,"The problem linked to the population is not the same as the problem given in argument.");
	}

	// Let's store some useful variables.
	const population::size_type pop_size = pop.size();

	// get the constraints dimension
	//constraint_vector c(m_original_problem->get_c_dimension(), 0.);
	problem::base::c_size_type prob_c_dimension = m_original_problem->get_c_dimension();
	problem::base::c_size_type number_of_eq_constraints =
			m_original_problem->get_c_dimension() -
			m_original_problem->get_ic_dimension();

	const std::vector<double> &c_tol = m_original_problem->get_c_tol();

	m_c_scaling.resize(m_original_problem->get_c_dimension());
	std::fill(m_c_scaling.begin(),m_c_scaling.end(),0.);

	// evaluates the scaling factor
	for(population::size_type i=0; i<pop_size; i++) {
		// updates the current constraint vector
		const population::individual_type &current_individual = pop.get_individual(i);

		const constraint_vector &c = current_individual.cur_c;

		// computes scaling with the right definition of the constraints (can be in base problem? currently used
		// by con2mo as well)
		for(problem::base::c_size_type j=0; j<number_of_eq_constraints; j++) {
			m_c_scaling[j] = std::max(m_c_scaling[j], std::max(0., (std::abs(c.at(j)) - c_tol.at(j))) );
		}
		for(problem::base::c_size_type j=number_of_eq_constraints; j<prob_c_dimension; j++) {
			m_c_scaling[j] = std::max(m_c_scaling[j], std::max(0., c.at(j) - c_tol.at(j)) );
		}
	}
}

/// Computes the solution infeasibility measure for the given constraint,
/// need the constraints scaling to be updated before calling this method.
/**
 * Updates the solution infeasibility vector with the population given.
 * @param[in] constraint_vector c.
 * @param[out] solution infeasibility.
 */
double cstrs_self_adaptive::compute_solution_infeasibility(const constraint_vector &c) const
{
	// get the constraints dimension
	problem::base::c_size_type prob_c_dimension = m_original_problem->get_c_dimension();
	problem::base::c_size_type number_of_eq_constraints =
			m_original_problem->get_c_dimension() -
			m_original_problem->get_ic_dimension();

	double solution_infeasibility = 0.;

	const std::vector<double> &c_tol = m_original_problem->get_c_tol();

	// computes solution infeasibility with the right definition of the constraints (can be in base problem? currently used
	// by con2mo as well)
	for(problem::base::c_size_type j=0; j<number_of_eq_constraints; j++) {
		// test needed otherwise the c_scaling can be 0, and division by 0 occurs
		if(m_c_scaling[j] > 0.) {
			solution_infeasibility += std::max(0.,(std::abs(c.at(j)) - c_tol.at(j))) / m_c_scaling[j];
		}
	}
	for(problem::base::c_size_type j=number_of_eq_constraints; j<prob_c_dimension; j++) {
		if(m_c_scaling[j] > 0.) {
			solution_infeasibility += std::max(0.,c.at(j) - c_tol.at(j)) / m_c_scaling[j];
		}
	}

	solution_infeasibility /= prob_c_dimension;

	return solution_infeasibility;
}
}}

//! @endcond

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::cstrs_self_adaptive)

