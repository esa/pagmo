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
#include "cstrs_co_evolution.h"

///Doxygen will ignore whatever is in //! @cond
//! @cond

namespace pagmo { namespace problem {

/**
 * Constructor of co-evolution problem from a constrained problem
 *
 * Note: This problem is not intended to be used by itself, only by the
 * pagmo::algorithm::cstrs_co_evolution algorithm 
 *
 * @param[in] problem base::problem to be modified to use a co-evolution
 * as constraints handling technique.
 * @param[in] method method_type to used for the co-evolution constraints
 * handling technique. Three posssibililties are available: SIMPLE,
 * SPLIT_NEQ_EQ and SPLIT_CONSTRAINTS. The simple one is the original
 * version of the Coello/He implementation (one penalty coefficient weights 
 * the sum of the constraints violation, one the number of violated constraints). 
 * The SPLIT_NEQ_EQ, splits the equalities and inequalities constraints in two different sets for the
 * penalty weigths, containing respectively inequalities and equalities
 * weigths. The SPLIT_CONSTRAINTS splits the constraints in M set of weigths
 * with M the number of constraints.
 *
 */
cstrs_co_evolution::cstrs_co_evolution(const base &problem, const algorithm::cstrs_co_evolution::method_type method):
	base_meta(
		 problem, 
		 problem.get_dimension(),
		 problem.get_i_dimension(),
		 problem.get_f_dimension(),
		 0,
		 0,
		 std::vector<double>()),
	m_penalty_coeff(),
	m_method(method),
	m_map_fitness(),
	m_map_constraint(),
	m_decision_vector_hash()
{
	if(m_original_problem->get_c_dimension() <= 0){
		pagmo_throw(value_error,"The original problem has no constraints.");
	}

	// check that the dimension of the problem is 1
	if (m_original_problem->get_f_dimension() != 1) {
		pagmo_throw(value_error,"The original fitness dimension of the problem must be one, multi objective problems can't be handled with co-evolution meta problem.");
	}

	std::fill(m_penalty_coeff.begin(),m_penalty_coeff.end(),0.);
}

cstrs_co_evolution::cstrs_co_evolution(const base &problem, const population& pop, const algorithm::cstrs_co_evolution::method_type method):
	base_meta(
		 problem, 
		 problem.get_dimension(),
		 problem.get_i_dimension(),
		 problem.get_f_dimension(),
		 0,
		 0,
		 std::vector<double>()),
	m_penalty_coeff(),
	m_method(method),
	m_map_fitness(),
	m_map_constraint(),
	m_decision_vector_hash()
{
	if(m_original_problem->get_c_dimension() <= 0){
		pagmo_throw(value_error,"The original problem has no constraints.");
	}

	// check that the dimension of the problem is 1
	if (m_original_problem->get_f_dimension() != 1) {
		pagmo_throw(value_error,"The original fitness dimension of the problem must be one, multi objective problems can't be handled with co-evolution meta problem.");
	}
	if(problem != pop.problem()) {
		pagmo_throw(value_error,"The problem linked to the population is not the same as the problem given in argument.");
	}

	std::fill(m_penalty_coeff.begin(),m_penalty_coeff.end(),0.);

	m_map_fitness.clear();
	m_map_constraint.clear();
	// store f and c in maps depending on the the x hash
	for(population::size_type i=0; i<pop.size(); i++) {
		const population::individual_type &current_individual = pop.get_individual(i);
		// m_map_fitness.insert(std::pair<std::size_t, fitness_vector>(m_decision_vector_hash(current_individual.cur_x),current_individual.cur_f));
		m_map_fitness[m_decision_vector_hash(current_individual.cur_x)]=current_individual.cur_f;
		m_map_constraint[m_decision_vector_hash(current_individual.cur_x)]=current_individual.cur_c;
	}

}

/// Clone method.
base_ptr cstrs_co_evolution::clone() const
{
	return base_ptr(new cstrs_co_evolution(*this));
}

/// Implementation of the objective function.
/// (Wraps over the original implementation)
/**
 *  Returns the penalized fitness if the decision vector penalize with the penalty
 *  coefficient given.
 */
void cstrs_co_evolution::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	std::map<std::size_t, fitness_vector>::const_iterator it_f;

	it_f = m_map_fitness.find(m_decision_vector_hash(x));
	if(it_f != m_map_fitness.end()) {
		f = it_f->second;
	} else {
		m_original_problem->objfun(f, x);
	}

	std::vector<double> sum_viol;
	std::vector<int> num_viol;

	compute_penalty(sum_viol,num_viol,x);

	// assuming minimization
	switch(m_method)
	{
	case algorithm::cstrs_co_evolution::SIMPLE:
	{
		f[0] += sum_viol.at(0) * m_penalty_coeff.at(0) + double(num_viol.at(0)) * m_penalty_coeff.at(1);
		break;
	}
	case algorithm::cstrs_co_evolution::SPLIT_NEQ_EQ:
	{
		f[0] += sum_viol.at(0) * m_penalty_coeff.at(0) + double(num_viol.at(0)) * m_penalty_coeff.at(1);
		f[0] += sum_viol.at(1) * m_penalty_coeff.at(2) + double(num_viol.at(1)) * m_penalty_coeff.at(3);
		break;
	}
	case algorithm::cstrs_co_evolution::SPLIT_CONSTRAINTS:
	{
		int c_dimension = m_original_problem->get_c_dimension();
		for(int i=0; i<c_dimension; i++) {
			f[0] += sum_viol.at(i) * m_penalty_coeff.at(i*2) + double(num_viol.at(i)) * m_penalty_coeff.at(1+i*2);
		}
		break;
	}
	}
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the type of constraint handling
 */
std::string cstrs_co_evolution::human_readable_extra() const
{
	std::ostringstream oss;
	oss << m_original_problem->human_readable_extra() << std::endl;
	oss << "\n\tConstraints handled with co-evolution method ";
	oss << std::endl;
	return oss.str();
}

std::string cstrs_co_evolution::get_name() const
{
	return m_original_problem->get_name() + " [cstrs_co_evolution]";
}

/// Updates the fitness information based on the population.
/**
 *  By calling this method, penalties coefficients and
 *  fitnesses for the whole given population are computed.
 */
void cstrs_co_evolution::set_penalty_coeff(const std::vector<double> &penalty_coeff)
{
	switch(m_method)
	{
	case algorithm::cstrs_co_evolution::SIMPLE:
	{
		if(penalty_coeff.size() != 2) {
			pagmo_throw(value_error,"The size of the penalty coefficient vector is not 2.");
		}
		break;
	}
	case algorithm::cstrs_co_evolution::SPLIT_NEQ_EQ:
	{
		if(penalty_coeff.size() != 4) {
			pagmo_throw(value_error,"The size of the penalty coefficient vector is not 4.");
		}
		break;
	}
	case algorithm::cstrs_co_evolution::SPLIT_CONSTRAINTS:
		if(penalty_coeff.size() != 2*m_original_problem->get_c_dimension()) {
			pagmo_throw(value_error,"The size of the penalty coefficient vector is not 2*number constraints");
		}
		break;
	default: 
		pagmo_throw(value_error,"The constraints co-evolutionary method is not valid.");
		break;
	}
	m_penalty_coeff = penalty_coeff;
}

/// Returns the size of the penalty coefficient the problem expects
/// depending on the method used.
int cstrs_co_evolution::get_penalty_coeff_size() {
	switch(m_method)
	{
	case algorithm::cstrs_co_evolution::SIMPLE:
	{
		return 2;
		break;
	}
	case algorithm::cstrs_co_evolution::SPLIT_NEQ_EQ:
	{
		return 4;
		break;
	}
	case algorithm::cstrs_co_evolution::SPLIT_CONSTRAINTS:
		return 2*m_original_problem->get_c_dimension();
		break;
	default: 
		pagmo_throw(value_error,"The constraints co-evolutionary method is not valid.");
		break;
	}
}

/// Computes the penalty depending on the provided penalty coefficient.
/**
 * Updates the solution infeasibility vector with the population given.
 * @param[in,out] std::vector<double> sum_viol sum of violation vector.
 * @param[in,out] std::vector<double> num_viol number of violation vector.
 * @param[in] decision_vector x.
 */
void cstrs_co_evolution::compute_penalty(std::vector<double> &sum_viol, std::vector<int> &num_viol, const decision_vector &x) const
{
	// get the constraints dimension
	problem::base::c_size_type prob_c_dimension = m_original_problem->get_c_dimension();
	constraint_vector c(prob_c_dimension, 0.);
	constraint_vector c_vio(prob_c_dimension, 0.);
	problem::base::c_size_type number_of_eq_constraints = prob_c_dimension - m_original_problem->get_ic_dimension();

	const std::vector<double> &c_tol = m_original_problem->get_c_tol();

	// updates the current constraint vector
	std::map<std::size_t, constraint_vector>::const_iterator it_c;

	it_c = m_map_constraint.find(m_decision_vector_hash(x));
	if(it_c != m_map_constraint.end()) {
		c = it_c->second;
	} else {
		m_original_problem->compute_constraints(c,x);
	}
	

	// sets the right definition of the constraints violation
	for(problem::base::c_size_type j=0; j<prob_c_dimension; j++) {
		if(j<number_of_eq_constraints){
			c_vio[j] = std::abs(c.at(j)) - c_tol.at(j);
		}
		else{
			c_vio[j] = c.at(j) - c_tol.at(j);
		}
	}
	for(problem::base::c_size_type j=0; j<prob_c_dimension; j++) {
		c_vio[j] = std::max(0., c_vio.at(j));
	}

	// updates the vectors depending on the method
	switch(m_method)
	{
	case algorithm::cstrs_co_evolution::SIMPLE:
	{
		sum_viol.resize(1);
		num_viol.resize(1);
		std::fill(sum_viol.begin(),sum_viol.end(),0.);
		std::fill(num_viol.begin(),num_viol.end(),0.);

		// update sum_num_viol
		for(problem::base::c_size_type j=0; j<prob_c_dimension; j++) {
			sum_viol[0] += c_vio.at(j);
		}

		for(problem::base::c_size_type j=0; j<prob_c_dimension; j++) {
			if(!m_original_problem->test_constraint(c, j)) {
				num_viol[0] += 1;
			}
		}
		break;
	}
	case algorithm::cstrs_co_evolution::SPLIT_NEQ_EQ:
	{
		sum_viol.resize(2);
		num_viol.resize(2);
		std::fill(sum_viol.begin(),sum_viol.end(),0.);
		std::fill(num_viol.begin(),num_viol.end(),0.);

		// update sum_num_viol
		for(problem::base::c_size_type j=0; j<number_of_eq_constraints; j++) {
			sum_viol[0] += c_vio.at(j);
		}
		for(problem::base::c_size_type j=number_of_eq_constraints; j<prob_c_dimension; j++) {
			sum_viol[1] += c_vio.at(j);
		}
		for(problem::base::c_size_type j=0; j<number_of_eq_constraints; j++) {
			if(!m_original_problem->test_constraint(c, j)) {
				num_viol[0] += 1;
			}
		}
		for(problem::base::c_size_type j=number_of_eq_constraints; j<prob_c_dimension; j++) {
			if(!m_original_problem->test_constraint(c, j)) {
				num_viol[1] += 1;
			}
		}
		break;
	}
	case algorithm::cstrs_co_evolution::SPLIT_CONSTRAINTS:
	{
		sum_viol.resize(prob_c_dimension);
		num_viol.resize(prob_c_dimension);

		std::fill(sum_viol.begin(),sum_viol.end(),0.);
		std::fill(num_viol.begin(),num_viol.end(),0.);

		// update sum_num_viol
		for(problem::base::c_size_type j=0; j<prob_c_dimension; j++) {
			sum_viol[j] += c_vio.at(j);
		}
		for(problem::base::c_size_type j=0; j<prob_c_dimension; j++) {
			if(!m_original_problem->test_constraint(c, j)) {
				num_viol[j] += 1;
			}
		}
		break;
	}
	}
}

/**
 * Constructor of co-evolution problem using initial constrained problem
 *
 * Note: This problem is not inteded to be used by itself. Instead use the
 * co-evolution algorithm if you want to solve constrained problems.
 *
 * @param[in] problem base::problem to be modified to use a co-evolution
 * as constraints handling technique.
 *
 */
cstrs_co_evolution_penalty::cstrs_co_evolution_penalty(const base &problem, int dimension, int size):
	base(dimension,
		 problem.get_i_dimension(),
		 1,
		 0,
		 0,
		 0.),
	m_original_problem(problem.clone()),
	m_pop_2_x_vector(size, decision_vector(0)),
	m_feasible_count_vector(size,0),
	m_feasible_fitness_sum_vector(size,0.0),
	m_max_feasible_fitness(0.),
	m_total_sum_viol(size,0.0),
	m_total_num_viol(size,0)
{
	if(m_original_problem->get_c_dimension() <= 0){
		pagmo_throw(value_error,"The original problem has no constraints.");
	}

	// check that the dimension of the problem is 1
	if (m_original_problem->get_f_dimension() != 1) {
		pagmo_throw(value_error,"The original fitness dimension of the problem must be one, multi objective problems can't be handled with co evolution meta problem.");
	}

	set_bounds(0.,10000.);
}

/// Copy Constructor. Performs a deep copy
cstrs_co_evolution_penalty::cstrs_co_evolution_penalty(const cstrs_co_evolution_penalty &prob):
	base(prob),
	m_original_problem(prob.m_original_problem->clone()),
	m_pop_2_x_vector(prob.m_pop_2_x_vector),
	m_feasible_count_vector(prob.m_feasible_count_vector),
	m_feasible_fitness_sum_vector(prob.m_feasible_fitness_sum_vector),
	m_max_feasible_fitness(prob.m_max_feasible_fitness),
	m_total_sum_viol(prob.m_total_sum_viol),
	m_total_num_viol(prob.m_total_num_viol)
{
	set_bounds(prob.get_lb(),prob.get_ub());
}

/// Clone method.
base_ptr cstrs_co_evolution_penalty::clone() const
{
	return base_ptr(new cstrs_co_evolution_penalty(*this));
}

/// Implementation of the objective function.
/// (Wraps over the original implementation)
/**
 *  Returns the averaged fitness penalized with feasabiltity informaitions.
 */
void cstrs_co_evolution_penalty::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	// search where the x vector is located
	int position=-1;
	for(std::vector<decision_vector>::size_type i=0; i<m_pop_2_x_vector.size(); i++) {
		if(m_pop_2_x_vector.at(i) == x) {
			position = i;
			break;
		}
	}

	if(position != -1) {
		// computes the fitness for population 2
		// case the solution containts at least one feasible individual
		if(m_feasible_count_vector.at(position) > 0) {
			f[0] = m_feasible_fitness_sum_vector[position] / double(m_feasible_count_vector[position]) - double(m_feasible_count_vector[position]);
		} else {
			f[0] = m_max_feasible_fitness +
					m_total_sum_viol.at(position)/double(m_total_num_viol.at(position)) -
					double(m_total_num_viol.at(position));
		}
	}
	else {
		// what to do?
		f[0] = m_max_feasible_fitness;
	}
}

/// Implementation of fitness vectors comparison.
/**
 * @brief compare_fitness_impl calls the compare_fitness method of the original problem.
 * @return true if v_f1 is dominating v_f2, false otherwise.
 */
bool cstrs_co_evolution_penalty::compare_fitness_impl(const fitness_vector &v_f1, const fitness_vector &v_f2) const
{
	return m_original_problem->compare_fitness(v_f1,v_f2);
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the type of constraint handling
 */
std::string cstrs_co_evolution_penalty::human_readable_extra() const
{
	std::ostringstream oss;
	oss << m_original_problem->human_readable_extra() << std::endl;
	oss << "\n\tConstraints handled with co-evolution method ";
	oss << std::endl;
	return oss.str();
}

std::string cstrs_co_evolution_penalty::get_name() const
{
	return m_original_problem->get_name() + " [cstrs_co_evolution_2]";
}

/// Updates the fitness information based on the population.
/**
 *  By calling this method, penalties coefficients and
 *  fitnesses for the whole given population are computed.
 */
void cstrs_co_evolution_penalty::update_penalty_coeff(population::size_type &index, const decision_vector &pop_2_x, const population  &pop_1)
{
	m_pop_2_x_vector.at(index) = pop_2_x;

	population::size_type pop_1_size = pop_1.size();

	m_max_feasible_fitness = 0.;
	// computes the total sum_viol, num_viol...

	m_total_sum_viol[index] = 0.;
	m_total_num_viol[index] = 0;

	// computes the number of feasible solutions and their sum for the current population
	int feasible_count = 0;
	double feasible_fitness_sum = 0.;

	double sum_viol_temp = 0.;
	int num_viol_temp = 0;
	for(population::size_type i=0; i<pop_1_size; i++) {
		const fitness_vector &current_f = pop_1.get_individual(i).cur_f;
		const decision_vector &current_c = pop_1.get_individual(i).cur_c;

		if(m_original_problem->feasibility_c(current_c)) {
			feasible_count++;
			feasible_fitness_sum += current_f[0];
			
			// computes max feasible fitness
			if(i==0) {
				m_max_feasible_fitness = current_f[0];
			} else {
				if(m_max_feasible_fitness < current_f[0]) {
					m_max_feasible_fitness = current_f[0];
				}
			}
		}

		compute_penalty(sum_viol_temp,num_viol_temp,current_c);
		m_total_sum_viol[index] += sum_viol_temp;
		m_total_num_viol[index] += num_viol_temp;
	}
	m_feasible_count_vector[index] = feasible_count;
	m_feasible_fitness_sum_vector[index] = feasible_fitness_sum;
}

/// Computes the penalty depending on the provided penalty coefficient.
/**
 * Updates the solution infeasibility vector with the population given.
 * @param[in,out] std::vector<double> sum_viol sum of violation vector.
 * @param[in,out] std::vector<double> num_viol number of violation vector.
 * @param[in] decision_vector x.
 */
void cstrs_co_evolution_penalty::compute_penalty(double &sum_viol, int &num_viol, const constraint_vector &c) const
{
	// get the constraints dimension
	problem::base::c_size_type prob_c_dimension = m_original_problem->get_c_dimension();
	problem::base::c_size_type number_of_eq_constraints = prob_c_dimension - m_original_problem->get_ic_dimension();

	const std::vector<double> &c_tol = m_original_problem->get_c_tol();

	// update sum_num_viol
	sum_viol = 0.;
	for(problem::base::c_size_type j=0; j<number_of_eq_constraints; j++) {
		sum_viol += std::max(0.,std::abs(c.at(j)) - c_tol.at(j));
	}

	for(problem::base::c_size_type j=number_of_eq_constraints; j<prob_c_dimension; j++) {
		sum_viol += std::max(0.,c.at(j));
	}

	num_viol = 0;
	for(problem::base::c_size_type j=0; j<prob_c_dimension; j++) {
		if(!m_original_problem->test_constraint(c, j)) {
			num_viol += 1;
		}
	}
}
}}

//! @endcond

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::cstrs_co_evolution)
BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::cstrs_co_evolution_penalty)
