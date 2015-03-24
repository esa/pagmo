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
#include "base.h"
#include "antibodies_problem.h"

namespace pagmo { namespace problem {
/**
 * Constructor of antibodies meta-problem
 *
 * Note: This problem is not intended to be used by itself. Instead use the
 * cstrs_immune_system algorithm if you want to solve constrained problems.
 *
 * @param[in] problem base::problem to be used to set up the boundaries
 * @param[in] method method_type to used for the distance computation.
 * Two posssibililties are available: HAMMING, EUCLIDEAN.
 */
antibodies_problem::antibodies_problem(const base &problem, const algorithm::cstrs_immune_system::distance_method_type method):
	base((int)problem.get_dimension(),
		 problem.get_i_dimension(),
		 problem.get_f_dimension(),
		 0,
		 0,
		 0.),
	m_original_problem(problem.clone()),
	m_pop_antigens(),
	m_method(method)
{
	if(m_original_problem->get_c_dimension() <= 0){
		pagmo_throw(value_error,"The original problem has no constraints.");
	}

	// check that the dimension of the problem is 1
	if(m_original_problem->get_f_dimension() != 1) {
		pagmo_throw(value_error,"The original fitness dimension of the problem must be one, multi objective problems can't be handled with co-evolution meta problem.");
	}

	// encoding for hamming distance
	m_bit_encoding = 25;
	m_max_encoding_integer = int(std::pow(2., m_bit_encoding));

	set_bounds(m_original_problem->get_lb(),m_original_problem->get_ub());
}

/// Copy Constructor. Performs a deep copy
antibodies_problem::antibodies_problem(const antibodies_problem &prob):
	base((int)prob.get_dimension(),
		 prob.get_i_dimension(),
		 prob.get_f_dimension(),
		 prob.get_c_dimension(),
		 prob.get_ic_dimension(),
		 prob.get_c_tol()),
	m_original_problem(prob.m_original_problem->clone()),
	m_pop_antigens(prob.m_pop_antigens),
	m_method(prob.m_method),
	m_bit_encoding(prob.m_bit_encoding),
	m_max_encoding_integer(prob.m_max_encoding_integer)
{
	set_bounds(m_original_problem->get_lb(),m_original_problem->get_ub());
}

/// Clone method.
base_ptr antibodies_problem::clone() const
{
	return base_ptr(new antibodies_problem(*this));
}

/// Implementation of the objective function.
/// (Wraps over the original implementation)
/**
 *  Returns the distance of the given decision_vector to the antigenes population.
 */
void antibodies_problem::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	// compute the fitness which is the distance to the current population
	f[0] = compute_distance(x);
}

/// Implementation of fitness vectors comparison.
/**
 * @brief compare_fitness_impl calls the compare_fitness method of the original problem.
 * @return true if v_f1 is dominating v_f2, false otherwise.
 */
bool antibodies_problem::compare_fitness_impl(const fitness_vector &v_f1, const fitness_vector &v_f2) const
{
	return m_original_problem->compare_fitness(v_f1,v_f2);
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the type of constraint handling
 */
std::string antibodies_problem::human_readable_extra() const
{
	std::ostringstream oss;
	oss << m_original_problem->human_readable_extra() << std::endl;
	oss << "\n\tWith antibodies method";
	oss << std::endl;
	return oss.str();
}

std::string antibodies_problem::get_name() const
{
	return m_original_problem->get_name() + " [antibodies_problem]";
}

/// Updates the antigens population used to compute the fitness.
/**
 *
 */
void antibodies_problem::set_antigens(const std::vector<decision_vector> &pop_antigens)
{
	if(pop_antigens.size() == 0) {
		pagmo_throw(value_error,"The size of the antigens population must be different from 0.");
	}
	m_pop_antigens = pop_antigens;
}

/// Computes the distance to the antibodies population.
/**
 *
 */
double antibodies_problem::compute_distance(const decision_vector &x) const {
	double distance = 0.;

	// hamming distance

	switch(m_method) {
	case(algorithm::cstrs_immune_system::HAMMING): {
		const decision_vector &lb = get_lb();
		const decision_vector &ub = get_ub();

		for(decision_vector::size_type i=0; i<x.size(); i++) {

			std::vector<int> current_binary_gene = double_to_binary(x.at(i), lb.at(i), ub.at(i));

			for(decision_vector::size_type j=0; j<m_pop_antigens.size(); j++) {

				std::vector<int> antigens_binary_gene = double_to_binary((m_pop_antigens.at(j)).at(i), lb.at(i), ub.at(i));

				for(std::vector<int>::size_type k=0; k<antigens_binary_gene.size(); k++) {
					distance += antigens_binary_gene.at(k) && current_binary_gene.at(k);
				}
			}
		}

		// we take the inverse of the distance as the measure we need is
		// how close the x is from the antigen population
		// which means that we need to maximize the ressemblance
		distance = - distance;
		break;
	}
	case(algorithm::cstrs_immune_system::EUCLIDEAN): {
		for(decision_vector::size_type j=0; j<m_pop_antigens.size(); j++) {

			double euclid = 0.;
			const decision_vector &antigen_decision = m_pop_antigens.at(j);

			for(decision_vector::size_type i=0; i<x.size(); i++) {
				euclid += std::pow(x.at(i) - antigen_decision.at(i),2);
			}
			distance += std::sqrt(euclid);
		}

		break;
	}
	}

	return distance;
}

/// Convert a number to its binary representation.
/**
 *
 * @param[in] number: number to convert in binary representation.
 * @param[in] lb: lower bound for the encoding.
 * @param[in] ub: upper bound for the encoding.
 * @param[out] chromosome containing the binary representation of the number.
 */
std::vector<int> antibodies_problem::double_to_binary(const double &number, const double &lb, const double &ub) const
{
	// initialize the vector of size m_bit_encoding
	std::vector<int> binary(m_bit_encoding, 0);

	// convert the current number into its binary representation considering the domain
	// available
	int temp_number = (number - lb) * (m_max_encoding_integer - 1) / (ub - lb) ;

	// store the binary number
	int position = 0;
	while (temp_number!=0)
	{
		if( (temp_number % 2) == 0 ) {
			binary[position] = 0;
		} else {
			binary[position] = 1;
		}
		temp_number = (int)std::floor(temp_number/2);
		position++;
	}
	// reverse the order as this algorithm provides the reverse binary reprentation
	std::reverse(binary.begin(),binary.end());

	return binary;
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::antibodies_problem)
