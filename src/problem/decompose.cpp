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
#include <boost/random/uniform_real.hpp>

#include "../exceptions.h"
#include "../types.h"
#include "../population.h"
#include "../rng.h"
#include "decompose.h"

namespace pagmo { namespace problem {


/**
 * Constructor
 *
 * @param[in] p base::problem to be decomposed
 * @param[in] method decomposition method (WEIGHTS, TCHEBYCHEFF, BI)
 * @param[in] weights the weight vector (by default is set to random weights)
 * @param[in] z ideal reference point (used in Tchebycheff and Boundary Intersection (BI) methods, by default it is set to 0)
 * @param[in] adapt_ideal if true it updates the ideal reference point each time the objective function is called checking if the computed fitness is better (assumes minimization)
 * @see For the uniform random generation of weights vector see Appendix 2 in "A. Jaszkiewicz - On the Performance of Multiple-Objective Genetic Local Search
on the 0/1 Knapsack Problem—A Comparative Experiment"
 * @see For the different decomposition methods see "Q. Zhang - MOEA/D: A Multiobjective Evolutionary Algorithm Based on Decomposition"

 */
decompose::decompose(const base & p, method_type method, const std::vector<double> & weights, const std::vector<double> & z, const bool adapt_ideal):
	base_meta(
		 p,
		 p.get_dimension(), // Ambiguous without the cast ...
		 p.get_i_dimension(),
		 1, //it transforms the problem into a single-objective problem
		 p.get_c_dimension(),
		 p.get_ic_dimension(),
		 p.get_c_tol()),
		 m_method(method),
		 m_weights(weights),
		 m_z(z),
		 m_adapt_ideal(adapt_ideal)
{

	//0 - Check whether method is implemented
	if(m_method != WEIGHTED && m_method != TCHEBYCHEFF && m_method != BI) {
		pagmo_throw(value_error,"non existing decomposition method");
	}

	if (p.get_f_dimension() == 1) {
		pagmo_throw(value_error,"decompose works only for multi-objective problems, you are trying to decompose a single objective one.");
	}

	//1 - Checks whether the weight vector has a dimension, if not, sets its default value
	if (m_weights.size() == 0) {
		//Initialise a random weight vector
		rng_double m_drng = rng_generator::get<rng_double>();
		m_weights = std::vector<double>((int)p.get_f_dimension(), 0.0);
		double sum = 0;
		for(std::vector<double>::size_type i = 0; i<m_weights.size(); ++i) {
			m_weights[i] = (1-sum) * (1 - pow(boost::uniform_real<double>(0,1)(m_drng), 1.0 / (m_weights.size() - i - 1)));
			sum += m_weights[i];
		}
	} else {

		//1.1 - Checks whether the weight has lenght equal to the fitness size
		if (m_weights.size() != p.get_f_dimension()) {
			pagmo_throw(value_error,"the weight vector must have length equal to the fitness size");
		}

		//1.2 - Checks whether the weight vector sums to 1
		double sum = 0.0;
		for (std::vector<double>::size_type i=0; i<m_weights.size(); ++i) {
			sum += m_weights[i];
		}
		if (fabs(sum-1.0) > 1E-8) {
			pagmo_throw(value_error,"the weight vector should sum to 1 with a tolerance of E1-8");
		}
		
		//1.4 - Checks that all weights are positive
		for (std::vector<double>::size_type i=0; i<m_weights.size(); ++i) {
			if (m_weights[i] < 0) {
				pagmo_throw(value_error,"the weight vector should contain only positive values");
			}
		}
	}

	//2 - Checks whether the reference point has a dimension, if not, sets its default value m_z = (0, ..., 0)
	if (m_z.size() == 0) {
		m_z = std::vector<double>((int)p.get_f_dimension(), 0.0); //by default m_z = (0, ..., 0)
	} else {
		//2.1 - Checks whether the reference point has lenght equal to the fitness size
		if (m_z.size() != p.get_f_dimension()) {
			pagmo_throw(value_error,"the the reference point vector must have equal length to the fitness size");
		}
	}
}

/// Clone method.
base_ptr decompose::clone() const
{
	return base_ptr(new decompose(*this));
}

///Implementation of the objective function
void decompose::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	fitness_vector fit(m_original_problem->get_f_dimension());
	compute_original_fitness(fit, x);
	compute_decomposed_fitness(f, fit,m_weights);
}

/// Gets the ideal point
fitness_vector decompose::get_ideal_point() const {
	return m_z;
}

/// Sets the ideal point
void decompose::set_ideal_point(const fitness_vector &f) {
	m_z = f;
}

/// Computes the original fitness
/**
 * Computes the original fitness of the multi-objective problem. It also updates the ideal point in case
 * m_adapt_ideal is true
 *
 * @param[out] f non-decomposed fitness vector
 * @param[in] x chromosome
 */
void decompose::compute_original_fitness(fitness_vector &f, const decision_vector &x) const {
	m_original_problem->objfun(f,x);
	if (m_adapt_ideal) {
		for (fitness_vector::size_type i=0; i<f.size(); ++i) {
			if (f[i] < m_z[i]) m_z[i] = f[i];
		}
	}
}

/// Computes the decomposed fitness
/**
 * Computes the decomposed fitness from the original multi-objective fitness using the data member m_weights
 *
 * @param f decomposed fitness vector
 * @param original_fit original multi-objective fitness vector
 */
void decompose::compute_decomposed_fitness(fitness_vector &f, const fitness_vector &original_fit) const
{
	compute_decomposed_fitness(f, original_fit,m_weights);
}

/// Computes the decomposed fitness
/**
 * Computes the decomposed fitness from the original multi-objective one and a weight vector
 *
 * @param[out] f decomposed fitness vector
 * @param[in] original_fit original multi-objective fitness vector
 * @param[in] weights weights vector
 */
void decompose::compute_decomposed_fitness(fitness_vector &f, const fitness_vector &original_fit, const fitness_vector &weights) const
{
	if ( (m_weights.size() != weights.size()) || (original_fit.size() != m_weights.size()) ) {
		pagmo_throw(value_error,"Check the sizes of input weights and fitness vector");
	}
	if(m_method == WEIGHTED) {
		f[0] = 0.0;
		for(base::f_size_type i = 0; i < m_original_problem->get_f_dimension(); ++i) {
			f[0]+= weights[i]*original_fit[i];
		}
	} else if (m_method == TCHEBYCHEFF) {
		f[0] = 0.0;
		double tmp,weight;
		for(base::f_size_type i = 0; i < m_original_problem->get_f_dimension(); ++i) {
			(weights[i]==0) ? (weight = 1e-4) : (weight = weights[i]); //fixes the numerical problem of 0 weights
			tmp = weight * fabs(original_fit[i] - m_z[i]);
			if(tmp > f[0]) {
				f[0] = tmp;
			}
		}
	} else { //BI method
		const double THETA = 5.0;
		double d1 = 0.0;
		double weight_norm = 0.0;
		for(base::f_size_type i = 0; i < m_original_problem->get_f_dimension(); ++i) {
			d1 += (original_fit[i] - m_z[i]) * weights[i];
			weight_norm += pow(weights[i],2);
		}
		weight_norm = sqrt(weight_norm);
		d1 = fabs(d1)/weight_norm;

		double d2 = 0.0;
		for(base::f_size_type i = 0; i < m_original_problem->get_f_dimension(); ++i) {
			d2 += pow(original_fit[i] - (m_z[i] + d1*weights[i]/weight_norm), 2);
		}
		d2 = sqrt(d2);

		f[0] = d1 + THETA * d2;
	}
}

std::string decompose::get_name() const
{
	return m_original_problem->get_name() + " [Decomposed]";
}

std::string decompose::human_readable_extra() const
{
	std::ostringstream oss;
	oss << m_original_problem->human_readable_extra() << std::endl;
	oss << "\n\tDecomposition method: ";
	switch(m_method){
		case WEIGHTED: {
			oss << "Weighted ";
			break;
		}
		case BI: {
			oss << "Boundary Interception ";
			break;
			}
		case TCHEBYCHEFF: {
			oss << "Tchebycheff ";
			break;
			}
	}
	oss << std::endl << "\tWeight vector: " << m_weights << std::endl;
	oss << "\tReference point: " << m_z << std::endl;
	return oss.str();
}


/**
 * Get the weights vector
 *
 * \return the weight vector
 */
const std::vector<double>& decompose::get_weights() const
{
	return m_weights;
}
}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::decompose)
