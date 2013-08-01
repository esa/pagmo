/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
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
#include "base.h"
#include "decompose.h"

namespace pagmo { namespace problem {


/**
 * Constructor
 *
 * @param[in] p base::problem to be decomposed
 * @param[in] method decomposition method (WEIGHTS, TCHEBYCHEFF, BI)
 * @param[in] weights the weight vector (by default is set to random weights)
 * @param[in] z reference point (used in Tchebycheff and Boundary Intersection (BI) methods, by default it is set to 0)
 *
 * @see For the uniform random generation of weights vector see Appendix 2 in "A. Jaszkiewicz - On the Performance of Multiple-Objective Genetic Local Search
on the 0/1 Knapsack Problem—A Comparative Experiment"
 * @see For the different decomposition methods see "Q. Zhang - MOEA/D: A Multiobjective Evolutionary Algorithm Based on Decomposition"

 */
decompose::decompose(const base & p, method_type method, const std::vector<double> & weights, const std::vector<double> & z):
	base((int)p.get_dimension(), // Ambiguous without the cast ...
		 p.get_i_dimension(),
		 1, //it transforms the problem into a single-objective problem
		 p.get_c_dimension(),
		 p.get_ic_dimension(),
		 p.get_c_tol()),
		 m_original_problem(p.clone()),
		 m_method(method),
		 m_weights(weights),
		 m_z(z)
{

	//0 - Check whether method is implemented
	if(m_method != WEIGHTED && m_method != TCHEBYCHEFF && m_method != BI) {
		pagmo_throw(value_error,"non existing decomposition method");
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
	}

	//2 - Checks whether the reference point has a dimension, if not, sets its default value
	if (m_z.size() == 0) {
		m_z = std::vector<double>((int)p.get_f_dimension(), 0.0); //by default m_z = (0, ..., 0)
	} else {
		//2.1 - Checks whether the reference point has lenght equal to the fitness size
		if (m_z.size() != p.get_f_dimension()) {
			pagmo_throw(value_error,"the the reference point vector must have equal length to the fitness size");
		}
	}

	//Setting the bounds according to the original problem
	set_lb(m_original_problem->get_lb());
	set_ub(m_original_problem->get_ub());
}

/// Copy Constructor. Performs a deep copy
decompose::decompose(const decompose &p):
	base((int)p.get_dimension(), // Ambiguous without the cast
		 p.get_i_dimension(),
		 p.get_f_dimension(),
		 p.get_c_dimension(),
		 p.get_ic_dimension(),
		 p.get_c_tol()),
		 m_original_problem(p.m_original_problem->clone()),
		 m_method(p.m_method),
		 m_weights(p.m_weights),
		 m_z(p.m_z)
		{set_lb(p.get_lb()); set_ub(p.get_ub());}

/// Clone method.
base_ptr decompose::clone() const
{
	return base_ptr(new decompose(*this));
}

/// Implementation of the objective function.
void decompose::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	fitness_vector fit(m_original_problem->get_f_dimension());
	m_original_problem->objfun(fit, x);

	if(m_method == WEIGHTED) {
		f[0] = 0.0;
		for(base::f_size_type i = 0; i < m_original_problem->get_f_dimension(); ++i) {
			f[0]+= m_weights[i]*fit[i];
		}
	} else if (m_method == TCHEBYCHEFF) {
		f[0] = m_weights[0] * fabs(fit[0] - m_z[0]);
		double tmp;
		for(base::f_size_type i = 0; i < m_original_problem->get_f_dimension(); ++i) {
			tmp = m_weights[i] * fabs(fit[i] - m_z[i]);
			if(tmp > f[0]) {
				f[0] = tmp;
			}
		}
	} else { //BI method
		const double THETA = 5.0;
		double d1 = 0.0;
		double weight_norm = 0.0;
		for(base::f_size_type i = 0; i < m_original_problem->get_f_dimension(); ++i) {
			d1 += (fit[i] - m_z[i]) * m_weights[i];
			weight_norm += pow(m_weights[i],2);
		}
		weight_norm = sqrt(weight_norm);
		d1 = fabs(d1)/weight_norm;

		double d2 = 0.0;
		for(base::f_size_type i = 0; i < m_original_problem->get_f_dimension(); ++i) {
			d2 += pow(fit[i] - (m_z[i] + d1*m_weights[i]/weight_norm), 2);
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
	oss << "\n\tDecomposition method: " << m_method << std::endl;
	oss << "\tWeight vector: " << m_weights << std::endl;
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

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::decompose);
