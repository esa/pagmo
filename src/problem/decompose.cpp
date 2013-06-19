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
 * @param[in] weights the weight vector (by default is set to random weights)
 *
 * @see problem::base constructors.
 */
decompose::decompose(const base & p, const std::vector<double> & weights ):
	base((int)p.get_dimension(), // Ambiguous without the cast ...
		 p.get_i_dimension(),
		 1, //it transforms the problem into a single-objective problem
		 p.get_c_dimension(),
		 p.get_ic_dimension(),
		 p.get_c_tol()),
		 m_original_problem(p.clone()),
		 m_weights(weights)
{

	//1 - Checks whether the weight vector has a dimension, if not, sets its default value 
	if (m_weights.size() == 0) {
		//Initialise a random weight vector  
		rng_double m_drng = rng_generator::get<rng_double>();
		m_weights = std::vector<double>((int)p.get_f_dimension(), 0.0);
		double sum = 0;
		for(std::vector<double>::size_type i = 0; i<m_weights.size(); ++i) {
			m_weights[i] = boost::uniform_real<double>(0,1)(m_drng);
			sum += m_weights[i];
		} 
		for(std::vector<double>::size_type i = 0; i<m_weights.size(); ++i) {
			m_weights[i] /= sum;
		} 
	}
	
	//2 - Checks whether the weight vector sums to 1
	double sum = 0.0;
	for (std::vector<double>::size_type i=0; i<m_weights.size(); ++i) {
		sum += m_weights[i];
	}
	if (fabs(sum-1.0) > 1E-8) {
		pagmo_throw(value_error,"the weight vector should sum to 1 with a tolerance of E1-8");
	}
	
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
		 m_weights(p.m_weights)
		 {}

/// Clone method.
base_ptr decompose::clone() const
{
	return base_ptr(new decompose(*this));
}

/// Implementation of the objective function.
/// (Wraps over the original implementation with translated input x)
void decompose::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	fitness_vector fit(m_original_problem->get_f_dimension());
	m_original_problem->objfun(fit, x);
	
	f[0] = 0;
	for(base::f_size_type i = 0; i < m_original_problem->get_f_dimension(); ++i) {
		f[0]+= m_weights[i]*fit[i];	
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
