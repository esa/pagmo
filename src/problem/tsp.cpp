/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
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

#include <boost/numeric/conversion/cast.hpp>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "tsp.h"

namespace pagmo { namespace problem {

static const double default_weights[5][5] = {
	{0,1,2,3,4},
	{1,0,2.236067,4,4.123105},
	{2,2.236067,0,3.605551,6},
	{3,4,3.605551,0,5},
	{4,4.123105,6,5,0}
};

/// Default constructor.
/**
 * This constructor will build a TSP instance with 5 cities in the following positions:
 *
 * - city 1: (0,0)
 * - city 2: (1,0)
 * - city 3: (0,2)
 * - city 4: (-3,0)
 * - city 5: (0,-4)
 */
tsp::tsp():base_aco(5,1,0)
{
	m_weights.resize(5);
	for (int i = 0; i < 5; ++i) {
		m_weights[i].resize(5);
		for (int j = 0; j < 5; ++j) {
			m_weights[i][j] = default_weights[i][j];
		}
	}
	set_lb(0);
	set_ub(4); //number of nodes in the graph -1 (we count from 0)
	set_heuristic_information_matrix();
}

/// Constructor from vectors and maximum weight.
/**
 * Initialize weights of the edges (city distances) from the matrix.
 *
 * @param[in] weights matrix of distances between cities.
 */
tsp::tsp(const std::vector<std::vector<double> > &weights):
	base_aco(boost::numeric_cast<int>(weights[0].size()),1,0), m_weights(weights) {
	
	//Check weights matrix
	for(problem::base_aco::size_type i = 0; i < m_weights.size(); ++i) {
		if (m_weights[i].size() != m_weights.size()) {
			pagmo_throw(value_error,"Weights matrix must be a square matrix!");		
		}
		if(m_weights[i][i] != 0) {
			pagmo_throw(value_error,"Weights matrix must have 0's on the diagonal!");
		}
		for (problem::base_aco::size_type j = 0; j < m_weights[i].size(); ++j) {
			if(m_weights[i][j] != m_weights[j][i]) {
				pagmo_throw(value_error,"Weights matrix must be a simmetric matrix!");	
			}
		}
	}

	set_lb(0);
	set_ub(weights[0].size()-1); //number of nodes in the graph -1 (we count from 0)
	set_heuristic_information_matrix();
}


/** For tsp eta[k][i][j] represents the cost of having the node j in position k of the path and the node i in position k+1. 
 *  this represents the weight of the edge between i and j (distance from city i and j) and doesn't depends from k.
 */
void tsp::set_heuristic_information_matrix() {
	//allocates the memory for eta.
	create_heuristic_information_matrix();

	for(std::vector<std::vector<std::vector<fitness_vector> > >::size_type k = 0; k < m_eta.size(); ++k) {
		for(std::vector<std::vector<fitness_vector> >::size_type i=0; i < m_eta[0].size(); ++i) {
			for(std::vector<fitness_vector>::size_type  j = 0; j < m_eta[0][0].size(); ++j) {
					m_eta[k][i][j][0] = m_weights[i][j];
			}
		}
	}

}
/*
 * Chek if a node appears two times in the solution. In that case the solution is not feasible
 */
bool tsp::check_partial_feasibility(const decision_vector &x) const{
	if (x.size() > get_i_dimension()) {
		pagmo_throw(value_error,"Invalid chromosome length for partial feasibility check.");
	}

	m_tmpDecisionVector = x;
	std::sort(m_tmpDecisionVector.begin(), m_tmpDecisionVector.end());

	return 	std::unique(m_tmpDecisionVector.begin(), m_tmpDecisionVector.end()) == m_tmpDecisionVector.end();

}
 

/// Clone method.
base_ptr tsp::clone() const
{
	return base_ptr(new tsp(*this));
}

/// Implementation of the objective function.
void tsp::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == 1);
	pagmo_assert(x.size() == get_dimension() && x.size() == m_weights[0].size());
	f[0] = 0;
	for (size_type i = 1; i < get_dimension(); ++i) {
			f[0] += m_weights[boost::numeric_cast<int>(x[i-1])][boost::numeric_cast<int>(x[i])];
	}
	f[0] += m_weights[boost::numeric_cast<int>(x[get_dimension()-1])][boost::numeric_cast<int>(x[0])];
}

/// Re-implement constraint computation,
//We check whether we have selected all the nodes (the decision vector has to be a permutation of the set of nodes).
//The constraint is positive (not satisfied) if we have selected more than once the same node or equivalently not all 
//the nodes have been selected
void tsp::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	if (check_partial_feasibility(x)) {
		c[0] = 0;
	}
	else {
		c[0] = 1;
	}
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the weights matrix.
 */
std::string tsp::human_readable_extra() const
{
	std::ostringstream oss;
	for(problem::base::size_type i=0; i < m_weights[0].size(); ++i) {
		oss << "Weights: " << m_weights[i] << '\n';
	}
	return oss.str();
}

std::string tsp::get_name() const
{
	return "Travelling Salesman Problem";
}

}
}
