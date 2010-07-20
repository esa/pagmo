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

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "tsp.h"

namespace pagmo { namespace problem {

/// Constructor from vectors and maximum weight.
/**
 * Initialize weights of the edges from the matrix.
 *
 * @param[in] weights matrix of weights.
 */
tsp::tsp(const std::vector<std::vector<double> > &weights):
	base_aco(boost::numeric_cast<int>(weights[0].size()),boost::numeric_cast<int>(weights[0].size()),1,1,0,0),
	m_weights(weights) {
	set_lb(0);
	set_ub(weights[0].size()-1); //number of nodes in the graph -1 (we count from 0)
}

/**
 * Read the weight matrix from file
 * @param[in] ifile file containing the weights matrix
 */

/*
tsp::tsp(std::ifstream &ifile){
//TO IMPLEMENT
}*/

/** For tsp eta[k][i][j] represents the cost of having the node j in position k of the path and the node i in position k+1. 
 *  this represents the weight of the edg between i and j (distance from city i and j) and doesn't depends from k.
 */
void tsp::get_heuristic_information_matrix(std::vector<std::vector<std::vector<fitness_vector> > > &eta) const {
	for(std::vector<std::vector<std::vector<fitness_vector> > >::size_type k = 0; k < eta.size(); ++k) {
		for(std::vector<std::vector<fitness_vector> >::size_type i=0; i < eta[0].size(); ++i) {
			for(std::vector<fitness_vector>::size_type  j = 0; j < eta[0][0].size(); ++j) {
					eta[k][i][j][0] = m_weights[i][j];
			}
		}
	}

}
/*
 * Using a set check if the same node appears two times in the solution. In that case the solution is not feasible
 */
bool tsp::check_partial_feasibility(const decision_vector x) const{
	std::set<int> nodes;
	int node;
	for (size_type i = 0; i < x.size(); ++i) {
		node = boost::numeric_cast<int>(x[i]);
		if (node < get_lb()[i] || node > get_ub()[i]) {
			return false;
		}
		else {
			nodes.insert(node);
		}
	}
	return nodes.size() == x.size();
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
	std::set<int> nodes;
	int node;
	pagmo_assert(c.size() == 1 && x.size() == m_weights[0].size());
	for (size_type i = 0; i < get_dimension(); ++i) {
		node = boost::numeric_cast<int>(x[i]);
		if (node < get_lb()[i] || node > get_ub()[i]) {
			c[0] = 1;
			return;
		}
		else {
			nodes.insert(node);
		}
	}
	if (nodes.size() == get_dimension()) {
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
