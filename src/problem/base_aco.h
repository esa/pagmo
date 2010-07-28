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

#ifndef PAGMO_PROBLEM_BASE_ACO_H
#define PAGMO_PROBLEM_BASE_ACO_H

#include "base.h"
#include <vector>

namespace pagmo{ namespace problem {

/// Base ACO.
/**
 *
 * All integer optimization problems must extend this class in order to be solved by Ant Colony Optimization.
 *
 * m_eta is the heuristic information matrix. It represent an a priori knowledge on the problem and has to be set in problem implementation.
 * m_eta[k][i][j] represents the cost of having the j-th value in position k of the chromosome
 * and the i-th value in position k+1. This type of info can later used by the algorithm to guide the optimization.
 * For example in the algorithm::aco it makes the ant prefer clever steps. In fact the probability for a particular step to be chosen
 * is proportional to the the product between the m_eta value of that step (heuristic information) and the amount of pheromone left
 * by previous ants on that step.
 * 
 *
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 */

class __PAGMO_VISIBLE base_aco : public base
{
	public:
		base_aco(int, int = 0, int = 0);

		/**
		 * Checks if a partial solution x is feasible. x.size() may be less than problem length.
		 * @returns true if there is at least one solution having x as a prefix that is feasible. False otherwise
		 */
		virtual bool check_partial_feasibility(const decision_vector &x) const;

		/**
		 * Gets the heuristic information matrix
		 * @returns const reference to m_eta: the heuristic information matrix
		 */
		const std::vector<std::vector<std::vector<fitness_vector> > > &get_heuristic_information_matrix() const;

	protected:
		/**
		 * Set the heuristic information matrix. This should be overridden by subclasses to set the proper heuristic information
		 * matrix for the problem
		 */
		virtual void set_heuristic_information_matrix();

		/**
		 * The heuristic information matrix
		 */
		std::vector<std::vector<std::vector<fitness_vector> > > m_eta;
		
		/**
		 * Allocate memory for the heuristic information matrix. That must be 
		 * called at the begining of each set_heuristic_information_matrix() implementation
		 */
		void create_heuristic_information_matrix();
};

}} //namespaces

#endif // PAGMO_PROBLEM_BASE_ACO_H
