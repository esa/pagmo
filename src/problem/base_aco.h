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
 * CR - More documentation is needed here. In particular what is m_eta? How should a user know? WHat
 * is check_partial_feasibility working? All that a user need to know in order to be able to use the calss goes here
 *
 * CR Example: m_eta[k][i][j] represents the cost of having the j-th value in position k of the chromosome
 * and the i-th value in position k+1. This type of info can later used by the algorithm to guide the optimization.
 * For example in the algorithm::aco it makes the ant prefer clever steps.....
 *
 *
 * All the problems that can be solved using Ant Colony Optimization should extend this class
 *
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 */

class __PAGMO_VISIBLE base_aco : public base
{
	public:
		//CR - being limited to integer programming single objective, global dim = int dim, ctol does not make sense
		//and nf = 1. The way you implemented it was highly risky as allowing n to be different to ni. Also the default values you
		//used did not make any sense
		base_aco(int, int = 0, int = 0);
		//check if a partial solution x is feasible. x.size() is less than problem length. It checks whether there is
		//at least one solution having x as a prefix that is feasible.

		//CR - const. ref. otherwise you waste memory allocation!!!
		virtual bool check_partial_feasibility(const decision_vector &x) const;

		// CR - All public and protected members need to be documented properly example:
		/// Gets the heuristic information matrix
		/**
		 * @returns const reference to m_eta: the heuristic information matrix
		 */

		 //CR - virtual de che? 2- perche allochi memoria????!?!?!?!?!? Ritorna const ref.
		const std::vector<std::vector<std::vector<fitness_vector> > > &get_heuristic_information_matrix() const;	
	protected:
		// CR - consistency in names is important
		virtual void set_heuristic_information_matrix();
		std::vector<std::vector<std::vector<fitness_vector> > > m_eta;
	private: // CR - questi metodi non vengono usati nella classe derivata dunque private.
		void create_heuristic_information_matrix();
};

}} //namespaces

#endif // PAGMO_PROBLEM_BASE_ACO_H
