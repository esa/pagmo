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
 * All the problems that can be solved using Ant Colony Optimization should extend this class
 *
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 */

class __PAGMO_VISIBLE base_aco : public base
{
	public:
		base_aco(int, int = 0, int = 1, int = 0, int = 0, const double & = 0);

		//write heuristic information matrix on T
		void get_heuristic_information_matrix(std::vector<std::vector<std::vector<fitness_vector> > > &eta);
		
		//check if a partial solution x is feasible. x.size() is less than problem length. It checks whether there is
		//at least one solution having x as a prefix that is feasible.
		bool check_partial_feasibility(decision_vector x); 
};

}} //namespaces

#endif // PAGMO_PROBLEM_BASE_ACO_H
