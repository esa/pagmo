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

#ifndef PAGMO_ALGORITHM_ACO_H
#define PAGMO_ALGORITHM_ACO_H

#include "../config.h"
#include "base.h"
#include "../problem/base.h"


namespace pagmo { namespace algorithm {

/// Ant Colony Optimization (ACO)
/**
 * Ant colony optimization (ACO) is a population-based metaheuristic that can be used to find approximate solutions to difficult combinatorial optimization problems. 
 *
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 *
 * @see http://www.scholarpedia.org/article/Ant_colony_optimization
 */

class __PAGMO_VISIBLE aco: public base
{
public:
	aco(int iter, double rho = 0.2);
	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;
protected:
	std::string human_readable_extra() const;
private:
	static void deposit_pheromone(std::vector<std::vector<std::vector<fitness_vector> > > &T, decision_vector &X, fitness_vector fit, double rho);
	static void selection_probability(std::vector<fitness_vector> probability, std::vector<int> &selection, const pagmo::problem::base &prob, population::size_type NP);
	// Number of iterations
	const double m_iter;
	const double m_rho;
};

}} //namespaces

#endif // PAGMO_ALGORITHM_ACO_H
