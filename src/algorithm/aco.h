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
#include "../problem/base_aco.h"
#include "../serialization.h"
#include "../types.h"
#include "base.h"

namespace pagmo { namespace algorithm {

/// Ant Colony Optimization (ACO)
/**
 * \image html ant.png "Ant Colony Optimization"
 * \image latex ant.png  "Ant Colony Optimization" width=3cm
 * Ant colony optimization (ACO) is a population-based metaheuristic that can be used to find approximate solutions to difficult combinatorial optimization problems. This implementation of ACO works on any constrained integer problem that extends the base_aco problem.
 *
 * NOTE: when called on mixed-integer problems ACO treats the continuous part as fixed and optimizes
 * the integer part.
 *
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 *
 * @see http://www.scholarpedia.org/article/Ant_colony_optimization
 */

class __PAGMO_VISIBLE aco: public base
{
public:
	aco(int iter = 1, double rho = 0.2);
	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;
protected:
	std::string human_readable_extra() const;
private:
	static void deposit_pheromone(std::vector<std::vector<std::vector<fitness_vector> > > &T, decision_vector &X, fitness_vector fit, double rho);
	static void selection_probability(std::vector<fitness_vector> &probability, std::vector<bool> &fComponents, std::vector<fitness_vector> &eta, std::vector<int> &selection, const pagmo::problem::base &prob);
	static void feasible_components(std::vector<bool> &fComponents,const pagmo::problem::base_aco &prob, decision_vector &X, problem::base::size_type xSize, double lb, double ub);
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<double &>(m_iter);
		ar & const_cast<double &>(m_rho);
	}
	// Number of iterations
	const double m_iter;
	const double m_rho;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::aco);

#endif // PAGMO_ALGORITHM_ACO_H
