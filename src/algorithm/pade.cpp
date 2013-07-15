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


#include <string>
#include <vector>
#include <algorithm>

#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../archipelago.h"
#include "../island.h"
#include "../population.h"
#include "../topology/unconnected.h"
#include "../problem/decompose.h"
#include "../types.h"
#include "base.h"
#include "pade.h"

namespace pagmo { namespace algorithm {
/// Constructor
 /**
 * Constructs a PaDe algorithm
 *
 * @param[in] gen Number of generations to evolve.
 * @param[in] solver the algorithm to solve the single objective problems.
 *
 * @throws value_error if gen is negative
 */
pade::pade(int gen, const pagmo::algorithm::base & solver):base(),m_gen(gen),m_solver(solver.clone())
{
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
    }
}

/// Clone method.
base_ptr pade::clone() const
{
    return base_ptr(new pade(*this));
}

/// Evolve implementation.
/**
 * Run the PaDe algorithm for the number of generations specified in the constructors.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void pade::evolve(population &pop) const
{
    // Get out if there is nothing to do.
    if (m_gen == 0) {
        return;
    }

    pagmo::archipelago arch;
    arch.set_topology(topology::unconnected());
    //std::vector<pagmo::population> decomposed_pop;

    //Each island in the archipelago solve a different single-objective problem
    for(pagmo::population::size_type i=0; i<pop.size();++i) {
        pagmo::problem::decompose prob(pop.problem());
        std::cout << "Vector " << i << ": "<< prob.get_weights() << std::endl;

        pagmo::population decomposed_pop(prob,0);

        //decomposed_pop.push_back(pagmo::population(prob,0));

        //Set the individual of the new population of the island
        for(pagmo::population::size_type j=0; j<pop.size();++j) {
            decomposed_pop.push_back(pop.get_individual(j).cur_x);
        }
        arch.push_back(pagmo::island(*m_solver,decomposed_pop));
    }
    arch.evolve(m_gen);
    arch.join();

    //The population is set to contain the best individual of each island
    for(pagmo::population::size_type i=0; i<pop.size();++i) {
        pop.set_x(i, arch.get_island(i)->get_population().champion().x);
        std::cout << "Champion " << i << ": "<< arch.get_island(i)->get_population().champion().x << std::endl;
    }
}

/// Algorithm name
std::string pade::get_name() const
{
    return "Parallel Decomposition (PaDe)";
}

/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string pade::human_readable_extra() const
{
	std::ostringstream s;
	s << "gen:" << m_gen << ' ';
    s << "solver:" << m_solver << ' ';
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::pade);
