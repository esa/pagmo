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

#include "aco.h"

namespace pagmo { namespace algorithm {

    /// Constructor.
    /**
     * @param[in] number of cycles. Each cycle, the m ants complete a tour of n cities
     */
     //TODO initialize members. Add getter functions. Add checks on cycle (!=0) and rho (0<=rho<1). 
     //     1-rho is the evaporation parameter. rho can be seen as the trail persistance
    aco::aco(int cycle, int ants, double rho):base() 
    {
    }

    /// Clone method.
    base_ptr aco::clone() const
    {
            return base_ptr(new aco(*this));
    }

    /// Evolve implementation.
    /**
     * Run the ACO algorithm for the number of generations specified in the constructors.
     *
     * @param[in,out] pop input/output pagmo::population to be evolved.
     */

    void aco::evolve(population &pop) const
    {
        // Let's store some useful variables.
        const problem::base_tsp &prob = dynamic_cast<const problem::base_tsp &>(pop.problem());
        const problem::base::size_type prob_i_dimension = prob.get_i_dimension();
        const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
        const population::size_type NP =  pop.size();

        //TODO: add check on the problem. The problem need to be an integer problem TSP like, single objective. The population size needs to be greater than 1 (at least 2 individuals are needed) 
        //      add the check that the number of ants need to be greater than the population size.

        //TODO: first task, we have to deposit the initial pheromone along the edges of the graph. The pheromone can be deposited or from the information given by the initial population, uniformly or random.
        //      Please prepare a protected function for this implementation. It needs to contain a switch between the three diffrent modes. 
        //      You can chose to represent the pheromones on the edges as a graph or as a matrix. I strongly suggest the second one. 
        //      You have to initialize a std::vector<std::vector<double> > containing the pheromone values. All these values must be positive  
        //      Uniform: initialize every element of the matrix (diagonal excluded) to a costant value c = M/C where M is the number of ants and C is the value is the length of the round-trip 
        //               obtained by applying the nearest-neighbour heuristic. This heuristic randomly selects a city and then creates a round-trip by always moving to the nearest city
        //      Random: initialize every element of the matrix (diagonal excluded) with random value between 0 and 1
        //      From population: use equation (1), (2), (3) from the paper enclosed in the email to update the elements of the matrix: tau_ij = rho*tau_ij + delta_tau_ij. Using
        //                       as tau_ij on the right hand side the constant c value of the uniform case. Delta_tau_ij is computed as in equation (2)-(3)


        // Main ACO loop: stopping condition either maximum number of cycle reached or all ants make the same tour
        for (int t = 0; t < m_cycle; ++t) {

        } // end of main ACO loop
    }



    /// Algorithm name
    std::string aco::get_name() const
    {
        return "Simple Ant System";
    }


    /// Extra human readable algorithm info.
    /**
     * Will return a formatted string displaying the parameters of the algorithm.
     */
    std::string aco::human_readable_extra() const
    {
        std::ostringstream s;
        s << "\nNumber of cycles: " << m_cycle
            << "\nNumber of ants: " << m_ants
            << "\nPheromone evaporation (Rho): " << m_rho << std::endl;
        return s.str();
    }

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::aco)