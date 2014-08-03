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
     * @param[in] cycle - number of cycles. Each cycle, the m ants complete a tour of n cities.
     * @param[in] ants - number of ants in the colony.
     * @param[in] rho - each ant leaves a trail of pheromone which evaporates according to this constant.
     */
    aco::aco(int cycle, int ants, double rho):base(), m_cycle(cycle), m_ants(ants), m_rho(rho)
    {
        if (cycle <= 0) pagmo_throw(value_error, "the number of cycles must be positive, non negative.");
        if (rho < 0 || rho >= 1) pagmo_throw(value_error, "the pheromone evaporation constant (rho) must be in [0, 1).");
        // if population must be at least 2 and ants > population => ants must be at least 3
        if (ants <= 2) pagmo_throw(value_error, "there must be at least three ants in the colony.");
    }
    
    /// Returns the number of cycles the algorithm is run
    /**
     * Returns the number of cycles the algorithm is run
     * @return int m_cycle
     */
    int aco::get_cycle() const
    {
        return m_cycle;
    }
    
    /// Returns the number of ants in the colony
    /**
     * Returns the number of ants in the colony
     * @return int m_ants
     */
    int aco::get_ants() const
    {
        return m_ants;
    }
    
    /// Returns Rho - the pheromone evaporation constant
    /**
     * Returns Rho - the pheromone evaporation constant
     * @return double m_rho
     */
    double aco::get_rho() const
    {
        return m_rho;
    }

    /// Clone method.
    base_ptr aco::clone() const
    {
            return base_ptr(new aco(*this));
    }
    
    /// Performs a greedy nearest neighbor traverse round-trip of a fully connected graph.
    /**
     * Returns a list of vertices (cities) in a graph, by always selecting 
     * the nearest neighboring city. We assume that the adjacency matrix 
     * is fully connected, e.g. there is a possible cycle starting from any city.
     * 
     * @param start - the initial position, has to also be the stop.
     * @param matrix - the adjacency matrix containing the graph.
     * @return a list of vertices according to the greedy traverse
     */
    std::vector<int> aco::nn_trip(int start, const std::vector<std::vector<double> > matrix) const
    {
        std::vector<int> visited;
        std::vector<double> current = matrix.at(start);
        while(visited.size() < matrix.size())
        {
            // get index of closest neighbor
            int idx = std::distance(current.begin(), std::min_element(current.begin(), current.end()));
            // push back in the list
            visited.push_back(idx);
            // select as current row in matrix, repeat
            current = matrix.at(idx);
        }
        return visited;
    }
    
    /// Initializes the pheromone levels randomly
    /**
     * Sets the pheromone levels of each edge from a random distribution (0,1)
     * @param[in] dimension - the number of vertices in the graph
     * @param[in] low - lower bound for random variable generation, default 0
     * @param[in] high - upper bound for random variable generation, default 1
     */
    std::vector<std::vector<double> > aco::initialize_pheromone(int dimension, double low=0, double high=1) const 
    {
        std::default_random_engine rengine(time(NULL)); // seed software PRNG
        std::uniform_real_distribution<double> distr(low, high); // range
     
        std::vector<std::vector<double> > pheromones(dimension, std::vector<double>(dimension, 0));
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                if (i==j) continue;
                pheromones[i][j] = distr(rengine);
            }
        }
        return pheromones;
    }
    
    /// Initializes the pheromone levels uniformly.
    /**
     * Sets the pheromones levels of the edges in the matrix as uniform.
     * Uniform: initialize every element of the matrix (diagonal excluded) 
     * to a constant value c = M/C where M is the number of ants and 
     * C is the length of the round-trip obtained by applying 
     * the nearest-neighbor heuristic which randomly selects a city and then 
     * creates a round-trip by always moving to the nearest city.
     * 
     * @param[in] dimension - the number of vertices in the graph
     * @param[in] matrix - the matrix containing the distances between vertices
     */
    std::vector<std::vector<double> > aco::initialize_pheromone(int dimension, const std::vector<std::vector<double> > matrix) const 
    {
        double C = 0;
        // make a tour of the graph, starting randomly using a greedy nn search
        std::vector<int> trip = nn_trip(rand() % dimension, matrix);
        // compute the cost of the round-trip
        for (int i = 0; i < (int)trip.size() - 1; ++i)
            C += matrix.at(i).at(i+1);
        C += matrix.at(trip.size()).at(0);
        
        // initialize pheromone matrix according to formula M/C
        std::vector<std::vector<double> > pheromones(dimension, std::vector<double>(dimension, 0));
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                if(i==j) continue;
                pheromones[i][j] = m_ants / C;
            }
        }
        return pheromones;
    }
    
    /// Initializes the pheromone levels in the matrix according to the population.
    /**
     * Sets the pheromone levels in the matrix according to the population.
     * 
     * From population: use equation (1), (2), (3) from the paper enclosed 
     * in the email to update the elements of the matrix: 
     *  tau_ij = rho * tau_ij + delta_tau_ij. 
     *  Using as tau_ij on the right hand side, the constant c value of the uniform case. 
     * Delta_tau_ij is computed as in equation (2)-(3)
     * 
     * @param[in] dimension - the number of vertices in the graph
     * @param[in] population
     */
    std::vector<std::vector<double> > aco::initialize_pheromone(int dimension, const population& pop) const 
    {
        std::vector<std::vector<double> > pheromones(dimension, std::vector<double>(dimension, 0));
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                if (i==j) continue;
                pheromones[i][j] = 1;
            }
        }
        return pheromones;
    }

    /// Evolve implementation.
    /**
     * Runs the ACO algorithm for the number of generations specified in the constructor.
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

        //TODO: add check on the problem. The problem need to be an integer problem TSP like, single objective.
        if (NP <= 1) pagmo_throw(value_error, "population size needs to be greater than one.");
        if (m_ants <= NP) pagmo_throw(value_error, "the number of ants needs to be greater than the population size.");

        //TODO: first task, we have to deposit the initial pheromone along the edges of the graph.

        // Main ACO loop: stopping condition either maximum number of cycle reached or all ants make the same tour
        for (int t = 0; t < m_cycle; ++t) {
            break;
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