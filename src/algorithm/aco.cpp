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
    
    /// Returns the number of cycles the algorithm is run.
    /**
     * Returns the number of cycles the algorithm is run
     * @return int m_cycle
     */
    int aco::get_cycle() const
    {
        return m_cycle;
    }
    
    /// Returns the number of ants in the colony.
    /**
     * Returns the number of ants in the colony
     * @return int m_ants
     */
    int aco::get_ants() const
    {
        return m_ants;
    }
    
    /// Returns Rho - the pheromone evaporation constant.
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
    
    /// Performs a greedy nearest neighbor traverse of a fully connected graph.
    /**
     * Returns a list of vertices (cities) in a graph, by always selecting 
     * the nearest neighboring city. We assume that the adjacency matrix 
     * is fully connected, e.g. there is a possible cycle starting from any city.
     * 
     * @param[in] start - the initial position, starting vertex
     * @param[in] matrix - the adjacency matrix of the graph
     * @return a list of vertices according to the greedy traverse
     */
    std::vector<size_t> aco::greedy_nn_trip(size_t start, const std::vector<std::vector<double> > matrix) const
    {
        size_t novertices = matrix.size();
        std::vector<size_t> visited;
        std::vector<double> current = matrix.at(start);
        while(visited.size() < novertices)
        {
            // get index of closest neighbor
            size_t idx = std::distance(current.begin(), std::min_element(current.begin(), current.end()));
            // mark as visited
            visited.push_back(idx);
            // go to next vertex, repeat
            current = matrix.at(idx);
        }
        return visited;
    }
    
    std::vector<bool> aco::list2decision_vector(const std::vector<size_t> trip) 
    {
        size_t n = *std::max_element(trip.begin(), trip.end());
        // by default there are no connections
        std::vector<bool> decision_vector(n*(n-1), 0);
        
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if (i==j) continue;
//                decision_vector.at( problem::tsp.compute_idx(i, j, n) ) = true;
            }
        }
        return decision_vector;
    }
    
    /// Initializes the pheromone levels randomly.
    /**
     * Sets the pheromone levels of each edge from a random distribution (0,1)
     * @param[in] dimension - the number of vertices in the graph
     * @param[in] low - lower bound for random variable generation, default 0
     * @param[in] high - upper bound for random variable generation, default 1
     */
    std::vector<std::vector<double> > aco::initialize_pheromone(size_t dim) const
    {
        std::default_random_engine rengine(time(NULL));
        std::uniform_real_distribution<double> distr(std::numeric_limits<double>::min(), std::numeric_limits<double>::max());
        std::vector<std::vector<double> > pheromones (dim, std::vector<double>(dim, 0));
        for (size_t i = 0; i < dim; ++i) {
            for (size_t j = 0; j < dim; ++j) {
                if (i==j) continue;
                pheromones.at(i).at(j) = distr(rengine);
            }
        }
        return pheromones;
    }
    
        /// Initializes the pheromone levels uniformly (absolute).
    /**
     * Sets the pheromones levels of the edges in the matrix as uniform
     * to a small constant c;
     * 
     * @param[in] dimension - the number of vertices in the graph
     * @param[in] c - a small constant 0 < c < 1
     */
    std::vector<std::vector<double> > aco::initialize_pheromone(size_t dim, double c) const
    {   
        if(c >= 1 || c <=0 ) pagmo_throw(value_error, "the initial pheromone level must be in (0, 1).");
        std::vector<std::vector<double> > pheromones (dim, std::vector<double>(dim, 0));
        for (size_t i = 0; i < dim; ++i) {
            for (size_t j = 0; j < dim; ++j) {
                if(i==j) continue;
                pheromones.at(i).at(j) = c;
            }
        }
        return pheromones;
    }
    
    /// Initializes the pheromone levels, according to a greedy traverse from each vertex (uniform).
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
    std::vector<std::vector<double> > aco::initialize_pheromone(size_t dim, const std::vector<std::vector<double> > matrix) const
    {   
        std::vector<std::vector<double> > pheromones (dim, std::vector<double>(dim, 0));
        // traverse starting from all vertices
        for (size_t v = 0; v < dim; ++v) {
            // perform a greedy traverse starting from this vertex
            std::vector<size_t> trip = greedy_nn_trip(v, matrix);
            // compute the cost of a one-way trip
            double cost = 0;
            for (size_t k = 0; k < trip.size() - 1; ++k)
                cost += matrix.at(k).at(k+1);
            // deposit pheromone according to #ants / round-trip distance
            for (size_t i = 0; i < trip.size() - 1; ++i)
                pheromones.at(i).at(i+1) = m_ants / cost*2;
        }
        return pheromones;
    }
    
    /// Initializes the pheromone levels in the matrix according to the population.
    /**
     * Sets the pheromone levels in the matrix according to the population.
     *  tau_ij = rho * tau_ij + delta_tau_ij
     * 
     * @param[in] dimension - the number of vertices in the graph
     * @param[in] population
     */
//    void aco::initialize_pheromone(int dimension, const population& pop) 
//    {
//        for (int i = 0; i < dimension; ++i) {
//            for (int j = 0; j < dimension; ++j) {
//                if (i==j) continue;
//                // compute distance to nearest neighbor
////                double C = *std::min_element(x.begin(), x.end());
//                // initialize pheromones
//                m_pheromone[i][j] = m_rho;
//            }
//        }
//    }

    
    
    /// Enforces a rule for the starting position on the tours taken by the ants.
    /**
     * Selects the minimum index of a vertex and rotates elements around that 
     * item in order to produce consistency of the same trip, by always starting
     * from the vertex with the lowest index. e.g 3->4->1->2 becomes 1->2->3->4
     * Complexity is linear O(n)
     * 
     * @param[in/out] trip - a vector of vertices (e.g. list of cities)
     */
    void aco::make_tour_consistent(std::vector<size_t>& trip) const 
    {
        std::rotate(trip.begin(), std::min_element(trip.begin(), trip.end()), trip.end());
    }

    /// Evolve implementation.
    /**
     * Runs the ACO algorithm for the number of generations specified in the constructor.
     *
     * @param[in,out] pop input/output pagmo::population to be evolved.
     */
    void aco::evolve(population &pop) const
    {
        // Configuration parameters, maybe put somewhere else? make members, set params for constructor?
        const double alpha = 1;
        const double beta = 5;
        const double c = 0.00001;
        const double Q = 100;
        
        // Let's store some useful variables.
        const problem::tsp &tsp_prob = dynamic_cast<const problem::tsp &>(pop.problem());
//        const std::vector<std::vector<double> > weights (tsp_prob.get_weights());
        const std::vector<std::vector<double> > &weights = tsp_prob.get_weights();
        const problem::base::size_type no_vertices = tsp_prob.get_i_dimension();
        const decision_vector &lb = tsp_prob.get_lb(), &ub = tsp_prob.get_ub();
        const population::size_type NP = pop.size();

        //TODO: add check on the problem. The problem needs to be an integer problem TSP like, single objective.
        if (NP <= 1) pagmo_throw(value_error, "population size needs to be greater than one.");
        if (m_ants <= static_cast<int>(NP)) pagmo_throw(value_error, "the number of ants needs to be greater than the population size.");

        // keep a counter for duplicate items, we might want to terminate early and/or reinitialize
        size_t duplicates = 0;
        
        // for each cycle, store a map of the lowest cost and the cities visited
        // we get solutions ordered for free, though maps have O(log(n)) insertion/delete and lookup time
        std::map<double, std::vector<size_t> > winners;
        
        // Main ACO loop: stops when either maximum number of cycles reached or (? all ants make the same tour )
        for (int t = 0; t < m_cycle; ++t) {
            // deposit the initial pheromone along the edges of the graph
            std::vector<std::vector<double> > tau(initialize_pheromone(no_vertices, c));
            // set delta_tau to 0
            std::vector<std::vector<double> > delta_tau(no_vertices, std::vector<double>(no_vertices, 0));
            
            // initialize the cost vector to 0
            std::vector<double> cost(m_ants, 0);
            
            // initialize the starting point using a random uniform distribution
            std::default_random_engine generator(time(NULL));
            std::uniform_int_distribution<size_t> dist(0, no_vertices);
            std::vector<std::vector<size_t> > tabu(m_ants);
            for (int k = 0; k < m_ants; ++k)
                tabu.at(k).at(0) = dist(generator);

            /**
             *                      tau_{i,j}^alpha * eta_{i,j}^beta            (1)
             * p_{i,j}^k = ----------------------------------------------------
             *              sum_{tabu_k} (tau_{i,k}^alpha * eta_{i,k}^beta)     (2)
             */
            // compute shortest paths and costs
            for (int ant = 0; ant < m_ants; ++ant) {
                // keep a sum of probabilities for all ants (2)
                // initialize to smallest possible positive value to avoid division by zero
                double prob = std::numeric_limits<double>::min();
                
                // get the starting position according to initialization
                size_t current = tabu.at(ant).at(0);
                
                // keep moving until all vertices are visited
                while (tabu.at(ant).size() < no_vertices) {
                    // for each round, keep only maximum probability and it's position
                    double max = 0; 
                    size_t next = 0;
                    
                    // compute transition probabilities for all valid vertices
                    for (size_t possible = 0; possible < no_vertices; ++possible) {
                        // already visited, skip vertex
                        if (std::find(tabu.at(ant).begin(), tabu.at(ant).end(), possible) != tabu.at(ant).end()) continue;
                        // not a valid transition, skip vertex (checking weights for 0 and NaNs)
                        // this might be useless since problem::tsp.check_matrix does this already
                        if (weights.at(current).at(possible) == 0 || !weights.at(current).at(possible) == weights.at(current).at(possible)) continue;
                        
                        // otherwise compute probability
                        double prob_next = ( pow(tau.at(current).at(possible), alpha) * pow(1/weights.at(current).at(possible), beta) ) / prob;
                        // keep track of maximum and it's position
                        if (max < prob_next)  {
                            max = prob_next;
                            next = possible;
                        }
                    } // done searching for next transition
                    
                    if (!next) { // we haven't found a possible next step for this ant
                        // option 1. kill the ant (can't do it since m_ants is a member and this method is const)
                        // option 2. reinitialize ant (though others might already be several steps ahead)
                        // option 3. do something smart with this path (mark it or something)
                        // option 4. stop searching paths for this ant, consider it stuck:
                        break; //use goto?
                    }
                    
                    // remember visited vertices for each ant
                    tabu.at(ant).push_back(next);
                    // remember the total cost for each ant
                    cost.at(ant) += weights.at(current).at(next)/Q;
                    // remember sum of probabilities (2)
                    prob += max;
                    
                    // update \delta \tau_{i,j} (this sums over all ants)
                    delta_tau.at(current).at(next) += cost.at(ant)/Q;
                    
                    // move to next vertex
                    current = next;
                } // done with tabu list, we should now have a complete trip
                
                // update \delta \tau_{i,j} (this sums over all ants)
//                for (size_t i = 0; i < no_vertices - 1; ++i)
//                    delta_tau.at( tabu.at(ant).at(i) ).at( tabu.at(ant).at(i+1) ) += cost.at(ant)/Q;

            } // finished with all ants, a new cycle can start from here
            
            // update pheromone trail: \tau(t+1) = \tau(t) * rho + \delta \tau_{i,j}
            for (size_t i = 0; i < no_vertices; ++i)
                for(size_t j = 0; j < no_vertices; ++j)
                    tau.at(i).at(j) = tau.at(i).at(j) * m_rho + delta_tau.at(i).at(j);
            
            // reset \delta \tau_{i,j} to 0
            std::vector<std::vector<double> > temp(no_vertices, std::vector<double>(no_vertices, 0));
            delta_tau = temp;
            
            // store cost and path for the winner ant
            std::vector<double>::iterator low_it = std::min_element(cost.begin(), cost.end());
            double lowest_cost = *low_it;
            std::vector<size_t> shortest_path = tabu.at( std::distance(cost.begin(), low_it) );
            
            // optional step, linear takes O(n), rotate around minimum to have consistent paths
            make_tour_consistent(shortest_path);
            
            // if not true, it's not inserted, we have a duplicate key (hence cost, hence trip)
            if(!winners.insert(std::make_pair(lowest_cost, shortest_path)).second)
                ++duplicates;
            
        } // end of main ACO loop (cycles)
        
        // do stuff with population & winners here
        // convert to decision vector
        // convert costs according to lower and upper bounds
        (void)lb;
        (void)ub;
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