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

#include <deque>

#include "aco.h"

// debug
#define USE_DEBUG

#ifdef USE_DEBUG
#define Debug( x ) std::cout << x
#else
#define Debug( x ) 
#endif

namespace pagmo { namespace algorithm {
    
    /// Constructor.
    /**
     * @param[in] cycle - number of cycles. Each cycle, the m ants complete a tour of n cities.
     * @param[in] ants - number of ants in the colony.
     * @param[in] rho - each ant leaves a trail of pheromone which evaporates according to this constant.
     */
    aco::aco(int cycle, int ants, double rho):base(), m_cycles(cycle), m_ants(ants), m_rho(rho)
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
    int aco::get_cycles() const
    {
        return m_cycles;
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
    
    // Sets the number of cycles the algorithm runs for.
    /**
     * Setter for the m_cycle property
     * @param cycles - int, the number of cycles
     */
    void aco::set_cycles(int cycles)
    {
        m_cycles = cycles;
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
    
    /// Enforces a rule for the starting position on the tours taken by the ants.
    /**
     * Selects the minimum index of a vertex and rotates elements around that 
     * item in order to produce consistency of the same trip, by always starting
     * from the vertex with the lowest index. e.g 3->4->1->2 becomes 1->2->3->4
     * Complexity is linear O(n)
     * 
     * @param[in/out] trip - a vector of vertices (e.g. list of cities)
     */
    void aco::make_tour_consistent(std::vector<size_t>& trip) 
    {
        std::rotate(trip.begin(), std::min_element(trip.begin(), trip.end()), trip.end());
    }
    
    /// Converts a tsp tour to a pagmo chromosome
    /**
     * Converts a tsp tour to a binary pagmo chromosome.
     * @param[in] trip a vector representing the sequence of vertices visited
     * @return decision vector - binary vector, see tsp::compute_idx
     */
    decision_vector aco::tour2chromosome(std::vector<size_t> trip) 
    {
        // make tour consistent (smallest element idx always first)
        make_tour_consistent(trip);
        
        // figure out the size of the matrix
        size_t n = *std::max_element(trip.begin(), trip.end()) + 1;
        // by default there are no connections
        decision_vector chromosome(n*(n-1), 0);
        
        // iterate through the vector and create chromosome
        for (size_t k = 0; k < trip.size() - 1; ++k) {
            chromosome.at( problem::tsp::compute_idx(trip.at(k), trip.at(k+1), n) ) = 1;
        }
        chromosome.at( problem::tsp::compute_idx(trip.back(), trip.front(), n) ) = 1;
        
        return chromosome;
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

    
    
    
    /// Evolve implementation.
    /**
     * Runs the ACO algorithm for the number of generations specified in the constructor.
     *
     * @param[in,out] pop input/output pagmo::population to be evolved.
     */
    void aco::evolve(population &pop) const
    {
        // Configuration parameters, maybe put them somewhere else?
        // Values taken from paper
        const double alpha = 1;
        const double beta = 5;
        const double c = 0.00001;
        const double Q = 100;
        
        // Let's store some useful variables.
        const problem::tsp &tsp_prob = dynamic_cast<const problem::tsp &>(pop.problem());
        const std::vector<std::vector<double> > &weights = tsp_prob.get_weights();
        const problem::base::size_type dim = tsp_prob.get_i_dimension();
        const decision_vector &lb = tsp_prob.get_lb(), &ub = tsp_prob.get_ub();
        const population::size_type NP = pop.size();
        const problem::base::size_type no_vertices = tsp_prob.get_n_vertices();
        
        // random engine for ant init
        std::default_random_engine generator(time(NULL));
        std::uniform_int_distribution<size_t> dist(0, no_vertices - 1);

        //TODO: add check on the problem. The problem needs to be an integer problem TSP like, single objective.
        if (NP <= 1) pagmo_throw(value_error, "population size must be greater than one.");
        if (m_ants <= static_cast<int>(NP)) pagmo_throw(value_error, "the number of ants needs to be greater than the population size.");

        // keep a counter for duplicate items, we might want to terminate early and/or reinitialize
        size_t duplicates = 0;
        
        // for each cycle, store a map of the lowest cost and the cities visited
        // we get solutions ordered for free, though maps have O(log(n)) insertion/delete and lookup time
        std::map<double, std::vector<size_t> > winners;
        
        // Main ACO loop: stops when either maximum number of cycles reached or (? all ants make the same tour )
        for (int t = 0; t < m_cycles; ++t) {
            // deposit the initial pheromone along the edges of the graph
            std::vector<std::vector<double> > tau(initialize_pheromone(no_vertices, c));
            
            // set delta_tau to 0
            std::vector<std::vector<double> > delta_tau(no_vertices, std::vector<double>(no_vertices, 0));
            
            // initialize the tour length to 0
            std::vector<double> tour_length(m_ants, 0);
            
            // set the starting point for each ant using a random uniform distribution
            std::vector<std::vector<size_t> > tabu(m_ants);
            for (int k = 0; k < m_ants ; ++k)
                tabu.at(k).push_back(dist(generator));
                       
            /**
             *                      tau_{i,j}^alpha * eta_{i,j}^beta            (1)
             * p_{i,j}^k = ----------------------------------------------------
             *              sum_{tabu_k} (tau_{i,k}^alpha * eta_{i,k}^beta)     (2)
             * 
             * where eta_{i,j} = 1/weights_{i,j}
             */
            // compute shortest paths and costs
            for (int ant = 0; ant < m_ants; ++ant) {
//                Debug("\n ant #" << ant << ": ");
                
                // get the starting position according to initialization
                size_t current = tabu.at(ant).at(0);
//                Debug(" start -> " << current);

                // keep a sum of probabilities for all ants (2)
                // initialize to smallest possible positive value to avoid division by zero
                double prob_sum = std::numeric_limits<double>::min();
                
                // keep moving until all vertices are visited
                while (tabu.at(ant).size() < no_vertices) {
                    // for each round, keep only maximum probability and it's position
                    double prob_max = 0; 
                    size_t next = std::numeric_limits<size_t>::max();
                    
                    // compute transition probabilities for all valid vertices
                    for (size_t possible = 0; possible < no_vertices; ++possible) {
                        // already visited, skip vertex
                        if (std::find(tabu.at(ant).begin(), tabu.at(ant).end(), possible) != tabu.at(ant).end()) continue;
                        // not a valid transition, skip vertex (checking weights for 0 and NaNs)
                        // this might be useless since problem::tsp::check_matrix does this already
                        if (weights.at(current).at(possible) == 0 || !weights.at(current).at(possible) == weights.at(current).at(possible)) continue;
                        
                        // otherwise compute probability (1)
                        double prob_next = ( pow(tau.at(current).at(possible), alpha) * pow(1/weights.at(current).at(possible), beta) ) / prob_sum;
                        // keep track of maximum and it's position
                        if (prob_max < prob_next) {
                            prob_max = prob_next;
                            next = possible;
                        }
                    } // done searching for next transition
                    
                    if (next == std::numeric_limits<size_t>::max()) { // we haven't found a possible next step for this ant
                        // stop searching for this ant, consider it stuck (set cost to infinity)
                        tour_length.at(ant) = std::numeric_limits<double>::max();
                        break; //use goto?
                    }
//                    Debug(next << " ");
                    
                    // remember visited vertices for each ant
                    tabu.at(ant).push_back(next);
                    
                    // remember the total cost for each ant
                    tour_length.at(ant) += weights.at(current).at(next);
                    
                    // update \delta \tau_{i,j} (this sums over all ants)
                    delta_tau.at(current).at(next) += Q/tour_length.at(ant);
                    
                    // remember sum of probabilities (2)
                    prob_sum += prob_max;
                    
                    // move to next vertex
                    current = next;
                } // done with tabu list, we only need to make the round trip
                
                // check to see if it's a valid transition
                size_t last = tabu.at(ant).back();
                size_t first = tabu.at(ant).front();
                if (weights.at(last).at(first) == 0 || !weights.at(last).at(first) == weights.at(last).at(first)) {
                    // set cost to infinity if we can't do a round-trip
                    tour_length.at(ant) = std::numeric_limits<double>::max();
                    continue;
                } else { // add the cost from last to first (init)
                    tour_length.at(ant) += weights.at(last).at(first);
                }
                
//                Debug("stop => cost: " << tour_length.at(ant));
            } // finished with all ants, a new cycle can start from here
            
            // update pheromone trail: \tau(t+1) = \tau(t) * rho + \delta \tau_{i,j}
            for (size_t i = 0; i < no_vertices; ++i)
                for(size_t j = 0; j < no_vertices; ++j)
                    tau.at(i).at(j) = tau.at(i).at(j) * m_rho + delta_tau.at(i).at(j);
            
            // reset \delta \tau_{i,j} to 0
            std::vector<std::vector<double> > temp(no_vertices, std::vector<double>(no_vertices, 0));
            delta_tau = temp;
            
            // store cost and path for the winner ant
            std::vector<double>::iterator low_it = std::min_element(tour_length.begin(), tour_length.end());
            double lowest_cost = *low_it;
            std::vector<size_t> shortest_path = tabu.at( std::distance(tour_length.begin(), low_it) );
            
            Debug(print_histogram(tour_length));               

            // optional step, linear takes O(n), rotate around minimum to have consistent paths
            make_tour_consistent(shortest_path);
            
            Debug("\n Finished cycle " << t << "\t Score: " << lowest_cost << " Tour: " << shortest_path);
            
            // if not true, it's not inserted, we have a duplicate key (hence cost, hence trip)
            if(!winners.insert(std::pair<double, std::vector<size_t> >(lowest_cost, shortest_path)).second) {
                ++duplicates;
            } else {
                // convert to decision vector and 
                // Save the shortest path (winning ant) for the last n cycles in the population. 
                // Where n is the number of individuals in the population. 
                pop.set_x(t, tour2chromosome(shortest_path));
            }
            
            // stop cycles when there are 10% duplicates
//            if (duplicates > static_cast<size_t>(m_ants * 0.1))
//                break;
        } // end of main ACO loop (cycles)
        
        Debug("\n # duplicates: " << duplicates << " and # winners: " << winners.size());
        Debug("\n Lowest Cost: " << winners.begin()->first << " Tour: " << winners.begin()->second);
        
        // do something with these?
        (void)lb;
        (void)ub;
        (void)dim;
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
        s << "\nNumber of cycles: " << m_cycles
            << "\nNumber of ants: " << m_ants
            << "\nPheromone evaporation (Rho): " << m_rho << std::endl;
        return s.str();
    }

    /// Prints a histogram of a fitness tour
    /**
     * Prints a histogram of the fitness tour.
     * Each line is a unique fitness (cost) value for a tour.
     * @param tour - the tour vector containing costs
     * @return a string of the output
     */
    std::string aco::print_histogram(std::vector<double> tour) const
    {
        std::stringstream out;
        out << "\n\nHistogram (no tours / cost bin)\n";
        std::sort(tour.begin(), tour.end());
        std::vector<double>::iterator it;
        it = std::unique(tour.begin(), tour.end());
        std::vector<double> uniq(tour.begin(), it);        
        std::vector<int> histo(uniq.size(), 0);
        
        for (size_t i = 0; i < uniq.size(); ++i) {
            for (size_t j = 0 ; j < tour.size(); ++j)
                if (uniq[i] == tour[j])
                    histo[i]++;
        }
        
        for (size_t k = 0; k < histo.size(); ++k)
            out << uniq[k] << " " <<  std::string(histo[k], '#') << std::endl;
        
        return out.str();
    }
    
}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::aco)