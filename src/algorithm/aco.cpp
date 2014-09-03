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
    aco::aco(int cycle, int ants, double rho):base(), m_cycles(cycle), m_ants(ants), m_rho(rho)
    {
        if (cycle <= 0) pagmo_throw(value_error, "the number of cycles must be positive, non negative.");
        if (rho < 0 || rho >= 1) pagmo_throw(value_error, "the pheromone evaporation constant (rho) must be in [0, 1).");
        // if population must be at least 2 and ants > population => ants must be at least 3
        if (ants <= 2) pagmo_throw(value_error, "there must be at least three ants in the colony.");
        m_alpha = 1;
        m_beta = 2;
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
    
    /// Returns the vector of lambdas for each cycle.
    std::vector<double> aco::get_lambda() const
    {
        return m_lambda;
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
        const double c = 0.1;
        const double Q = 1;
        
        // Let's store some useful variables.
        const problem::tsp &tsp_prob = dynamic_cast<const problem::tsp &>(pop.problem());
        const std::vector<std::vector<double> > &weights = tsp_prob.get_weights();
        const problem::base::size_type no_vertices = tsp_prob.get_n_vertices();
        size_t NP = pop.size();
        
        // random engine for ant init
        std::default_random_engine generator(time(NULL));
        std::uniform_int_distribution<size_t> dist(0, no_vertices - 1);

        //TODO: add check on the problem. The problem needs to be an integer problem TSP like, single objective.
        if (NP <= 1) pagmo_throw(value_error, "population size must be greater than one.");
        if (m_ants <= static_cast<int>(NP)) pagmo_throw(value_error, "the number of ants needs to be greater than the population size.");
        
        // deposit the initial pheromone along the edges of the graph
        std::vector<std::vector<double> > tau(initialize_pheromone(no_vertices, c));
        
        // Main ACO loop: stops when either maximum number of cycles reached or all ants make the same tour
        for (int t = 0; t < m_cycles; ++t) {    
            // (re)set delta_tau to 0
            std::vector<std::vector<double> > delta_tau(no_vertices, std::vector<double>(no_vertices, 0));
            
            // keep all the ant tours for one cycle here
            std::vector<aco_tour> ant_tours(m_ants);

            // compute shortest paths and costs for all ants
            for (int ant = 0; ant < m_ants; ++ant) {                
                // launch an ant and make a round-trip from a random starting position
                ant_tours.at(ant) = forage(dist(generator), weights, tau);
                
                // update delta_tau
                for (size_t i = 0; i < ant_tours.at(ant).tour.size() - 1; ++i) {
                    size_t from = ant_tours.at(ant).tour.at(i), to = ant_tours.at(ant).tour.at(i+1);
                    delta_tau.at(from).at(to) += Q/ant_tours.at(ant).length;
                }
                delta_tau.at( ant_tours.at(ant).tour.back() ).at( ant_tours.at(ant).tour.front() ) += Q/ant_tours.at(ant).length;
                
            } // finished with all ants, cycle computations can be performed
            
            // search for lowest cost and get corresponding tour after this cycle
            std::vector<aco_tour>::iterator low_it = std::min_element(ant_tours.begin(), ant_tours.end());
            std::vector<size_t> shortest_path = ant_tours.at( std::distance(ant_tours.begin(), low_it) ).tour;
            
            // update pheromone matrix
            for (size_t i = 0; i < tau.size(); ++i) 
                for (size_t j = 0; j < tau.size(); ++j)
                    tau.at(i).at(j) = tau.at(i).at(j) * m_rho + delta_tau.at(i).at(j);
            
            // make tour consistent
            make_tour_consistent(shortest_path);
            //convert to decision vector and 
            // save the shortest path or the last n cycles in the population. 
            // where n is the number of individuals in the population. 
            pop.set_x( NP > 0 ? --NP : pop.size()-1 , tour2chromosome(shortest_path));
            
            // store lambdas
            m_lambda.push_back( get_l_branching(0.5, tau) );
            
            // stop cycles if we're close to three and f'(x) -> 0
            if (t > m_cycles/4 && m_lambda.at(t) < 4 && abs(m_lambda.at(t) - m_lambda.at(t-2)) )
                break;
            
//            system("clear");
//            std::cout << print_histogram(ant_tours);
//            std::cout << print_tau(tau);
//            std::cout << "\n Finished cycle " << t+1 << "\t Score: " << (*low_it).length << " Tour: " << shortest_path;
        } // end of main ACO loop (cycles)
        
    }
    
    /// Computes the average lambda branching factor
    /**
     * Computes and returns the average lambda branching factor
     * 
     * avg lambda = sum( lambda_ij )
     * lambda_ij = sum ( tau_ij >= min(tau_i) + lambda ( max(tau_i) - min(tau_i) ) )
     * 
     * @param[in] tau - the pheromone matrix
     * @return - the average lambda branching factor
     */
    double aco::get_l_branching(double lambda, const std::vector<std::vector<double> >& tau) const
    {
        double avg_lambda = 0;
        size_t no_vertices = tau.size();
        for (size_t i = 0; i < no_vertices; ++i) {
            auto min_max = std::minmax_element(tau.at(i).begin(), tau.at(i).end());
            double min = *min_max.first, max = *min_max.second;
            
            avg_lambda += std::count_if(tau.at(i).begin(), tau.at(i).end(), std::bind1st(std::less<double>(), min + lambda * (max - min)));
        }
        return avg_lambda / static_cast<double>(no_vertices);
    }
    
    /// Performs a round-trip from a given starting vertex
    /**
     * Makes one round-trip given a starting position and a weight matrix
     * 
     *                      tau_{i,j}^alpha * eta_{i,j}^beta            (1)
     * p_{i,j}^k = ----------------------------------------------------
     *              sum_{allowed k} (tau_{i,k}^alpha * eta_{i,k}^beta)  (2)
     * 
     * where eta_{i,j} = 1/weights_{i,j}
     * 
     * @param[in] start - the starting vertice
     * @param[in] weights - the weights matrix
     * @return aco_tour
     */
    aco_tour aco::forage(const size_t start, const std::vector<std::vector<double> >& weights, const std::vector<std::vector<double> >& tau) const
    {       
        // init tour with starting position
        std::vector<size_t> tour(1, start);
        // init length to max (infinity)
        double tour_length = 0;
        
        // get number of vertices
        size_t no_vertices = weights.size();
        
        size_t current = tour.front();
        // keep moving until all vertices are visited
        while (tour.size() < no_vertices) {
            // for each round, keep only maximum probability and it's position
            double prob_max = 0; 
            size_t next = std::numeric_limits<size_t>::max();

            // compute transition probabilities for all valid vertices
            for (size_t possible = 0; possible < no_vertices; ++possible) {
                // already visited, skip vertex
                if ( std::find(tour.begin(), tour.end(), possible) != tour.end() ) continue;
                // not a valid transition, skip vertex (checking weights for 0)
                if ( !weights.at(current).at(possible) ) continue;

                // otherwise compute probability, denominator in formula is irrelevant
                double prob_next = pow(tau.at(current).at(possible), m_alpha) * pow(1/weights.at(current).at(possible), m_beta);

                // keep track of maximum and it's position
                if (prob_max < prob_next) {
                    prob_max = prob_next;
                    next = possible;
                }
            } // done searching for next transition

            // if we haven't found a possible next step for this ant, return max (infinity)
            if ( next == std::numeric_limits<size_t>::max() )
                return aco_tour(std::numeric_limits<double>::max(), tour);

            // remember visited vertices for each ant
            tour.push_back(next);

            // remember the total cost for each ant
            tour_length += weights.at(current).at(next);

            // move to next vertex
            current = next;
        } // done with tour list
        
        // compute the round trip cost, if round-trip possible
        if ( weights.at(tour.back()).at(tour.front()) != 0 )         
            tour_length += weights.at(tour.back()).at(tour.front());
        else
            return aco_tour(std::numeric_limits<double>::max(), tour);
        
        // return cost and tour
        return aco_tour(tour_length, tour);
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
    
    /// Prints a pheromone matrix
    std::string aco::print_tau(const std::vector<std::vector<double> >& tau) const 
    {
        std::stringstream out;
        out << "\n";
        for (size_t i = 0; i < tau.size(); ++i) {
            out << "\n";
            for(size_t j = 0; j < tau.at(i).size(); ++j)
                out << std::scientific << std::setprecision(1) << tau.at(i).at(j) << " ";
        }
        out << "\n";
        return out.str();
    }

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::aco)