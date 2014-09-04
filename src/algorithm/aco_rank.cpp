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

#include "aco_rank.h"

namespace pagmo { namespace algorithm {
    
    /// Constructor.
    /**
     * @param[in] cycle - number of cycles. Each cycle, the m ants complete a tour of n cities.
     * @param[in] ants - number of ants in the colony.
     * @param[in] rho - each ant leaves a trail of pheromone which evaporates according to this constant.
     * @param[in] e - each elite ant leaves a stronger trail of pheromone. this constant influences the quantity.
     */
    aco_rank::aco_rank(int cycles, int ants, double rho, double alpha, double beta, double e): m_cycles(cycles), m_ants(ants), m_rho(rho), m_alpha(alpha), m_beta(beta), m_e(e)
    {
        init_checks();
    }
    
    /// Clone method.
    base_ptr aco_rank::clone() const
    {
        return base_ptr(new aco_rank(*this));
    }
    
    /// Evolve implementation.
    /**
     * Runs the ACO algorithm for the number of generations specified in the constructor.
     *
     * @param[in,out] pop input/output pagmo::population to be evolved.
     */
    void aco_rank::evolve(population &pop) const
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
            for (int ant = 0; ant < m_ants; ++ant)
                ant_tours.at(ant) = forage(dist(generator), weights, tau);
            
            // sort lowest to highest
            std::set<aco_tour> unique_sorted (ant_tours.begin(), ant_tours.end());
            double no_ranks = unique_sorted.size();
            
            // update and reinforce delta_tau according to rank (decreasing)
            for (int ant = 0; ant < m_ants; ++ant) {
                std::set<aco_tour>::iterator foundit = unique_sorted.find(ant_tours.at(ant).length);
                double weight = (no_ranks - std::distance(unique_sorted.begin(), foundit))/no_ranks;
                
                for (size_t i = 0; i < ant_tours.at(ant).tour.size() - 1; ++i) {
                    size_t from = ant_tours.at(ant).tour.at(i), to = ant_tours.at(ant).tour.at(i+1);
                    delta_tau.at(from).at(to) += (weight * m_e * Q/ant_tours.at(ant).length);
                }
                delta_tau.at( ant_tours.at(ant).tour.back() ).at( ant_tours.at(ant).tour.front() ) += ( weight * m_e * Q/ant_tours.at(ant).length );
            }
            
            /** reinforce delta_tau with elite ant -  (shortest tour) AGAIN
            aco_tour shortest_path = *unique_sorted.begin();
            for (size_t i = 0; i < shortest_path.tour.size() - 1; ++i) {
                size_t from = shortest_path.tour.at(i), to = shortest_path.tour.at(i+1);
                delta_tau.at(from).at(to) += (m_e * Q/shortest_path.length);
            }
            delta_tau.at(shortest_path.tour.back()).at(shortest_path.tour.front()) += (m_e * Q/shortest_path.length);
            */
            
            // update pheromone matrix
            for (size_t i = 0; i < tau.size(); ++i) 
                for (size_t j = 0; j < tau.size(); ++j)
                    tau.at(i).at(j) = tau.at(i).at(j) * m_rho + delta_tau.at(i).at(j);
            
            // make tour consistent
            make_tour_consistent(ant_tours.front().tour);
            //convert to decision vector and 
            // save the shortest path or the last n cycles in the population. 
            // where n is the number of individuals in the population. 
            pop.set_x( NP > 0 ? --NP : pop.size()-1 , tour2chromosome(ant_tours.front().tour));
            
            // store lambdas
            m_lambda.push_back( get_l_branching(0.5, tau) );
            
            // stop cycles if we're close to three and f'(x) -> 0
            if (t > m_cycles/4 && m_lambda.at(t) < 4 && abs(m_lambda.at(t) - m_lambda.at(t-2)) )
                break;
        } // end of main ACO loop (cycles)
        
    }

    /// Algorithm name
    std::string aco_rank::get_name() const
    {
        return "Rank-based Ant System";
    }

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::aco_rank)