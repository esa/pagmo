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

//#include <boost/graph/graph_concepts.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>

#include "../exceptions.h"
#include "../types.h"
#include "tsp.h"
#include <boost/graph/graphviz.hpp>

namespace pagmo { namespace problem {

    static vector2D<double> const default_weights = {
        {0, 1,              2,              3,              4},
        {1, 0,              2.236067,       4,              4.123105},
        {2, 2.236067,       0,              3.605551,       6},
        {3, 4,              3.605551,       0,              5},
        {4, 4.123105,       6,              5,              0}
    };
    /// Default constructor.
    /**
     * This constructor will build a TSP instance with 5 cities in the following positions:
     * - city 1: ( 0,    0)
     * - city 2: ( 1,    0)
     * - city 3: ( 0,    2)
     * - city 4: (-3,    0)
     * - city 5: ( 0,   -4)
     */
    tsp::tsp(): base_tsp(5) {        
        set_graph(default_weights);
        m_weights = default_weights;
        set_lb(0);
        set_ub(5); //number of nodes/vertices in the graph
    }
    
    /// Constructor from boost graph object
    /**
     * Initialize the number of vertices and the BGL tsp_graph
     * @param[in] boost graph object containing weights between vertices.
     */
    tsp::tsp(tsp_graph const& graph): base_tsp(boost::num_vertices(graph)) {
        set_graph(graph);
        base_tsp::convert_graph_to_vector2D(graph, m_weights);
        set_lb(0);
        set_ub(boost::num_vertices(graph)); //number of nodes in the graph
    }

    /// Constructor from vectors and maximum weight.
    /**
     * Initialize weights between vertices (city distances) from the matrix.
     * @param[in] weights vector of distances between cities
     */
    tsp::tsp(vector2D<double> const& weights): base_tsp(get_no_vertices(weights)) {
        set_graph(weights);
        m_weights = weights;
        set_lb(0);
        set_ub(boost::num_vertices(m_graph));
    }

    /// Clone method.
    base_ptr tsp::clone() const
    {
        return base_ptr(new tsp(*this));
    }

    /// Implementation of the objective function
    /**
     * Computes the fitness vector associated to a decision vector.
     * The fitness is defined as Sum_ij(w_ij * x_ij) 
     * where w_ij are the weights defining the distances between the cities
     * @param[out] f fitness vector
     * @param[in] x decision vector (a permutation of the vertices)
     */
    void tsp::objfun_impl(fitness_vector &f, decision_vector const& x) const {
        size_t dimension = get_dimension();
        
        pagmo_assert(f.size() == 1);
        pagmo_assert(x.size() == dimension && 
                boost::numeric_cast<int>(x.size()) == get_no_vertices(m_weights));
        
        f[0] = 0;
        for (size_t i = 1; i < dimension; ++i)
            f[0] += m_weights[boost::numeric_cast<int>(x[i-1])][boost::numeric_cast<int>(x[i])];
        f[0] += m_weights[boost::numeric_cast<int>(x[dimension-1])][boost::numeric_cast<int>(x[0])];        
    }

    /// Constraint computation.
    /**
     * The decision vector has to be a permutation of the set of nodes.
     * The constraint is positive when not satisfied:
     *  - if we have selected more than once the same node or 
     *  - equivalently not all the nodes have been selected.
     */
    void tsp::compute_constraints_impl(constraint_vector &c, decision_vector const& x) const {
        c[0] = 0; // assume all is well
        
        size_t dimension = get_dimension();
        
        // Checks if the length is the same
        if (x.size() != dimension ||  
                boost::numeric_cast<int>(x.size()) == get_no_vertices(m_weights)) {
            c[0] = 1;
            return;
        }
        
        // Checks if there are duplicate items
        decision_vector temp(x);
        std::sort(temp.begin(), temp.end());
        if ( std::unique(temp.begin(), temp.end()) != temp.end() ) {
            c[0] = 1;
            return;
        }
        
        // Checks if there actually is an edge between two adjacent vertices
        for (size_t i = 1; i < dimension; ++i) {
            if (!m_weights[boost::numeric_cast<int>(x[i-1])][boost::numeric_cast<int>(x[i])]) {
                c[0] = 1;
                return;
            }
        }
        if (!m_weights[boost::numeric_cast<int>(x[dimension-1])][boost::numeric_cast<int>(x[0])])
            c[0] = 1; // finish to start
    }

    /// Extra human readable info for the problem.
    /**
     * Will return a list of vertices and edges
     */
    std::string tsp::human_readable_extra() const
    {
        std::ostringstream oss;
        oss << "The Boost Graph (Adjacency List): \n";// << m_graph << std::endl;
//        boost::write_graphviz(oss, m_graph, boost::make_edge_attributes_writer( boost::get(boost::edge_weight_t(), m_graph) ) );
        
        tsp_vertex_map_const_index vtx_idx = boost::get(boost::vertex_index_t(), m_graph);
        tsp_edge_map_const_weight weights = boost::get(boost::edge_weight_t(), m_graph);

        oss << "Vertices = { ";
        
        tsp_vertex_range_t v_it;
        for (v_it = boost::vertices(m_graph); v_it.first != v_it.second; ++v_it.first)
                oss << vtx_idx[*v_it.first] <<  " ";
        oss << "}" << std::endl;
        
        oss << "Edges (Source, Target) = Weight : " << std::endl;
        
        tsp_edge_range_t e_it;
        for (e_it = boost::edges(m_graph); e_it.first != e_it.second; ++e_it.first) {
            int i = vtx_idx[boost::source(*e_it.first, m_graph)];
            int j = vtx_idx[boost::target(*e_it.first, m_graph)];
            oss << "(" << i << ", " << j<< ") = " << weights[*e_it.first] << std::endl;
        }
        oss << std::endl;
        
        return oss.str();
    }

    std::string tsp::get_name() const 
    {
        return "Traveling Salesman Problem";
    }

}} //namespaces

//BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::tsp);