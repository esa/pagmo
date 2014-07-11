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
    tsp::tsp(vector2D<double> const& weights): base_tsp(base_tsp::count_vertices(weights)) {
        set_graph(weights);
        m_weights = weights;
        set_lb(0);
        set_ub(boost::num_vertices(get_graph()));
    }

    /// Clone method.
    base_ptr tsp::clone() const
    {
        return base_ptr(new tsp(*this));
    }
    
    /**
     * Converts a binary decision vector to a binary adjacency matrix
     * @param[in] x decision vector
     * @param[out] mat two dimensional boolean vector
     */
//    void tsp::convert_decision_vector_to_vector2D(decision_vector const& x, vector2D<bool>& mat) {
//        size_t n = x.size();
//        n = n / (n + 1); // number of vertices
//        
//        size_t k = 0;
//        for (size_t i = 0; i < n; ++i) {
//            for (size_t j = 0; j < n; ++j) {
//                if(i==j) { ++k; continue; }
//
//                size_t index = i*n + j - k;
//                mat[i][j] = x[index];
//            }
//        }
//    }
    
    /**
     * Converts a boolean adjacency matrix to decision vector (chromosome)
     * @param[in] mat two dimensional boolean vector
     * @param[out] x decision vector
     */
//    void tsp::convert_vector2D_to_decision_vector(vector2D<bool> const& mat, decision_vector& x) {
//        size_t n = 5;//base_tsp::count_vertices<bool>(mat);
//        size_t k = 0;
//        for (size_t i = 0; i < n ; ++i) {
//            for (size_t j = 0 ; j < n; ++j) {
//                if (i == j) continue; // we don't add diagonal elements
//                x[k++] = mat[i][j];
//            }
//        }
//    }

    /// Implementation of the objective function
    /**
     * Computes the fitness vector associated to a decision vector.
     * The fitness is defined as Sum_ij(w_ij * x_ij) 
     * where w_ij are the weights defining the distances between the cities
     * The decision vector x_ij is a concatenated binary adjacency matrix, 
     * with the diagonal elements skipped since they're always zero because
     * in a route you can't go from one vertex to itself.
     * @param[out] f fitness vector
     * @param[in] x decision vector
     */
    void tsp::objfun_impl(fitness_vector &f, decision_vector const& x) const {
        size_t n = get_no_vertices();
        // dim = n^2 - n (matrix - diagonal)
        
        pagmo_assert(f.size() == 1);
        pagmo_assert(x.size() == get_dimension() && x.size() == n * (n - 1));
        
        f[0]= 0;
        size_t k = 0; 
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                if(i==j) { ++k; continue; }
                
                size_t index = i * n + j - k;
                f[0] += m_weights[i][j] * x[index];
            }
        }       
    }

    /// Constraint computation.
    /**
     * The sum for rows and columns of the binary adjacency matrix
     * has to be 0 for all rows and columns to be valid
     * @param[out] c constraint_vector
     * @param[in] x decision_vector
     */
    void tsp::compute_constraints_impl(constraint_vector &c, decision_vector const& x) const {
        size_t n = get_no_vertices();
        
        size_t row = 0;
        size_t k = 0;
        for (size_t i = 0; i < n-1; i++) {
            for (size_t j = 0; j < n*2-1; j += (n-1) ) {
                row++;

                // row sums
                c[k] += x[j / (n-1) + i * n];
                if ( row % (n-1) == 0 ) ++k;

                // column sums
                c[i+n] += x[i+j];
            }
        }
        // pagmo needs c(x)-1
        for (size_t i = 0; i < 2*n-1; ++i)
            c[i] -= 1; // figure out how to do this in previous loops
        
        // total size of c is n*2 - 1 (minus diagonal)
    }

    /// Extra human readable info for the problem.
    /**
     * Will return a std::string containing a list of vertices and edges
     */
    std::string tsp::human_readable_extra() const {
        tsp_graph const the_graph = get_graph();
        
        std::ostringstream oss;
        oss << "The Boost Graph (Adjacency List): \n";// << the_graph << std::endl;
//        boost::write_graphviz(oss, m_graph, boost::make_edge_attributes_writer( boost::get(boost::edge_weight_t(), m_graph) ) );
        
        tsp_vertex_map_const_index vtx_idx = boost::get(boost::vertex_index_t(), the_graph);
        tsp_edge_map_const_weight weights = boost::get(boost::edge_weight_t(), the_graph);

        oss << "Vertices = { ";
        
        tsp_vertex_range_t v_it;
        for (v_it = boost::vertices(the_graph); v_it.first != v_it.second; ++v_it.first)
                oss << vtx_idx[*v_it.first] <<  " ";
        oss << "}" << std::endl;
        
        oss << "Edges (Source, Target) = Weight : " << std::endl;
        
        tsp_edge_range_t e_it = boost::edges(the_graph);
        for (e_it = boost::edges(the_graph); e_it.first != e_it.second; ++e_it.first) {
            int i = vtx_idx[boost::source(*e_it.first, the_graph)];
            int j = vtx_idx[boost::target(*e_it.first, the_graph)];
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