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

    /**
     * The default constructor, calls default for base_tsp
     */
    tsp::tsp(): base_tsp(), m_weights(graph2matrix(get_graph()))
    {
        check_matrix(m_weights);
    }

    /**
     * Constructor from adjacency matrix weights
     * @param[in] weights
     */
    tsp::tsp(const std::vector<std::vector<double> >& weights): base_tsp(matrix2graph(weights)), m_weights(weights)
    {
        check_matrix(m_weights);
    }

    /**
     * Constructor from boost graph object
     * @param[in] graph - tsp_graph
     */
    tsp::tsp(const tsp_graph& graph): base_tsp(graph), m_weights(graph2matrix(get_graph())) 
    {
        check_matrix(m_weights);
    }

    /// Clone method.
    base_ptr tsp::clone() const
    {
        return base_ptr(new tsp(*this));
    }

    /// Returns a reference to the adjacency matrix
    /**
     * After constructing the tsp object, the weights matrix is initialized.
     * @return reference to the adjacency matrix
     */
    const std::vector<std::vector<double> >& tsp::get_weights() const 
    {
        return m_weights;
    }
    
    /// Converts a matrix to a boost graph
    /**
     * Converts an adjacency matrix (square two dimensional) to a boost graph
     * object of type tsp_graph, setting the internal property for the weights
     * to the values defined in the matrix.
     * @param matrix std::vector<std::vector<double>> square matrix
     * @return tsp_graph boost graph object with weights initialized from matrix
     */
    base_tsp::tsp_graph tsp::matrix2graph(const std::vector<std::vector<double> > matrix) 
    {
        tsp_graph retval;
        tsp_edge_map_weight weights = boost::get(boost::edge_weight_t(), retval);
        tsp_vertex from, to;
        tsp_edge link;
        
        // add vertices first 
        /* Checking if a vertex exists with no vertices inserted causes segfault
         * so we have to iterate 1st to get total number of vertices
         * then iterate again to insert them ... bummer
         * Couldn't figure out how to do it all in 2 for loops
         */
        size_t no_vertices = matrix.size();
        for (size_t v = 0; v < no_vertices; ++v)
            boost::add_vertex(v, retval);
        
        // add edges and weights
        for (size_t i = 0; i < no_vertices; ++i) {

            /* uncomment this and it's segfault
             * don't do this check and the logic is wrong
             */ 
            from = boost::vertex(i, retval);
//                if (from == tsp_graph::null_vertex())
//                    from = boost::add_vertex(i, the_graph);
            
            for (size_t j = 0 ; j < no_vertices; ++j) {
                // we don't allow connections to self
                if(i == j) continue;
                
                to = boost::vertex(j, retval);
                // create destination vertex only if not existent
                // for some reason this works, but is not enough
//                    if (to == tsp_graph::null_vertex())
//                        to = boost::add_vertex(j, the_graph);
                
                // create an edge connecting those two vertices
                link = (boost::add_edge(from, to, retval)).first;
                // add weight property to the edge
                weights[link] = matrix.at(i).at(j);
            }
        }
        return retval;
    }
    
    /// Converts a boost graph to a matrix
    /**
     * Converts a boost graph to a matrix (two dimensional std::vector of doubles).
     * The returned matrix is square and has the diagonal elements zeroed out.
     * @param graph boost graph of type tsp_graph
     * @return a two dimensional std::vector of doubles
     */
    std::vector<std::vector<double> > tsp::graph2matrix(const tsp_graph graph) 
    {
        size_t n_vertex = boost::num_vertices(graph);
        tsp_vertex_map_const_index vtx_idx = boost::get(boost::vertex_index_t(), graph);
        tsp_edge_map_const_weight weights = boost::get(boost::edge_weight_t(), graph);
        tsp_edge_range_t e_it = boost::edges(graph);
        
        std::vector<std::vector<double> > retval(n_vertex, std::vector<double>(n_vertex, 0.0));
        for (e_it = boost::edges(graph); e_it.first != e_it.second; ++e_it.first) {
            int i = vtx_idx[boost::source(*e_it.first, graph)];
            int j = vtx_idx[boost::target(*e_it.first, graph)];
            retval[i][i] = 0;
            retval[i][j] = weights[*e_it.first];
        }
        return retval;
    }

    /// Computes the index 
    /** 
     * Returns the index in the decision vector corresponding to the concatenated 
     * rows of a square matrix, which contains the edges between vertices.
     * @param[in] i row
     * @param[in] j column
     */
    size_t tsp::compute_idx(const size_t i, const size_t j, const size_t n) 
    {
        pagmo_assert( i!=j && i<n && j<n );
        return i*(n-1) + j - (j>i? 1:0);
    };
    
    /// Checks if we can instantiate a TSP or ATSP problem
    /**
     * Checks if a matrix (std::vector<std::vector<double>>) 
     * is square or bidirectional (e.g. no one way links between vertices).
     * If none of the two conditions are true, we can not have a tsp problem.
     * @param matrix - the adjacency matrix (two dimensional std::vector)
     * @throws pagmo_throw - matrix is not square and/or graph is not bidirectional
     */
    void tsp::check_matrix(const std::vector<std::vector<double> > &matrix) const 
    {   
        size_t n_cols = matrix.size();
        
        for (size_t i = 0; i < n_cols; ++i) {
            size_t n_rows = matrix.at(i).size();
            // check if the matrix is square
            if (n_rows != n_cols)
                pagmo_throw(value_error, "adjacency matrix is not square");
            
            for (size_t j = 0; j < n_rows; ++j) {
                if (i == j && matrix.at(i).at(j) != 0)
                    pagmo_throw(value_error, "main diagonal elements must all be zeros.");
                if (i != j && !matrix.at(i).at(j)) // fully connected
                    pagmo_throw(value_error, "adjacency matrix contains zero values.");
                if (i != j && !matrix.at(i).at(j) == matrix.at(i).at(j)) // fully connected
                    pagmo_throw(value_error, "adjacency matrix contains NaN values.");                    
            }
        }
    }
    
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
    void tsp::objfun_impl(fitness_vector &f, const decision_vector& x) const 
    {
        size_t n = get_n_vertices();        
        pagmo_assert(x.size() == n * (n - 1));
        
        f[0]= 0;
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                if(i==j) continue;
                f[0] += m_weights[i][j] * x[compute_idx(i, j, n)];
            }
        }       
    }

    /// Constraint computation.
    /**
     * Computes the equality and inequality constraints for a decision vector
     * and returns the |c| = n(n-1)+2 constraint vector with the concatenated
     * equality constraints (n-1)(n-2) and inequality constraints.
     * The equality constraints are ordered by row, then by sum.
     * For pagmo, the final sum is set to -1.
     * @param[out] c constraint_vector
     * @param[in] x decision_vector
     */
    void tsp::compute_constraints_impl(constraint_vector &c, const decision_vector& x) const 
    {
        size_t n = get_n_vertices();

        // 1 - We set the equality constraints
        for (size_t i = 0; i < n; i++) {
            c[i] = 0;
            c[i+n] = 0;
            for (size_t j = 0; j < n; j++) {
                if(i==j) continue; // ignoring main diagonal
                decision_vector::size_type rows = compute_idx(i, j, n);
                decision_vector::size_type cols = compute_idx(j, i, n);
                c[i] += x[rows];
                c[i+n] += x[cols];
            }
            c[i] = c[i]-1;
            c[i+n] = c[i+n]-1;
        }

        //2 - We set the inequality constraints
        //2.1 - First we compute the uj (see http://en.wikipedia.org/wiki/Travelling_salesman_problem#Integer_linear_programming_formulation)
        //      we start always out tour from the first city, without loosing generality
        size_t next_city = 0,current_city = 0;
        std::vector<int> u(n);
        for (size_t i = 0; i < n; i++) {
            u[current_city] = i+1;
            for (size_t j = 0; j < n; j++) 
            {
                if (current_city==j) continue;
                if (x[compute_idx(current_city, j, n)] == 1) 
                {
                    next_city = j;
                    break;
                }
            }
            current_city = next_city;
        }
        int count=0;
        for (size_t i = 1; i < n; i++) {
            for (size_t j = 1; j < n; j++) 
            {
                if (i==j) continue;
                c[2*n+count] = u[i]-u[j] + (n+1) * x[compute_idx(i, j, n)] - n;
                count++;
            }
        }
    }

    std::string tsp::get_name() const 
    {
        return "Traveling Salesman Problem (TSP and ATSP)";
    }

}} //namespaces

//BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::tsp);
