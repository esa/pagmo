/*****************************************************************************
 *   Copyright (C) 2004-2014 The PaGMO development team,                     *
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

#ifndef PAGMO_PROBLEM_BASE_TSP_H
#define PAGMO_PROBLEM_BASE_TSP_H

#include "base.h"
#include <vector>
#include <boost/graph/graphviz.hpp>
//#include "graph_helper.hpp"

namespace pagmo{ namespace problem {

/// Boost graph type (could also be added to types.h, for now it only has scope here.
/**
 * std::vector (vecS) are fastest for iterators: http://www.boost.org/doc/libs/1_55_0/libs/graph/doc/using_adjacency_list.html
 * External properties can also be added: http://www.boost.org/doc/libs/1_55_0/libs/graph/doc/quick_tour.html
 * for hashmaps or other types, which are better for other operations.
 */
/* Graph Internal Properties */    
typedef boost::property<boost::vertex_index_t, int 
    /*, boost::property<boost::vertex_name_t, std::string> */> tsp_vertex_properties;
typedef boost::property<boost::edge_index_t, int, // used for external properties
        boost::property<boost::edge_weight_t, double> > tsp_edge_properties;

/* Graph Type, internal containers for edges and lists */
typedef boost::adjacency_list<  
            boost::setS, // disallow parallel edges
            boost::listS, // vertex container
            boost::directedS, // directed graph
            tsp_vertex_properties, 
            tsp_edge_properties
        > tsp_graph;

/* Graph specific typedefs */
typedef boost::graph_traits<tsp_graph>::vertex_descriptor tsp_vertex;
typedef boost::graph_traits<tsp_graph>::edge_descriptor tsp_edge;

typedef boost::property_map<tsp_graph, boost::vertex_index_t>::type tsp_vertex_map_index;
typedef boost::property_map<tsp_graph, boost::vertex_index_t>::const_type tsp_vertex_map_const_index;
typedef boost::property_map<tsp_graph, boost::edge_weight_t>::type tsp_edge_map_weight;
typedef boost::property_map<tsp_graph, boost::edge_weight_t>::const_type tsp_edge_map_const_weight;
typedef boost::property_map<tsp_graph, boost::edge_index_t>::type tsp_edge_map_index;
typedef boost::property_map<tsp_graph, boost::edge_index_t>::const_type tsp_edge_map_const_index;

typedef boost::graph_traits<tsp_graph>::vertex_iterator tsp_vertex_iter;
typedef boost::graph_traits<tsp_graph>::edge_iterator tsp_edge_iter;

typedef std::pair<tsp_edge, tsp_edge> tsp_edge_pair;
typedef std::pair<tsp_vertex_iter, tsp_vertex_iter> tsp_vertex_range_t;
typedef std::pair<tsp_edge_iter, tsp_edge_iter> tsp_edge_range_t;

typedef boost::iterator_property_map<double*, tsp_edge_map_index, double, double&> tsp_ext_weight_iterator;

/// Shortened version of the 2D vector type
template <typename T>
using vector2D = typename std::vector<std::vector<T> >;

/// Base TSP.
/**
 * This is a base class for Traveling Salesman Problems. Encoding of the chromosome is as in
 * http://en.wikipedia.org/wiki/Travelling_salesman_problem#Integer_linear_programming_formulation
 * m_graph is of type tsp_graph, defined above.
 * It is an adjacency list representing symmetric or assymmetric graphs, with edges representing cost/weight between vertices.
 * The internal properties of the boost graph are defined above in the tsp_edge_properties.
 * The only internal property for a vertice is its index.
 * All problems are considered to be directed graphs, since for undirected it's just a matter of storage space.
 * Most of the TSPLIB problems are dense (e.g. fully connected graphs).
 *
 * @author Florin Schimbinschi (florinsch@gmail.com)
 */

class __PAGMO_VISIBLE base_tsp: public base
{
    public:
            /**
             * The default constructor
             * This constructs a 3-cities symmetric problem (naive TSP) 
             * with matrix [[0,1,1][1,0,1][1,1,0]]
             * base(n*(n-1), n*(n-1), 1, 2*n, (n-1)*(n-2), 0.0)
             */
            base_tsp(): base(6, 6, 1, 6, 2, 0.0), m_n_vertices(3) 
            {            
                tsp_vertex from, to;
                
                boost::add_vertex(0, m_graph);
                boost::add_vertex(1, m_graph);
                boost::add_vertex(2, m_graph);
                
                for (size_t i = 0; i < 3; ++i) {
                    from = boost::vertex(i, m_graph);
                    for (size_t j = 0 ; j < 3; ++j) {
                        if(i == j) continue; // no connections from a vertex to self
                        to = boost::vertex(j, m_graph);
                        // create an edge connecting those two vertices
                        boost::add_edge(from, to, m_graph);
                    }
                }
                set_lb(0);
                set_ub(1);
            };
            
            /**
             * Constructor from a tsp_graph object
             * @param[in] tsp_graph
             */
            base_tsp(tsp_graph const& graph): 
                base(
                        boost::num_vertices(graph)*(boost::num_vertices(graph)-1), 
                        boost::num_vertices(graph)*(boost::num_vertices(graph)-1), 
                        1, 
                        2*boost::num_vertices(graph),
                        (boost::num_vertices(graph)-1)*(boost::num_vertices(graph)-2),
                        0.0
                    ),
                    m_n_vertices(boost::num_vertices(graph)),
                    m_graph(graph)
            {
                set_lb(0);
                set_ub(1);
            };
            
            /**
             * Constructor from a vector2D
             * @param[in] vector2D
             */
            base_tsp(vector2D<double> const& weights): 
                base(
                        count_vertices(weights)*(count_vertices(weights)-1), 
                        count_vertices(weights)*(count_vertices(weights)-1), 
                        1, 
                        2*count_vertices(weights), 
                        (count_vertices(weights)-1)*(count_vertices(weights)-2), 
                        0.0
                    ), 
                    m_n_vertices(count_vertices(weights)) 
            {
                set_graph(weights);
                set_lb(0);
                set_ub(1);
            };
            
            /**
             * Getter for the m_graph
             * @return reference to the m_graph of type tsp_graph
             */
            tsp_graph const& get_graph() const;
            
            /**
             * Setter for the m_graph
             * @param[in] tsp_graph
             */
            void set_graph(tsp_graph const&);
            
            /**
             * Setter for the m_graph
             * @param[in] 2D vector of doubles
             */
            void set_graph(vector2D<double> const&);
            
            /**
             * Returns the number of vertices in the graph
             */
            size_t const& get_n_vertices() const;
            
            /**
             * Converts a 2D vector to a boost graph type tsp_graph
             * @param[in] the vector2D of doubles
             * @param[out] the tsp_graph adjacency list
             */
            static void convert_vector2D_to_graph(vector2D<double> const&, tsp_graph&);
            
            /**
             * Converts a graph back to a vector2D
             * @param[in] tsp_graph object
             * @param[out] vector2D of doubles
             * @return vector2D
             */
            static void convert_graph_to_vector2D(tsp_graph const&, vector2D<double>&);
            
            /**
             * Checks the maximum dimensions for both width and height.
             * The 'matrix' might not be square.
             * Returns the number of vertices in the 2D vector,
             * which is either max(row, col) for sparse matrices.
             * @return number of vertices
             */
            static size_t count_vertices(vector2D<double> const&);
           
    private:
            friend class boost::serialization::access;
            template <class Archive>
            void serialize(Archive &ar, const unsigned int)
            {
                    ar & boost::serialization::base_object<base>(*this);
                    ar & const_cast<tsp_graph& >(m_graph);
                    ar & const_cast<size_t& >(m_n_vertices);
            }
    
    private:
            /**
             * The number of vertices in the graph.
             * Stored here for not recomputing
             */
            const size_t m_n_vertices;
            /**
             * NOTE FLORIN: can't have this as const, need to convert it
             *                  from a matrix do a tsp_graph
             * The boost graph, an adjacency list.
             * Derived classes inherit this property
             */
            tsp_graph m_graph;
};

}} //namespaces

#endif // PAGMO_PROBLEM_BASE_TSP_H