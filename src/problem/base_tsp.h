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
#include "../serialization.h"

#include <vector>
#include <boost/graph/graphviz.hpp>



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

/// Base TSP.
/**
 * This is a base class for Traveling Salesman Problems. Encoding of the chromosome is as in
 * http://en.wikipedia.org/wiki/Travelling_salesman_problem#Integer_linear_programming_formulation
 * m_graph is of type tsp_graph, defined in base_tsp.h
 * It is an adjacency list representing symmetric or assymmetric graphs, with edges 
 * having internal property of weights, representing cost between vertices.
 * The only internal property for a vertex is its index.
 *
 * @author Florin Schimbinschi (florinsch@gmail.com)
 */

class __PAGMO_VISIBLE base_tsp: public base
{
    public:
            base_tsp();
            base_tsp(const tsp_graph&);    

            const tsp_graph& get_graph() const;
            const size_t& get_n_vertices() const;               
            std::string human_readable_extra() const;

    private:
            friend class boost::serialization::access;
            template <class Archive>
            void serialize(Archive &ar, const unsigned int)
            {
                    ar & boost::serialization::base_object<base>(*this);
                    ar & m_graph;
                    m_n_vertices = boost::num_vertices(m_graph);
            }
    
    private:
            size_t m_n_vertices;
            tsp_graph m_graph;
};

}} //namespaces

BOOST_SERIALIZATION_ASSUME_ABSTRACT(pagmo::problem::base_tsp)

#endif // PAGMO_PROBLEM_BASE_TSP_H