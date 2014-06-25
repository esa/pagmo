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
//#include "base.h"
#include "tsp.h"

namespace pagmo { namespace problem {

    static double const default_weights[5][5] = {
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
        vector2D<double> tmp;
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j)
                tmp[i][j] = default_weights[i][j];
        set_graph(tmp);
        
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
        set_lb(0);
        set_ub(boost::num_vertices(graph)); //number of nodes in the graph
    }

    /// Constructor from vectors and maximum weight.
    /**
     * Initialize weights between vertices (city distances) from the matrix.
     * @param[in] weights matrix of distances between cities.
     */
    tsp::tsp(vector2D<double> const& weights): base_tsp(1 /*!?!?! isit ???*/ ) {
        //TODO: figure out how to set problem bounds after converting
        
        set_graph(weights); 
        
        set_lb(0);
        set_ub(get_no_vertices()); // number of vertices after set_graph is run
    }

    /// Clone method.
    base_ptr tsp::clone() const
    {
        return base_ptr(new tsp(*this));
    }

    /// Implementation of the objective function
    /**
     * Computes the fitness vector associated to a decision vector.
     * @param f fitness vector
     * @param x decision vector (a permutation of the vertices)
     */
    void tsp::objfun_impl(fitness_vector &f, const decision_vector &x) const {
        pagmo_assert(f.size() == 1);
        pagmo_assert(x.size() == get_dimension() && x.size() == get_no_vertices());
        //TODO: figure this out later
//        f[0] = 0;
//        for (size_type i = 1; i < get_dimension(); ++i) {
//                        f[0] += m_weights[boost::numeric_cast<int>(x[i-1])][boost::numeric_cast<int>(x[i])];
//        }
//        f[0] += m_weights[boost::numeric_cast<int>(x[get_dimension()-1])][boost::numeric_cast<int>(x[0])];
        (void)x;
        (void)f;
    }

    /// Constraint computation.
    /**
     * The decision vector has to be a permutation of the set of nodes.
     * The constraint is positive when not satisfied:
     *  if we have selected more than once the same node or 
     *  equivalently not all the nodes have been selected.
     */
    void tsp::compute_constraints_impl(constraint_vector &c, decision_vector const& x) const
    {   //TODO: figure out if this is really needed
        
        // Checks if the length is the same
        if ( x.size() != get_no_vertices() )
            c[0] = 1;
        return;
        
        // Checks if there are duplicate items
        std::sort(x.begin(), x.end());
        if ( std::unique(x.begin(), x.end()) != x.end() )
            c[0] = 1;
        return;
        
        c[0] = 0; // we're good
    }

    /// Extra human readable info for the problem.
    /**
     * Will return a list of vertices and edges
     */
    std::string tsp::human_readable_extra() const
    {//TODO: There must be a better way of doing this..
        std::ostringstream oss;
        oss << "Adjacency List: " << std::endl;
        
        oss << m_graph; // might just work

        typedef boost::property_map<tsp_graph, boost::vertex_index_t>::type vertex_index_map;
        vertex_index_map index = get(boost::vertex_index, m_graph);

        oss << "Vertices = ";
        typedef boost::graph_traits<tsp_graph>::vertex_iterator vertex_iter;
        std::pair<vertex_iter, vertex_iter> vp;
        for (vp = boost::vertices(m_graph); vp.first != vp.second; ++vp.first)
                oss << index[*vp.first] <<  " ";
        oss << std::endl;
        
        oss << "Edges = ";
        boost::graph_traits<tsp_graph>::edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = boost::edges(m_graph); ei != ei_end; ++ei)
                oss << "(" << index[boost::source(*ei, m_graph)] 
                    << "," << index[boost::target(*ei, m_graph)] << ") ";
        oss << std::endl;
        
        
//        typedef typename boost::property_map<tsp_graph, boost::edge_weight_t>::const_type edge_map;
//        edge_map weight = boost::get(boost::edge_weight, m_graph);
//
//        typedef typename boost::property_traits<edge_map>::value_type edge_type;
//
//        edge_type source;
//
//            source = boost::get(name, source(*first, G));
        
        return oss.str();
    }

    std::string tsp::get_name() const 
    {
        return "Traveling Salesman Problem";
    }

}} //namespaces

//BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::tsp);