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

#include "base_tsp.h"
#include "tsp.h"

namespace pagmo { namespace problem {

    /**
     * The default constructor
     * This constructs a 3-cities symmetric problem (naive TSP) 
     * with weight matrix [[0,1,1][1,0,1][1,1,0]]
     */
    base_tsp::base_tsp(): base(6, 6, 1, 8, 2, 0.0), m_n_vertices(3), m_graph()
    {
        tsp_vertex from, to;
        tsp_edge link;
        tsp_edge_map_weight weights = boost::get(boost::edge_weight_t(), m_graph);
        
        boost::add_vertex(0, m_graph);
        boost::add_vertex(1, m_graph);
        boost::add_vertex(2, m_graph);
        
        for (size_t i = 0; i < 3; ++i) {
            from = boost::vertex(i, m_graph);
            for (size_t j = 0 ; j < 3; ++j) {
                if(i == j) continue; // no connections from a vertex to self
                to = boost::vertex(j, m_graph);
                // create an edge connecting those two vertices
                link = (boost::add_edge(from, to, m_graph)).first;
                // add weight property to the edge
                weights[link] = 1.23;
            }
        }
        set_lb(0);
        set_ub(1);
    }

    /**
     * Constructor from a tsp_graph object
     * @param[in] tsp_graph
     */
    base_tsp::base_tsp(const tsp_graph& graph): 
        base(
            boost::num_vertices(graph)*(boost::num_vertices(graph)-1), 
            boost::num_vertices(graph)*(boost::num_vertices(graph)-1), 
            1, 
            boost::num_vertices(graph)*(boost::num_vertices(graph)-1)+2,
            (boost::num_vertices(graph)-1)*(boost::num_vertices(graph)-2),
            0.0
        ),
        m_n_vertices(boost::num_vertices(graph)),
        m_graph(graph)
    {
        set_lb(0);
        set_ub(1);
    }

    /**
     * Getter for the m_graph
     * @return reference to the m_graph of type tsp_graph
     */
    const tsp_graph& base_tsp::get_graph() const 
    { 
        return m_graph; 
    }
    
    /**
     * Getter for the m_n_vertices
     * @return reference to the number of vertices in the graph
     */
    const size_t& base_tsp::get_n_vertices() const
    { 
        return m_n_vertices; 
    }

    /// Extra human readable info for the problem.
    /**
     * @return a std::string containing a list of vertices and edges
     */
    std::string base_tsp::human_readable_extra() const 
    {
        tsp_graph const the_graph = get_graph();
        
        std::ostringstream oss;
        oss << "\n\tThe Boost Graph (Adjacency List): \n";// << the_graph << std::endl;

        tsp_vertex_map_const_index vtx_idx = boost::get(boost::vertex_index_t(), the_graph);
        tsp_edge_map_const_weight weights = boost::get(boost::edge_weight_t(), the_graph);

        oss << "\tVertices = { ";
        
        tsp_vertex_range_t v_it;
        for (v_it = boost::vertices(the_graph); v_it.first != v_it.second; ++v_it.first)
                oss << vtx_idx[*v_it.first] <<  " ";
        oss << "}" << std::endl;
        
        oss << "\tEdges (Source, Target) = Weight : " << std::endl;
        
        int count = 0;
        for (tsp_edge_range_t e_it = boost::edges(the_graph); e_it.first != e_it.second; ++e_it.first) {
            int i = vtx_idx[boost::source(*e_it.first, the_graph)];
            int j = vtx_idx[boost::target(*e_it.first, the_graph)];
            oss << "\t(" << i << ", " << j<< ") = " << weights[*e_it.first] << std::endl;
            count++;
            if (count > 5) break;
        }
        oss << std::endl;
        
        return oss.str();
    }

    
}} //namespaces