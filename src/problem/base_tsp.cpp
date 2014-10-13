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

#include <algorithm> 

#include "base_tsp.h"

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

    /// Transforms a tsp chromosome into the sequence of city indexes
    /**
     * @param[in] x the chromosome that represents a city tour
     * @return a vector containing the indices of the visited cities in the encoded order
     */
    std::vector<pagmo::population::size_type> base_tsp::chromosome2cities(const pagmo::decision_vector &x) const
    {
        if (feasibility_x(x)) 
        {
            std::vector<pagmo::population::size_type> retval(m_n_vertices,0);
            pagmo::population::size_type next_city,cur_city = 0;
            retval[0]=cur_city;
            for(size_t j = 1; j < m_n_vertices; j++){
                next_city = std::find(x.begin() + cur_city*(m_n_vertices-1), x.begin() + (cur_city+1)*(m_n_vertices-1),1) - (x.begin() + cur_city*(m_n_vertices-1));
                std::cout << cur_city << ", " << next_city << ", " << (next_city <= cur_city) << std::endl;
                next_city = next_city  + ( (next_city >= cur_city) ? 1:0 );
                cur_city=next_city;
                retval[j] = next_city;
            }
            return retval;
        } else {
            pagmo_throw(value_error,"chromosome is not compatible with the problem");
        }
    }

    /// Transforms a permutation of city indexes into a tsp chromosome
    /**
     * @param[in] vities the chromosome that represents a city tour
     * @return a vector containing the indices of the visited cities in the encoded order
     */
    pagmo::decision_vector base_tsp::cities2chromosome(const std::vector<population::size_type> &x) const
    {
        if (x.size() != m_n_vertices) 
        {
            pagmo_throw(value_error,"city indexes are of incompatible length");
        }
        std::vector<population::size_type> range(m_n_vertices);
        std::iota(range.begin(),range.end(),0);
        if (!std::is_permutation(x.begin(),x.end(),range.begin()) )
        {
            pagmo_throw(value_error,"city indexes are not a permutation of 0,1,2,3,....");
        }
        pagmo::decision_vector retval(m_n_vertices*(m_n_vertices-1),0);
        for (std::vector<population::size_type>::size_type i=0; i<x.size()-1; ++i)
        {
            retval.at( x[i]*(m_n_vertices-1) + x[i+1] - (x[i+1]>=x[i]?1:0) ) = 1;
        } 
        retval[ x.at(x.size()-1)*(m_n_vertices-1) + x[0] + (x[0]>=x[x.size()-1]?1:0) ] = 1;
        return retval;
    }

    /**
     * Getter for the m_graph
     * @return reference to the m_graph of type tsp_graph
     */
    const base_tsp::tsp_graph& base_tsp::get_graph() const 
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
        int count = 0;
        for (v_it = boost::vertices(the_graph); v_it.first != v_it.second; ++v_it.first)
        {
                oss << vtx_idx[*v_it.first] <<  " ";
                count++;
                if (count > 5) {
                oss << "... ";
                    break;
                }
        }
        oss << "}" << std::endl;
        
        oss << "\tEdges (Source, Target) = Weight: " << std::endl;
        
        count =0;
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