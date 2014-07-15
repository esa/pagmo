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

namespace pagmo { namespace problem {
        /// Public
        tsp_graph const& base_tsp::get_graph() const { return m_graph; }
        
//        void base_tsp::set_graph(tsp_graph const& new_graph) { m_graph = new_graph; }
        
        void base_tsp::set_graph(vector2D<double> const& matrix) {
//            tsp_graph new_graph;
            convert_vector2D_to_graph(matrix, m_graph);
//            m_graph = new_graph;
        }
        
        size_t const& base_tsp::get_n_vertices() const { return m_n_vertices; }
        
        void base_tsp::convert_vector2D_to_graph(vector2D<double> const& the_vector, tsp_graph& the_graph) {
            tsp_edge_map_weight weights = boost::get(boost::edge_weight_t(), the_graph);
            tsp_vertex from, to;
            tsp_edge link;
            
            // add vertices first 
            /* Checking if a vertex exists with no vertices inserted causes segfault
             * so we have to iterate 1st to get total number of vertices
             * then iterate again to insert them ... bummer
             * Couldn't figure out how to do it all in 2 for loops
             */
            int no_vertices = count_vertices(the_vector);
            for (int v = 0; v < no_vertices; ++v)
                boost::add_vertex(v, the_graph);
            
            // add edges and weights
            for (size_t i = 0; i < the_vector.size(); ++i) {

                /* uncomment this and it's segfault
                 * don't do this check and the logic is wrong
                 */ 
                from = boost::vertex(i, the_graph);
//                if (from == tsp_graph::null_vertex())
//                    from = boost::add_vertex(i, the_graph);
                
                for (size_t j = 0 ; j < (the_vector.at(i)).size(); ++j) {
                    // we don't allow connections to self
                    if(i == j) continue;
                    
                    to = boost::vertex(j, the_graph);
                    // create destination vertex only if not existent
                    // for some reason this works, but is not enough
//                    if (to == tsp_graph::null_vertex())
//                        to = boost::add_vertex(j, the_graph);
                    
                    // create an edge connecting those two vertices
                    link = (boost::add_edge(from, to, the_graph)).first;
                    // add weight property to the edge
                    weights[link] = the_vector.at(i).at(j);
                }
            }
        }
        
        void base_tsp::convert_graph_to_vector2D(tsp_graph const& the_graph, vector2D<double>& the_vector) {
            tsp_vertex_map_const_index vtx_idx = boost::get(boost::vertex_index_t(), the_graph);
            tsp_edge_map_const_weight weights = boost::get(boost::edge_weight_t(), the_graph);
            tsp_edge_range_t e_it = boost::edges(the_graph);
            
            for (e_it = boost::edges(the_graph); e_it.first != e_it.second; ++e_it.first) {
                int i = vtx_idx[boost::source(*e_it.first, the_graph)];
                int j = vtx_idx[boost::target(*e_it.first, the_graph)];
                the_vector[i][i] = 0;
                the_vector[i][j] = weights[*e_it.first];
            }
        }
        
        size_t base_tsp::count_vertices(vector2D<double> const& the_vector) {
            unsigned int maxRow = the_vector.size();
            unsigned int maxCol = 0;
            vector2D<double>::const_iterator row;
            for (row = the_vector.begin(); row != the_vector.end(); ++row)
                if(row->size() > maxCol)
                    maxCol = row->size();
            
            return (int)(maxRow > maxCol ? maxRow : maxCol);
        }
        
        /// Protected
        
        /// Private

}} //namespaces
