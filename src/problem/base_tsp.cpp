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
	/**
	 * @param[in] n integer dimension of the problem (no vertices)
         * 3rd parameter from base is dimension of the fitness_vector (1)
         * Global and inequality constraints dimensions are equal to 0
         * Constraints tolerance is 0.0
	 */
	base_tsp::base_tsp(int n): base(n, n, 1, 0, 0, 0.0){ }
        
        /// Public
        tsp_graph const& base_tsp::get_graph() const { return m_graph; }
        
        void base_tsp::set_graph(tsp_graph const& new_graph) { m_graph = new_graph; }
        
        void base_tsp::set_graph(vector2D<double> const& new_graph) {
            vector2D_to_graph(new_graph, m_graph);
        }
        
        /// Protected
        void base_tsp::vector2D_to_graph(vector2D<double> const& the_vector, tsp_graph& the_graph) {
            tsp_edge_map_weight weights = boost::get(boost::edge_weight_t(), the_graph);
            tsp_vertex from, to;
            tsp_edge link;
            
            // add vertices first 
            /* Checking if a vertex exists with no vertices inserted causes segfault
             * so we have to iterate 1st to get total number of vertices
             * then iterate again to insert them ... bummer
             * Couldn't figure out how to do it all in 2 for loops
             */
            int no_vertices = get_no_vertices(the_vector);
            for (int v = 0; v < no_vertices; ++v)
                boost::add_vertex(v, the_graph);
            
            // add edges and weights
            for (size_t i = 0; i < the_vector.size(); ++i) {

                /* uncomment this and it's segfault
                 * don't do this check and the logic is wrong
                 */ 
//                from = boost::vertex(i, the_graph);
//                if (from == tsp_graph::null_vertex())
                    from = boost::add_vertex(i, the_graph);
                
                for (size_t j = 0 ; j < (the_vector.at(i)).size(); ++j) {
                    // we don't allow connections to self
                    if(i == j) continue;
                    
//                    to = boost::vertex(j, the_graph);
                    // create destination vertex only if not existent
                    // for some reason this works, but is not enough
//                    if (to == tsp_graph::null_vertex())
                        to = boost::add_vertex(j, the_graph);
                    
                    // create an edge connecting those two vertices
                    link = (boost::add_edge(from, to, the_graph)).first;
                    // add weight property to the edge
                    weights[link] = the_vector.at(i).at(j);
                }
            }
        }
        
        int base_tsp::get_no_vertices(vector2D<double> const& the_vector) const {
            unsigned int maxRow = the_vector.size();
            unsigned int maxCol = 0;
            vector2D<double>::const_iterator row;
            for (row = the_vector.begin(); row != the_vector.end(); ++row)
                if(row->size() > maxCol)
                    maxCol = row->size();
            
            return (int)(maxRow > maxCol ? maxRow : maxCol);
        }
        
        /// Private

}} //namespaces
