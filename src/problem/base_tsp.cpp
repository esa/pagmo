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
	 */
	base_tsp::base_tsp(int n) : base(n, n, 1, 0, 0, 0.0){ }
        
        /// Public
        tsp_graph const& base_tsp::get_graph() const { return m_graph; }
        
        void base_tsp::set_graph(tsp_graph const& new_graph) { m_graph = new_graph; }
        
        void base_tsp::set_graph(std::vector< std::vector<double> > const& new_graph) {
            list_to_graph(new_graph, m_graph);
        }
        
        /// Protected
        void base_tsp::list_to_graph(std::vector< std::vector<double> > const& the_list, tsp_graph& the_graph) const {            
            std::vector< std::vector<double> >::const_iterator row;
            std::vector<double>::const_iterator col;
            for (row = the_list.begin(); row != the_list.end(); ++row) {
                for (col = row->begin(); col != row->end(); ++col) {
                    // do stuff ... 
                    std::cout << *col; 
                    //boost::add_edge(std::distance(the_list.begin(),row), std::distance(row->begin(), col), *col, the_graph);
                }
            }
        }
        
        int base_tsp::get_no_vertices() const {
            return boost::num_vertices(m_graph);
        }
        
        /// Private

}} //namespaces
