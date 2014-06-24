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

namespace pagmo{ namespace problem {

/// Boost graph type (could also be added to types.h, for now it only has scope here.
/**
 * std::vector (vecS) are fastest for iterators: http://www.boost.org/doc/libs/1_55_0/libs/graph/doc/using_adjacency_list.html
 * External properties can also be added: http://www.boost.org/doc/libs/1_55_0/libs/graph/doc/quick_tour.html
 * for hashmaps or other types, which are better for other operations
 *
 * typedef boost::property<boost::vertex_index_t, int, boost::property<boost::vertex_name_t, std::string> > tsp_vertex_property;
 * We don't need vertice names, so an index will suffice for now:
 * boost::no_property automatically adds vertex_index_t, int for vecS
 */
typedef boost::property<boost::edge_weight_t, double> tsp_edge_properties;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, boost::no_property, tsp_edge_properties> tsp_graph;

/// Base TSP.
/**
 *
 * All integer optimization problems must extend this class in order to be solved by Ant Colony Optimization.
 *
 * m_graph is of type tsp_graph, defined above.
 * It is an adjacency list representing symmetric or assymmetric graphs, with edges representing cost/weight between vertices.
 * The internal properties of the boost graph are defined above in the tsp_edge_properties.
 * The only internal property for a vertice is its index
 * We consider all problems to be directed graphs, since for undirected it is just a matter of storage space.
 * Most of the TSPLIB problems are dense (e.g. fully connected graphs).
 *
 * @author Florin Schimbinschi (florinsch@gmail.com)
 */

class __PAGMO_VISIBLE base_tsp : public base
{
	public:
		base_tsp(int);
                /**
                 * @return the m_graph of type tsp_graph
                 */
                tsp_graph const& get_graph() const;
                /**
                 * Setters for the m_graph
                 * @param tsp_graph or list of list of doubles
                 */
                void set_graph(tsp_graph const&);
                void set_graph(std::vector< std::vector<double> > const&);

	protected:
                void list_to_graph(std::vector< std::vector<double> > const&, tsp_graph&) const;
                int get_no_vertices() const;
        private:                
		/**
		 * The boost graph, contains weights between nodes (vertices)
		 */
                tsp_graph m_graph;

//              friend class boost::serialization::access;
//		template <class Archive>
//		void serialize(Archive &ar, const unsigned int)
//		{
//			ar & boost::serialization::base_object<base>(*this);
//		}
    
};

}} //namespaces

#endif // PAGMO_PROBLEM_BASE_TSP_H