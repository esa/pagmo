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

#ifndef PAGMO_PROBLEM_GRAPH_HELPER_HPP
#define PAGMO_PROBLEM_GRAPH_HELPER_HPP

//#include "base.h"
#include "../config.h"
#include <iostream>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>

/**
 * This class is a template helper for the boost graph adjacency list.
 * It installs two "generic" property containers for the vertices and edges.
 * It allows easier adding/removing and operations with such structures.
 * 
 * Usage:
@verbatim
    struct VertexProperties {
        int someInt;
        std::string someName;
    };
    struct EdgeProperties {
        double someCost;
    };
    typedef graph_helper<VertexProperties, EdgeProperties> MyGraph;

    MyGraph g;
    VertexProperties vp;
    vp.someInt = 42;

    MyGraph::Vertex v = g.AddVertex(vp);
    g.properties(v).someInt = 23;
 @endverbatim
 * 
 * @author Florin Schimbinschi (florinsch@gmail.com)
 */

/* definition of generic boost::graph properties */
enum vertex_properties_t { vertex_properties };
enum edge_properties_t { edge_properties };

/* install properties to the boost namespace */
namespace boost {
    BOOST_INSTALL_PROPERTY(vertex, properties);
    BOOST_INSTALL_PROPERTY(edge, properties);
}

using namespace std;
using namespace boost;

namespace pagmo {

template < typename VP, typename EP >
class graph_helper
{
public:
        /**
        * std::vector (vecS) are fastest for iterators: http://www.boost.org/doc/libs/1_55_0/libs/graph/doc/using_adjacency_list.html
        * External properties can also be added: http://www.boost.org/doc/libs/1_55_0/libs/graph/doc/quick_tour.html
        *  for hashmaps or other types, which are better for other operations.
        */
        /* an adjacency_list */
	typedef adjacency_list<
		setS, // disallow parallel edges
		listS, // Vertex container
		directedS, // directed graph
		property<vertex_properties_t, VP>,
		property<edge_properties_t, EP>
	> graph_container;

	/* graph-specific typedefs */
	typedef typename graph_traits<graph_container>::vertex_descriptor Vertex;
	typedef typename graph_traits<graph_container>::edge_descriptor Edge;
	typedef std::pair<Edge, Edge> EdgePair;

	typedef typename graph_traits<graph_container>::vertex_iterator vertex_iter;
	typedef typename graph_traits<graph_container>::edge_iterator edge_iter;
	typedef typename graph_traits<graph_container>::adjacency_iterator adjacency_iter;
	typedef typename graph_traits<graph_container>::out_edge_iterator out_edge_iter;

	typedef typename graph_traits<graph_container>::degree_size_type degree_t;

	typedef std::pair<adjacency_iter, adjacency_iter> adjacency_vertex_range_t;
	typedef std::pair<out_edge_iter, out_edge_iter> out_edge_range_t;
	typedef std::pair<vertex_iter, vertex_iter> vertex_range_t;
	typedef std::pair<edge_iter, edge_iter> edge_range_t;


	/* constructors etc. */
	graph_helper()
	{}

	graph_helper(const graph_helper& g) :
		graph(g.graph)
	{}

	virtual ~graph_helper()
	{}


	/* structure modification methods */
	void Clear()
	{
		graph.clear();
	}

	Vertex AddVertex(const VP& prop)
	{
		Vertex v = add_vertex(graph);
		properties(v) = prop;
		return v;
	}

	void RemoveVertex(const Vertex& v)
	{
		clear_vertex(v, graph);
		remove_vertex(v, graph);
	}
        
        Edge AddEdge(const Vertex& fro, const Vertex& to, const EP& props)
	{
		Edge newEdge = add_edge(fro, to, graph).first;
		properties(newEdge) = props;

		return newEdge;
	}

	EdgePair AddEdge(const Vertex& fro, const Vertex& to, const EP& prop_fwd, const EP& prop_rev)
	{
		Edge addedEdgeFwd = add_edge(fro, to, graph).first;
		Edge addedEdgeRev = add_edge(to, fro, graph).first;

		properties(addedEdgeFwd) = prop_fwd;
		properties(addedEdgeRev) = prop_rev;

		return EdgePair(addedEdgeFwd, addedEdgeRev);
	}


	/* property access */
	VP& properties(const Vertex& v)
	{
		typename property_map<graph_container, vertex_properties_t>::type param = get(vertex_properties, graph);
		return param[v];
	}

	const VP& properties(const Vertex& v) const
	{
		typename property_map<graph_container, vertex_properties_t>::const_type param = get(vertex_properties, graph);
		return param[v];
	}

	EP& properties(const Edge& v)
	{
		typename property_map<graph_container, edge_properties_t>::type param = get(edge_properties, graph);
		return param[v];
	}

	const EP& properties(const Edge& v) const
	{
		typename property_map<graph_container, edge_properties_t>::const_type param = get(edge_properties, graph);
		return param[v];
	}


	/* selectors and properties */
	const graph_container& getGraph() const
	{
		return graph;
	}

	vertex_range_t getVertices() const
	{
		return vertices(graph);
	}

	adjacency_vertex_range_t getAdjacentVertices(const Vertex& v) const
	{
		return adjacent_vertices(v, graph);
	}

	int getVertexCount() const
	{
		return num_vertices(graph);
	}

	int getVertexDegree(const Vertex& v) const
	{
		return out_degree(v, graph);
	}


	/* operators */
	graph_helper& operator=(graph_helper const& grph)
	{
		graph = grph.graph;
		return *this;
	}
        
        friend ostream& operator<<(ostream &output, graph_helper const& grph)
        { 
                output << "test";
                return output;            
        }

protected:
	graph_container graph;
private:                
            friend class boost::serialization::access;
            template <class Archive>
            void serialize(Archive &ar, const unsigned int)
            {
                    ar & graph;
            }
};
}
#endif // PAGMO_PROBLEM_GRAPH_HELPER_HPP