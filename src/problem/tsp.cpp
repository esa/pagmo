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
//        std::vector<std::vector<double> > tmp_w;
//        for (int i = 0; i < 5; ++i)
//            for (int j = 0; j < 5; ++j)
//                tmp_w[i][j] = default_weights[i][j];
//        this->set_graph(tmp_w);
        
//        int myints[] = {16,2,77,29};
//        std::vector<int> fifth (myints, myints + sizeof(myints) / sizeof(int) );
        set_lb(0);
        set_ub(4); //number of nodes/vertices in the graph -1 (we count from 0)
    }

    /// Constructor from vectors and maximum weight.
    /**
     * Initialize weights between vertices (city distances) from the matrix.
     * @param[in] weights matrix of distances between cities.
     */
    tsp::tsp(const std::vector< std::vector<double> > &weights): base_tsp(1) {
        this->set_graph(weights);
        //this-> set problem dimension
            set_lb(0);
            set_ub(weights[0].size()-1); //number of nodes in the graph -1 (we count from 0)
    }

    /// Constructor from boost graph object
    /**
     * Initialize weights of the edges (city distances) from the boost graph.
     * @param[in] boost graph object containing weights between vertices.
     */
    tsp::tsp(const tsp_graph &graph) : base_tsp(boost::num_vertices(graph)), m_graph(graph) {
        set_lb(0);
        set_ub(boost::num_vertices(graph)-1); //number of nodes in the graph
    }

    /// Clone method.
    base_ptr tsp::clone() const
    {
            return base_ptr(new tsp(*this));
    }

    /// Implementation of the objective function.
    void tsp::objfun_impl(fitness_vector &f, const decision_vector &x) const
    {
//            pagmo_assert(f.size() == 1);
//            pagmo_assert(x.size() == get_dimension() && x.size() == m_weights[0].size());
//            f[0] = 0;
//            for (size_type i = 1; i < get_dimension(); ++i) {
//                            f[0] += m_weights[boost::numeric_cast<int>(x[i-1])][boost::numeric_cast<int>(x[i])];
//            }
//            f[0] += m_weights[boost::numeric_cast<int>(x[get_dimension()-1])][boost::numeric_cast<int>(x[0])];
    }

    /// Re-implement constraint computation,
    //We check whether we have selected all the nodes (the decision vector has to be a permutation of the set of nodes).
    //The constraint is positive (not satisfied) if we have selected more than once the same node or equivalently not all 
    //the nodes have been selected
//    void tsp::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
//    {
//            if (check_partial_feasibility(x)) {
//                    c[0] = 0;
//            }
//            else {
//                    c[0] = 1;
//            }
//    }
    
//    bool tsp::check_partial_feasibility(const decision_vector &x) const {
//            (void)x; //to avoid the  unused parameter ‘x’ warning by compiler
//            return true;
//    }

    /// Extra human readable info for the problem.
    /**
     * Will return a formatted string containing the weights matrix.
     */
    std::string tsp::human_readable_extra() const
    {
//            std::ostringstream oss;
//            oss << "\nWeights Matrix: " << std::endl;
//            for(problem::base::size_type i=0; i < m_weights[0].size(); ++i) {
//                    oss << "\t\t" << m_weights[i] << std::endl;
//            }
//            return oss.str();
    }

    std::string tsp::get_name() const {
            return "Travelling Salesman Problem";
    }

    }
}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::tsp)