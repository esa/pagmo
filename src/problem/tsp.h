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

#ifndef PAGMO_PROBLEM_TSP_H
#define PAGMO_PROBLEM_TSP_H

#include <boost/numeric/conversion/cast.hpp>
#include <cstddef>
#include <string>
#include <vector>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "base_tsp.h"

namespace pagmo { namespace problem {

/// Traveling salesman problem
/**
 *
 * Given a list of cities and their pairwise distances, the task is to find the 
 * shortest possible tour that visits each city exactly once.
 *
 * TSP can be modeled as a graph, such that cities are the graph's vertices, 
 * paths are the graph's edges, and a path's distance is the edge's length. 
 * A TSP tour becomes a Hamiltonian cycle, and the optimal TSP tour is the 
 * shortest Hamiltonian cycle.
 *
 * TODO: rewrite here description of the optimization problem
 *
 * @see http://en.wikipedia.org/wiki/Travelling_salesman_problem
 *
 * @author Florin Schimbinschi (florinsch@gmail.com)
 */
class __PAGMO_VISIBLE tsp: public base_tsp
{
    public:
            /*NOTE: here you can call the  relevant base_tsp constructor */
            tsp(): base_tsp() {};
            tsp(tsp_graph const& graph): base_tsp(graph) {};
            tsp(vector2D<double> const& weights): base_tsp(weights) {};
            base_ptr clone() const;
            std::string get_name() const;
            
    protected:
            /*Formal definition of the objective function and equality/inequality constraints must be taken from 
            * http://en.wikipedia.org/wiki/Travelling_salesman_problem#Integer_linear_programming_formulation
            * TODO: implement
            */
            void objfun_impl(fitness_vector &, decision_vector const&) const;
            void compute_constraints_impl(constraint_vector &, decision_vector const&) const;
            
            /*TODO: implement using i*(m_n_vertices-1) + (j>i)?j-1:j and use it in objfun_impl and compute_constraints_impl*/
            decision_vector::size_type compute_idx(const std::vector<int>::size_type i, const std::vector<int>::size_type j) const;
            
            std::string human_readable_extra() const;

            
    private:
            /*NOTE: serialization of const data member added*/
            friend class boost::serialization::access;
            template <class Archive>
            void serialize(Archive &ar, const unsigned int)
            {
                    ar & boost::serialization::base_object<base>(*this);
                    ar & const_cast<vector2D<double>& >(m_weights);
            }
    private:        
            const vector2D<double> m_weights;

};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::tsp)

#endif
