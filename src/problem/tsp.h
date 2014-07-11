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
 * This is a constrained integer single-objective problem.
 *
 * Given a list of cities and their pairwise distances, the task is to find the 
 * shortest possible tour that visits each city exactly once.
 *
 * TSP can be modeled as a graph, such that cities are the graph's vertices, 
 * paths are the graph's edges, and a path's distance is the edge's length. 
 * A TSP tour becomes a Hamiltonian cycle, and the optimal TSP tour is the 
 * shortest Hamiltonian cycle.
 *
 * In PaGMO's terminology, this problem has global and integer dimensions equal to N 
 * (the number of vertices of the graph), fitness dimension equal to 1 (cost of the path).
 * ??? global and inequality constraints dimensions equal to 1 (to check whether the solution is valid). 
 * A valid decision vector is a permutation of the edges. 
 * The optmial solution minimizes the cost of the path S(x)
 * 
 * \f[
 * 	\textnormal{minimize:} S(x) = w_{x_n,1} + \sum_{i=1}{n} w_{x_i, x_{i+1}} 
 * \f]
 *
 * @see http://en.wikipedia.org/wiki/Travelling_salesman_problem
 *
 * @author Florin Schimbinschi (florinsch@gmail.com)
 */
class __PAGMO_VISIBLE tsp: public base_tsp
{
    public:
            /**
             * The default constructor
             */
            tsp();
            /**
             * Constructor from a tsp_graph object
             * @param[in] tsp_graph
             */
            tsp(tsp_graph const&);
            /**
             * Constructor from a vector2D
             * @param[in] vector2D
             */
            tsp(vector2D<double> const&);
            base_ptr clone() const;
            std::string get_name() const;
            
    protected:
            void objfun_impl(fitness_vector &, decision_vector const&) const;
            void compute_constraints_impl(constraint_vector &, decision_vector const&) const;
            std::string human_readable_extra() const;
//            void convert_decision_vector_to_vector2D(decision_vector const&, vector2D<bool> &);
//            void convert_vector2D_to_decision_vector(vector2D<bool> const&, decision_vector &);
            
    private:
            vector2D<double> m_weights;
            
            friend class boost::serialization::access;
            template <class Archive>
            void serialize(Archive &ar, const unsigned int)
            {
                    ar & boost::serialization::base_object<base>(*this);
                    ar & m_weights;
            }
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::tsp)

#endif
