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

#include <boost/array.hpp>
#include <vector>

#include "base.h"
#include "../serialization.h"

namespace pagmo{ namespace problem {

/// Base TSP.
/**
 * This is a base class for Traveling Salesman Problems. Encoding of the chromosome is as in
 * http://en.wikipedia.org/wiki/Travelling_salesman_problem#Integer_linear_programming_formulation
 *
 * It is an adjacency list representing symmetric or assymmetric graphs, with edges 
 * having internal property of weights, representing cost between vertices.
 * The only internal property for a vertex is its index.
 *
 * @author Dario Izzo (dario.izzo@gmail.com)
 * @author Annalisa Riccardi
 * @author Florin
 */

class __PAGMO_VISIBLE base_tsp: public base
{
    public:
        /// Mechanism used to transform the input problem
        enum encoding {
            RANDOMKEYS = 0, ///< The objective function of the original problem is used as objective function of the transformed problem
            FULL = 1,       ///< The sum of the constraint violation is used as objective function of the transformed problem
            CITIES = 2      ///< The sum of the constraint violation is used as objective function of the transformed problem
        };
        base_tsp();
        base_tsp(const std::vector<std::vector<double> >&, const encoding &); 
        base_ptr clone() const;

        const std::vector<std::vector<double> >& get_weights() const;
        const decision_vector::size_type& get_n_cities() const;
        const encoding& get_encoding() const;  

        pagmo::decision_vector full2cities(const pagmo::decision_vector &) const;
        pagmo::decision_vector cities2full(const pagmo::decision_vector &) const;
        pagmo::decision_vector randomkeys2cities(const pagmo::decision_vector &) const;
        pagmo::decision_vector cities2randomkeys(const pagmo::decision_vector &, const pagmo::decision_vector &) const;

        std::string get_name() const;
        std::string human_readable_extra() const;

    private:
        static boost::array<int, 6> compute_dimensions(decision_vector::size_type n_cities, encoding);
        void check_weights(const std::vector<std::vector<double> >&) const;

        void objfun_impl(fitness_vector&, const decision_vector&) const;
        void compute_constraints_impl(constraint_vector&, const decision_vector&) const;

        friend class boost::serialization::access;
        template <class Archive>
        void serialize(Archive &ar, const unsigned int)
        {
            ar & boost::serialization::base_object<base>(*this);
            ar & m_weights;
            ar & const_cast<decision_vector::size_type &>(m_n_cities);
            ar & const_cast<encoding &>(m_encoding);
        }
    
    private:
        const decision_vector::size_type m_n_cities;
        std::vector<std::vector<double> > m_weights;
        const encoding m_encoding;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::base_tsp)

#endif // PAGMO_PROBLEM_BASE_TSP_H