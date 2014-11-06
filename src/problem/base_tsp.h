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

#include <vector>

#include "./base.h"
#include "../serialization.h"
#include "../population.h"

namespace pagmo { namespace problem {

/// Base TSP (Travelling Salesman Problem).
/**
 * All pagmo::problem that are TSP variants must derive from this class
 * Algorithms such as pagmo::algorithm::inverover and pagmo::aco can solve problem deriving
 * from this class as they make use of the base_tsp::distance and the base_tsp::get_encoding methods
 *
 * The virtual method base_tsp::distance is pure and must be reimplemented by the user in the derived class
 * returning the distance between two cities
 *
 * The sequence of cities visited can be encoded in one of the following ways:
 *
 * 1-CITIES
 * This encoding represents the ids of the cities visited directly in the chromosome. e.g. [3,1,0,2]
 *
 * 2-RANDOMKEYS
 * This encoding, first introduced in the paper
 * Bean, J. C. (1994). Genetic algorithms and random keys for sequencing and optimization. ORSA journal on computing, 6(2), 154-160.
 * It essentially represents the tour as a sequence of doubles bounded in [0,1].
 * The tour is reconstructed by the argsort of the sequence. (e.g. [0.34,0.12,0.76,0.03] -> [3,1,0,2])
 *
 * 3-FULL
 * The full encoding encodes the city tour in a matrix as detailed in
 * http://en.wikipedia.org/wiki/Travelling_salesman_problem#Integer_linear_programming_formulation
 * It is used to create TSP problems that are integer linear programming problems. (e.g. [0,1,0,1,0,0,0,0,1,0,1,0] -> [0,2,3,1])
 *
 * @author Dario Izzo (dario.izzo@gmail.com)
 */

class __PAGMO_VISIBLE base_tsp: public base
{
    public:
        /// Mechanism used to encode the sequence of vertices to be visited
        enum encoding_type {
            RANDOMKEYS = 0,  ///< As a vector of doubles in [0,1].
            FULL = 1,        ///< As a matrix with ones and zeros
            CITIES = 2       ///< As a sequence of cities ids.
        };

        base_tsp(int n_cities, int nc, int nic, encoding_type = CITIES);

        /** @name Getters.*/
        //@{
        encoding_type get_encoding() const;
        decision_vector::size_type get_n_cities() const;
        //@}

        /** @name Converters between encodings.*/
        //@{
        pagmo::decision_vector full2cities(const pagmo::decision_vector &) const;
        pagmo::decision_vector cities2full(const pagmo::decision_vector &) const;
        pagmo::decision_vector randomkeys2cities(const pagmo::decision_vector &) const;
        pagmo::decision_vector cities2randomkeys(const pagmo::decision_vector &, const pagmo::decision_vector &) const;
        //@}

        // Pure virtual method returning the distance between cities
        virtual double distance(decision_vector::size_type, decision_vector::size_type) const = 0;

    private:
        friend class boost::serialization::access;
        template <class Archive>
        void serialize(Archive &ar, const unsigned int)
        {
            ar & boost::serialization::base_object<base>(*this);
            ar & const_cast<encoding_type &>(m_encoding);
            ar & const_cast<pagmo::decision_vector::size_type &>(m_n_cities);
        }

    private:
        const encoding_type m_encoding;
        const pagmo::decision_vector::size_type m_n_cities;
};

}}  //namespaces

BOOST_SERIALIZATION_ASSUME_ABSTRACT(pagmo::problem::base_tsp)

#endif  //PAGMO_PROBLEM_BASE_TSP_H