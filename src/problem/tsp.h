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

#ifndef PAGMO_PROBLEM_tsp_H
#define PAGMO_PROBLEM_tsp_H

#include <boost/array.hpp>
#include <vector>

#include "base.h"
#include "../serialization.h"

namespace pagmo{ namespace problem {

/// Base TSP.
/**
 * This is a class representing Travelling Salesman Problems. 
 * The problem encoding can be of three different types
 *
 * 1-CITIES
 * This encoding represents the ids of the cities visited directly in the chromosome. It will
 * thus create a constrained problem as only permutation of the cities ids are valid (e.g. [0,2,1,5,0] is not
 * a valid chromosome)
 *
 * 2-RANDOMKEYS
 * This encoding, first introduced in the paper
 * Bean, J. C. (1994). Genetic algorithms and random keys for sequencing and optimization. ORSA journal on computing, 6(2), 154-160.
 * creates a box constrained problem without any constraint. It essentially represents the tour as a sequence of doubles bounded in [0,1].
 * The tour is reconstructed by the argsort of the sequence. (e.g. [0.34,0.12,0.76,0.03] -> [3,1,0,2])
 *
 * 3-FULL
 * In the full encoding the TSP is represented as a integer linear programming problem. The details can be found in
 * http://en.wikipedia.org/wiki/Travelling_salesman_problem#Integer_linear_programming_formulation
 *
 *
 * @author Dario Izzo (dario.izzo@gmail.com)
 * @author Annalisa Riccardi
 */

class __PAGMO_VISIBLE tsp: public base
{
    public:
        /// Mechanism used to transform the input problem
        enum encoding {
            RANDOMKEYS = 0, ///< The city tour is encoded as a vector of doubles in [0,1]. Results in an unconstrained box-bounded problem
            FULL = 1,       ///< The TSP is encoded ads a linear integer programming problem
            CITIES = 2      ///< The city teour is directly encoded in the chromosome as a sequence of cities ids.
        };
        tsp();
        tsp(const std::vector<std::vector<double> >&, const encoding & = FULL); 
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

BOOST_CLASS_EXPORT_KEY(pagmo::problem::tsp)

#endif // PAGMO_PROBLEM_tsp_H