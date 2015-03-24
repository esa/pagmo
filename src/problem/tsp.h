/*****************************************************************************
 *   Copyright (C) 2004-2014 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
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

#include <boost/array.hpp>
#include <vector>
#include <string>

#include "./base_tsp.h"
#include "../serialization.h"

namespace pagmo { namespace problem {

/// A static Travelling Salesman Problem
/**
 * This is a class representing the classic Travelling Salesman Problem. The problem
 * is that of finding the shortest Hamiltonian path in a weighted, bidirectional graph.
 *
 * The base_tsp::distance is thus defined as the (i,j) element of a matrix represented as
 * a std::vector<std::vector<double> >
 *
 * @author Dario Izzo (dario.izzo@gmail.com)
 * @author Annalisa Riccardi
 */

class __PAGMO_VISIBLE tsp: public base_tsp
{
    public:

        tsp();
        tsp(const std::vector<std::vector<double> >&, const base_tsp::encoding_type & = CITIES);

        /// Copy constructor for polymorphic objects (deep copy)
        base_ptr clone() const;

        const std::vector<std::vector<double> >& get_weights() const;

        /** @name Implementation of virtual methods*/
        //@{
        std::string get_name() const;
        std::string human_readable_extra() const;
        double distance(decision_vector::size_type, decision_vector::size_type) const;
        //@}

    private:
        static boost::array<int, 2> compute_dimensions(decision_vector::size_type n_cities, base_tsp::encoding_type);
        void check_weights(const std::vector<std::vector<double> >&) const;
        size_t compute_idx(const size_t i, const size_t j, const size_t n) const;

        void objfun_impl(fitness_vector&, const decision_vector&) const;
        void compute_constraints_impl(constraint_vector&, const decision_vector&) const;

        friend class boost::serialization::access;
        template <class Archive>
        void serialize(Archive &ar, const unsigned int)
        {
            ar & boost::serialization::base_object<base_tsp>(*this);
            ar & m_weights;
        }

    private:
        std::vector<std::vector<double> > m_weights;
};

}}  //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::tsp)

#endif  //PAGMO_PROBLEM_TSP_H