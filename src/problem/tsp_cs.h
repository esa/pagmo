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

#ifndef PAGMO_PROBLEM_TSP_CS_H
#define PAGMO_PROBLEM_TSP_CS_H

#include <boost/array.hpp>
#include <vector>
#include <string>

#include "./base_tsp.h"
#include "../serialization.h"

namespace pagmo { namespace problem {

/// The City-Selection Travelling Salesman Problem
/**
 * This is a class representing a new variant of the TSP which was firstly introduced (by us)
 * in the context of asteroid/space debris selection problems. 
 * 
 * Given a fully connected graph (V,E), where E is a list of weighted edges and V a list of valued vertices, the problem is that
 * of finding a path of length \f$ l<l^* \f$ accumulating the largest possible vertex value.
 * We encode a solution to this problem as an Hamiltonian path and evaluate its merit finding, along it, 
 * the best sequence satisfying the \f$ l<l^* \f$ constraint. Denoting its cumulated vertex value with \f$ P \f$ we also compute, 
 * along it, \f$ \epsilon = |l-l^*| \in [0,1] \f$. 
 * The problem will then be that of finding the Hamiltonian path maximizing \f$ P \f$ and having the largest \f$ \frac{\epsilon}{l^*} \f$.
 *
 * NOTE: if \f$ l^* \f$ is larger than than the shortest Hamiltonian path (i.e. the solution to the static TSP), the solution
 * to this problem is the solution to the static TSP. If \f$ l^* \f$ is, instead, smaller than the minimum weight across edges,
 * the solution to this problem is trivial and its global minimum is exactly \f$ \max V \f$ (corrsponding to \f$ P=\max V  \f$ and \f$ \epsilon = l^* \f$).
 *
 * @author Dario Izzo (dario.izzo@gmail.com)
 */

class __PAGMO_VISIBLE tsp_cs: public base_tsp
{
    public:

        /// Constructors
        tsp_cs();
        tsp_cs(const std::vector<std::vector<double> >&, const std::vector<double>&, const double, const base_tsp::encoding_type & = CITIES);

        /// Copy constructor for polymorphic objects
        base_ptr clone() const;

        /** @name Getters*/
        //@{
        const std::vector<std::vector<double> >& get_weights() const;
        const std::vector<double>& get_values() const;
        double get_max_path_length() const;
        //@}

        /** @name Implementation of virtual methods*/
        //@{
        std::string get_name() const;
        std::string human_readable_extra() const;
        double distance(decision_vector::size_type, decision_vector::size_type) const;
        //@}

        void find_subsequence(const decision_vector &, double &, double &, decision_vector::size_type &, decision_vector::size_type &) const;

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
            ar & const_cast<double &>(m_max_path_length);
	    ar & m_min_value;
	    ar & m_values;
        }

    private:
        std::vector<std::vector<double> > m_weights;
        std::vector<double> m_values ;
        const double m_max_path_length;
        double m_min_value;
};

}}  //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::tsp_cs)

#endif  //PAGMO_PROBLEM_TSP_CS_H
