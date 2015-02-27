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

#ifndef PAGMO_PROBLEM_TSP_DS_H
#define PAGMO_PROBLEM_TSP_DS_H

#include <vector>
#include <string>

#include "../types.h"
#include "./base_tsp.h"
#include "../serialization.h"
#include "../keplerian_toolbox/planets/planet.h"
#include "../keplerian_toolbox/planets/planet_ss.h"
#include "../keplerian_toolbox/astro_constants.h"

namespace pagmo { namespace problem {

/// The TSP - Debris Selection Problem
/**
 * This is a class representing a new variant of the TSP where the cities are orbiting
 * objects and must be visited with a predefined schedule. 
 * 
 * The cost of transferring from one orbiting object to the next is computed by means of a three impulse approximation
 * that uses the orbital element at epoch.
 *
 * The distance virtual method is instead providing the same computations using the initial orbital elements (i.e.
 * not propagating them)
 * 
 * The information about the schedule is contained in the data member m_epochs (MJD2000) 
 *
 * @author Dario Izzo (dario.izzo@gmail.com)
 */

class __PAGMO_VISIBLE tsp_ds: public base_tsp
{
    public:

        /// Constructor
        tsp_ds(
            const std::vector<kep_toolbox::planet_ptr>& planets = {kep_toolbox::planet_ss("venus").clone(), kep_toolbox::planet_ss("earth").clone(), kep_toolbox::planet_ss("mars").clone()}, 
            const std::vector<double>& values = {1.,1.,1.},
            const double max_DV = 30000, 
            const std::vector<double>&  epochs = {1200, 1550, 1940}, 
            const base_tsp::encoding_type & encoding = CITIES
        );

        /// Copy constructor for polymorphic objects
        base_ptr clone() const;

        /// Given an hamiltonian path finds the best subtour
        void find_subsequence(const decision_vector &, double &, double &, decision_vector::size_type &, decision_vector::size_type &, const bool = false) const;

        /** @name Getters*/
        //@{
        const std::vector<kep_toolbox::planet_ptr>& get_planets() const;
        const std::vector<double>& get_values() const;
        double get_max_DV() const;
        const decision_vector& get_epochs() const;
        //@}

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

        double three_impulses(double, double, double, double, double, double, double, double) const;
        void precompute_ephemerides() const;
        void compute_DVs(const decision_vector&, bool = false) const;
        double distance_3_imp(const decision_vector::size_type, const decision_vector::size_type, const size_t) const;

        friend class boost::serialization::access;
        template <class Archive>
        void serialize(Archive &ar, const unsigned int)
        {
            ar & boost::serialization::base_object<base_tsp>(*this);
            ar & const_cast<std::vector<kep_toolbox::planet_ptr> &>(m_planets);
            ar & const_cast<std::vector<double> &>(m_values);
            ar & const_cast<double &>(m_max_DV);
            ar & const_cast<std::vector<double> &>(m_epochs);
            precompute_ephemerides();
        }

    private:
        const std::vector<kep_toolbox::planet_ptr> m_planets;
        const std::vector<double> m_values;
        const double m_max_DV;
        const decision_vector m_epochs;
        const double m_mu;

        // These are to pre-allocate memory
        mutable std::vector<double> m_DV;
        mutable std::vector<std::vector<kep_toolbox::array3D> >m_eph_r;
        mutable std::vector<std::vector<kep_toolbox::array3D> >m_eph_v;
        mutable std::vector<std::vector<kep_toolbox::array6D> >m_eph_el;
};

}}  //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::tsp_ds)

#endif  //PAGMO_PROBLEM_TSP_DS_H