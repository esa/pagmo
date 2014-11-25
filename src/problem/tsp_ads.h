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

#ifndef PAGMO_PROBLEM_TSP_ADS_H
#define PAGMO_PROBLEM_TSP_ADS_H

#include <vector>
#include <string>

#include "../types.h"
#include "./base_tsp.h"
#include "../serialization.h"
#include "../keplerian_toolbox/planet.h"
#include "../keplerian_toolbox/planet_ss.h"

namespace pagmo { namespace problem {

/// The TSP - Asteroids-Debris Selection Problem
/**
 * This is a class representing a new variant of the TSP which refines the pagmo::problem::tsp_cs problem adding time.
 * While in the City Selection Problem (pagmo::problem::tsp_cs)
 * the objective function computes the cost of transferring between cities using a static weight matrix, here it is
 * computed using a Lambert's orbital transfer model where the transfer times are stored in the m_times data members. 
 * The distance virtual method is instead providing  a time independent approximation to such a distance.
 * 
 * The idea is that this problem can be used in a bi-level optimization where m_times is optimized in the outer layer 
 *
 * The information about arrival epochs is contained in the data member m_epochs (MJD2000) 
 * 
 * The possibility to use a waiting time \f$ WT \f$ is also provided via the data member m_waiting_time (in days)
 *
 * @author Dario Izzo (dario.izzo@gmail.com)
 */

class __PAGMO_VISIBLE tsp_ads: public base_tsp
{
    public:

        /// Constructor
        tsp_ads(
            const std::vector<kep_toolbox::planet_ptr>& planets = {kep_toolbox::planet_ss("venus").clone(), kep_toolbox::planet_ss("earth").clone(), kep_toolbox::planet_ss("venus").clone()}, 
            const std::vector<double>& values = {1.,1.,1.},
            const double max_DV = 3000, 
            const std::vector<double>&  epochs = {1200, 1550, 1940}, 
            const double waiting_time = 0., 
            const base_tsp::encoding_type & encoding = CITIES
        );

        /// Copy constructor for polymorphic objects
        base_ptr clone() const;

        /// Given an hamiltonian path finds the best subtours (accounting for the Lambert's problem)
        void find_best_selection(const decision_vector &, double &, double &, decision_vector::size_type &, decision_vector::size_type &) const;

        /** @name Getters*/
        //@{
        const std::vector<kep_toolbox::planet_ptr>& get_planets() const;
        const std::vector<double>& get_values() const;
        double get_max_DV() const;
        const decision_vector& get_epochs() const;
        double get_waiting_time() const;
        //@}

        /** @name Setters*/
        //@{
        void set_epochs(const decision_vector &);
        //@}

        /** @name Implementation of virtual methods*/
        //@{
        std::string get_name() const;
        std::string human_readable_extra() const;
        double distance(decision_vector::size_type, decision_vector::size_type) const;
        //@}
        double distance_lambert(decision_vector::size_type, decision_vector::size_type) const;

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
            ar & const_cast<std::vector<kep_toolbox::planet_ptr> &>(m_planets);
            ar & const_cast<std::vector<double> &>(m_values);
            ar & const_cast<double &>(m_max_DV);
            ar & m_epochs;
            ar & const_cast<double &>(m_waiting_time);
            ar & m_min_value;
        }

    private:
        const std::vector<kep_toolbox::planet_ptr> m_planets;
        const std::vector<double> m_values;
        const double m_max_DV;
        decision_vector m_epochs;
        const double m_waiting_time;

        // this data member is set in the constructor as the minimum in m_values and is used in the objecive function
        double m_min_value;


};

}}  //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::tsp_ads)

#endif  //PAGMO_PROBLEM_TSP_ADS_H