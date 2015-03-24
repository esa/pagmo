/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
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

#ifndef PAGMO_ALGORITHM_CSTRS_IMMUNE_SYSTEM_H
#define PAGMO_ALGORITHM_CSTRS_IMMUNE_SYSTEM_H

#include <string>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"
#include "jde.h"
#include "sga.h"

namespace pagmo { namespace algorithm {

/// Immune system constraints handling meta-algorithm
/**
 *
 * Immune constraints handling is a meta-algorithm that allow to solve constrained optimization
 * problems. The key idea of this constraint handling technique is to simulate an immune system
 * to drive the unfeasible individuals (called antibodies) towards the feasible individuals
 * (called the antigens).
 *
 * This meta-algorithm is based on the problem antibodies_problem.
 *
 * @see Hajela, P., & Lee, J. (1996). Constrained genetic search via schema adaptation: an immune
 * network solution. Structural optimization, 12(1), 11-15.
 * @see Coello, C. A. C., Cortés, N. C., San, C., & Zacatenco, P. (2001). Use of emulations of
 * the immune system to handle constraints in evolutionary algorithms.
 *
 * @author Jeremie Labroquere (jeremie.labroquere@gmail.com)
 */
		
class __PAGMO_VISIBLE cstrs_immune_system: public base
{
public:
	/// Type of immune system antibody selection method.
	/**
	* Definition of two types of antibody selections method: BEST_ANTIBODY, and INFEASIBILITY.
	* The BEST_ANTIBODY, is the selection defined by COELLO, where the antigen population
	* contains only one antigen which is the individual with lowest constraints violation.
	* INFEASIBILITY is an extension where the antigen population contains selected individuals
	* based on their infeasibility.
	*/
	// Immune system antibody selection method: best antibody, infeasibility
	enum select_method_type {BEST_ANTIBODY = 0, INFEASIBILITY = 1};

	/// Type of immune system antibody injection method.
	/**
	* Definition of two types of antibody injections method: CHAMPION, and BEST25.
	* The CHAMPION method reinject the best antibody after the immune system simulation
	* into the main population.
	* The BEST25 method reinject the best 25% antibody after the immune system simulation
	* into the main population.
	*/
	enum inject_method_type {CHAMPION = 0, BEST25 = 1};

	/// Type of antibodies problem distance method.
	/**
	 * Two distances are implemented: the hamming and euclidean distances: HAMMING, and EUCLIDEAN.
	*/
	// hamming distance, euclidean distance
	enum distance_method_type {HAMMING = 0, EUCLIDEAN = 1};

	cstrs_immune_system(const base & = jde(1), const base & = sga(), int gen = 1,
						select_method_type = BEST_ANTIBODY,
						inject_method_type = CHAMPION,
						distance_method_type = EUCLIDEAN,
						double = 0.5,
						double = 0.5,
						double = 1./3.,
						double = 1e-15, double = 1e-15);
	cstrs_immune_system(const cstrs_immune_system &);
	base_ptr clone() const;

public:
	void evolve(population &) const;
	std::string get_name() const;
	base_ptr get_algorithm() const;
	void set_algorithm(const base &);
	base_ptr get_algorithm_immune() const;
	void set_algorithm_immune(const base &);

protected:
	std::string human_readable_extra() const;

private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & m_original_algo;
		ar & m_original_algo_immune;
		ar & const_cast<int &>(m_gen);
		ar & m_select_method;
		ar & m_inject_method;
		ar & m_distance_method;
		ar & const_cast<double &>(m_phi);
		ar & const_cast<double &>(m_gamma);
		ar & const_cast<double &>(m_sigma);
		ar & const_cast<double &>(m_ftol);
		ar & const_cast<double &>(m_xtol);
	}
	base_ptr m_original_algo;
	base_ptr m_original_algo_immune;
	//Number of generations
	const int m_gen;
	// problem associated to population penalties variables
	select_method_type m_select_method;
	inject_method_type m_inject_method;
	distance_method_type m_distance_method;
	// algorithm constants
	const double m_phi;
	const double m_gamma;
	const double m_sigma;
	// tolerance
	const double m_ftol;
	const double m_xtol;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::cstrs_immune_system)

#endif // PAGMO_ALGORITHM_CSTRS_IMMUNE_SYSTEM_H
