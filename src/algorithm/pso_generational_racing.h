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

#ifndef PAGMO_ALGORITHM_PSO_GENERATIONAL_RACING_H
#define PAGMO_ALGORITHM_PSO_GENERATIONAL_RACING_H

#include "../config.h"
#include "../serialization.h"
#include "base.h"
#include "../util/race_pop.h"
#include <limits>

namespace pagmo { namespace algorithm {

/// Particle Swarm optimization generational with racing
/**
 * Compared to the generational version of Particle Swarm Optimization, this
 * PSO is further extended with racing mechanisms for individuals.
 *
 */

class __PAGMO_VISIBLE pso_generational_racing: public base
{
public:
	pso_generational_racing(int gen=1, double omega = 0.7298, double eta1 = 2.05, double eta2 = 2.05, double vcoeff = 0.5, int variant = 5, int neighb_type = 2, int neighb_param = 4, unsigned int nr_eval_per_x = 5, unsigned int max_fevals = std::numeric_limits<unsigned int>::max());
	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;
protected:	
	std::string human_readable_extra() const;
private:
	decision_vector particle__get_best_neighbor( population::size_type pidx, std::vector< std::vector<int> > &neighb, const std::vector<decision_vector> &lbX, const std::vector<fitness_vector> &lbfit, const problem::base &prob) const;
	void initialize_topology__gbest( const population &pop, decision_vector &gbX, fitness_vector &gbfit, std::vector< std::vector<int> > &neighb ) const;
	void initialize_topology__lbest( std::vector< std::vector<int> > &neighb ) const;
	void initialize_topology__von( std::vector< std::vector<int> > &neighb ) const;
	void initialize_topology__adaptive_random( std::vector< std::vector<int> > &neighb ) const;

	// Handles the averaging of stochastic fitness
	void particle__average_fitness(const problem::base &, fitness_vector &, const decision_vector &) const;
	//void particle__average_fitness_set_best(const problem::base &, std::vector<fitness_vector &, const decision_vector &) const;
	void particle__average_fitness_set_best(const problem::base &base, std::vector<fitness_vector> &F, population::size_type& best_idx, fitness_vector& best_fit, const std::vector<decision_vector> &X) const;

	// Helper routines for racing related features
	decision_vector particle__racing_get_best_neighbor( population::size_type pidx, std::vector< std::vector<int> > &neighb, const std::vector<decision_vector> &lbX, util::racing::race_pop& ) const;
	void racing__construct_race_environment( util::racing::race_pop & race_structure, const problem::base& prob, const std::vector<decision_vector> &x_list1, const std::vector<decision_vector> &x_list2 ) const;
	std::pair<population::size_type, unsigned int> racing__race_for_winner( util::racing::race_pop &race_structure, int idx1, int idx2, unsigned int max_fevals) const;

private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<int &>(m_gen);
		ar & const_cast<double &>(m_omega);
		ar & const_cast<double &>(m_eta1);
		ar & const_cast<double &>(m_eta2);
		ar & const_cast<double &>(m_vcoeff);
		ar & const_cast<int &>(m_variant);
		ar & const_cast<int &>(m_neighb_type);
		ar & const_cast<int &>(m_neighb_param);
		ar & m_fevals;
		ar & const_cast<unsigned int&>(m_max_fevals);
	}
	// Number of generations
	const int m_gen;
	// Particle Inertia weight, or alternatively the constriction coefficient
	const double m_omega;
	// magnitude of the force, applied to the particle's velocity, in the direction of its previous best position
	const double m_eta1;
	// magnitude of the force, applied to the particle's velocity, in the direction of the best position in its neighborhood
	const double m_eta2;
	// Velocity coefficient: velocity values will range in [ - var_range * m_vcoeff, var_range * m_vcoeff ]
	const double m_vcoeff;
	// Velocity update formula
	const int m_variant;
	// Swarm topology
	const int m_neighb_type;
	// parameterization of the swarm topology
	const int m_neighb_param;
	// Parameter controlling how many times a stochastic objective function
	// will be called for averaging
	const unsigned int m_nr_eval_per_x;
	// Incurred objective function evaluation
	mutable unsigned int m_fevals;
	// Maximum allowable fevals before algo terminates
	const unsigned int m_max_fevals;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::pso_generational_racing)

#endif // PAGMO_ALGORITHM_PSO_GENERATIONAL_RACING_H
