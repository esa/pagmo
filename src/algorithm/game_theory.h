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

#ifndef PAGMO_ALGORITHM_GAME_THEORY_H
#define PAGMO_ALGORITHM_GAME_THEORY_H

#include "../config.h"
#include "../serialization.h"
#include "base.h"
#include "jde.h"
#include "../problem/decompose.h"
#include "../population.h"
#include "../archipelago.h"

namespace pagmo { namespace algorithm {

//! Typedef for decomposition weights
typedef std::vector< double > weights_type;

//! Typedef for vector of decomposition weights
typedef std::vector< weights_type > weights_vector_type;

//! Typedef for population pointer
typedef boost::shared_ptr< population > pop_ptr;

/// Game Theory
/**
 *
 * TODOThis class implement a multi-objective optimization algorithm based on parallel decomposition.
 * For each element of the population a different single objective problem is generated using
 * a decomposition method. Those single-objective problems are thus solved in parallel.
 * At the end of the evolution the population is set as the best individual for each single-objective problem.
 *
 * Game Theory assumes all the objectives need to be minimized.
 *
 * @author Jacco Geul (jacco.geul@gmail.com)
 **/

class __PAGMO_VISIBLE game_theory: public base
{
  public:
	/// Mechanism used to generate the weight vectors
	enum weight_generation_type {
		UNIFORM,
		RANDOM,
		TCHEBYCHEFF,
		TCHEBYCHEFF_ADAPTIVE,
		ADAPTIVE
	};
	/// Different types of downscaling methods
	enum downscaling_type {
		NOSCALING,
		SINGULARSCALING,
		THRESHOLDSCALING,
		RANDOMSCALING
	};
	game_theory(
	        int = 10,
		unsigned int = 0,
		unsigned int = 1,
		const pagmo::algorithm::base & = pagmo::algorithm::jde(8),
		const weights_vector_type & = weights_vector_type(),
		const weights_vector_type & = weights_vector_type(),
		weight_generation_type = ADAPTIVE,
		downscaling_type = THRESHOLDSCALING,
		const std::vector< double > & = std::vector< double >(1,1e-6),
		const std::vector< double > & = std::vector< double >(1,1e-6)
	);
	game_theory(const game_theory &);

	base_ptr clone() const;
	void evolve(population &) const;
	std::string get_name() const;
	weights_vector_type generate_weights(const unsigned int, const unsigned int, const bool, const bool ) const;
	weights_vector_type get_var_weights() const;
	weights_vector_type get_obj_weights() const;
  protected:
	std::string human_readable_extra() const;

  private:
	friend class boost::serialization::access;
	template <class Archive>
		void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<int &>(m_gen);
		ar & const_cast<unsigned int &>(m_dim);
		ar & const_cast<unsigned int &>(m_init_dim);
		ar & const_cast<int &>(m_fevals);
		ar & const_cast<int &>(m_cevals);
		ar & m_arch;
		ar & m_z;
		ar & m_best_vector;
		ar & const_cast<pop_ptr &>(m_pop);
		ar & const_cast<unsigned int &>(m_threads);
		ar & const_cast<base_ptr &>(m_solver);
		ar & m_var_weights;
		ar & m_obj_weights;
		ar & m_init_var_weights;
		ar & m_init_obj_weights;
		ar & const_cast<weight_generation_type &>(m_weight_generation);
		ar & const_cast<downscaling_type &>(m_downscaling);
		ar & m_relative_tolerance;
		ar & m_absolute_tolerance;
	}
	//Number of generations
	const int m_gen;
	mutable unsigned int m_dim;
	mutable unsigned int m_init_dim;
	mutable int m_fevals;
	mutable int m_cevals;
	mutable archipelago m_arch;
	mutable fitness_vector m_z;
	mutable std::vector< fitness_vector > m_best_vector;
	mutable pop_ptr m_pop;
	const unsigned int m_threads;
	const base_ptr m_solver;
	mutable weights_vector_type m_var_weights;
	mutable weights_vector_type m_obj_weights;
	mutable weights_vector_type m_init_var_weights;
	mutable weights_vector_type m_init_obj_weights;
	mutable weight_generation_type m_weight_generation;
	mutable downscaling_type m_downscaling;
	std::vector< double > m_relative_tolerance;
	std::vector< double > m_absolute_tolerance;
	
	// Private functions
	void initialise_weights( const problem::base::f_size_type, 
		const problem::base::size_type ) const;
	void decompose( population ) const;
	void update_population( ) const;
	fitness_vector compute_best_vector( ) const;
	void adapt_obj_weights( ) const;
	void downscale( ) const;
  	template <typename T>
		std::vector<T> sum_of_vec(const std::vector<T>& a, const std::vector<T>& b) const;
	template <typename T>
		std::vector<T> had_of_vec(const std::vector<T>& a, const std::vector<T>& b) const;
	template <typename T>
		std::vector<T> inv_of_vec(const std::vector<T>& a) const;
	template <typename T>
		bool solution_within_tolerance(const std::vector<T>& a, const std::vector<T>& b) const;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::game_theory)

#endif // PAGMO_ALGORITHM_GAME_THEORY_H
