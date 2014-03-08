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


#include <string>
#include <vector>
#include <algorithm>

#include <boost/math/special_functions/binomial.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../archipelago.h"
#include "../island.h"
#include "../population.h"
#include "../topology/unconnected.h"
#include "../problem/decompose.h"
#include "../util/discrepancy.h"
#include "../util/neighbourhood.h"
#include "../types.h"
#include "base.h"
#include "game_theory.h"

namespace pagmo { namespace algorithm {

/// Constructor
 /**
 * Constructs a Game Theory algorithm
 *
 * @param[in] gen Number of generations to evolve.
 * @param[in] threads the amounts of threads that will be used
 * @param[in] solver the algorithm to solve the single objective problems.
 * @param[in] var_weights the weights for the decision variables.
 * @param[in] obj_weights the decomposition weights for the objective functions.
 * @param[in] weight_generation type of weight generation.
 * @param[in] relative_tolerance for determining convergence.
 * @param[in] absolute_tolerance for determining convergence.
 * @throws value_error if gen is negative, weight_generation is not sane.
 */
game_theory::game_theory(int gen,
	unsigned int dim,
	unsigned int threads,
	const pagmo::algorithm::base & solver,
	const weights_vector_type &var_weights,
	const weights_vector_type &obj_weights,
	weight_generation_type weight_generation,
	const std::vector< double > &relative_tolerance,
	const std::vector< double > &absolute_tolerance )
	:base(),
	 m_gen(gen),
	 m_dim(dim),
	 m_threads(threads),
	 m_solver(solver.clone()),
	 m_var_weights(var_weights),
	 m_obj_weights(obj_weights),
	 m_weight_generation(weight_generation),
	 m_relative_tolerance(relative_tolerance),
	 m_absolute_tolerance(absolute_tolerance)
{
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if(m_weight_generation != UNIFORM && m_weight_generation != RANDOM && 
		m_weight_generation != TCHEBYCHEFF && m_weight_generation != ADAPTIVE && 
	m_weight_generation != TCHEBYCHEFF_ADAPTIVE ) {
		pagmo_throw(value_error,"non existing weight generation method.");
	}
}

/// Copy constructor. Performs a deep copy. Necessary as a pointer to a base algorithm is here contained
game_theory::game_theory(const game_theory &algo):
	base(algo),
	m_gen(algo.m_gen),
	m_dim(algo.m_dim),
	m_threads(algo.m_threads),
	m_solver(algo.m_solver->clone()),
	m_var_weights(algo.m_var_weights),
	m_obj_weights(algo.m_obj_weights),
	m_weight_generation(algo.m_weight_generation),
	m_relative_tolerance(algo.m_relative_tolerance),
	m_absolute_tolerance(algo.m_absolute_tolerance)
{}

/// Clone method.
base_ptr game_theory::clone() const
{
	return base_ptr(new game_theory(*this));
}

// Sum of two vectors
template <typename T>
std::vector<T> game_theory::sum_of_vec(const std::vector<T>& a, const std::vector<T>& b) const
{
	assert(a.size() == b.size());

	std::vector<T> result;
	result.reserve(a.size());

	std::transform(a.begin(), a.end(), b.begin(), 
		std::back_inserter(result), std::plus<T>());
	return result;
}

// Hadamard product of two vectors
template <typename T>
std::vector<T> game_theory::had_of_vec(const std::vector<T>& a, const std::vector<T>& b) const
{
	assert(a.size() == b.size());

	std::vector<T> result;
	result.reserve(a.size());

	std::transform(a.begin(), a.end(), b.begin(), 
		std::back_inserter(result), std::multiplies<T>());
	return result;
}

// Inverse of a vector
template <typename T>
std::vector<T> game_theory::inv_of_vec(const std::vector<T>& a) const
{
	std::vector<T> b;
	for( int i = 0; i < a.size(); ++i )
		b.push_back(fabs( a[i] - 1 ));
	return b;
}

// Absolute difference between two vectors
template <typename T>
bool game_theory::solution_within_tolerance(const std::vector<T>& a, const std::vector<T>& b) const
{
	// Check if sizes match
	assert(m_relative_tolerance.size() == m_absolute_tolerance.size());
	assert(a.size() == b.size());

	// Create temporary vectors
	std::vector< double > rel_tol;
	std::vector< double > abs_tol;

	// Check if single tolerance value is used or one for each
	// variable.
	if( m_relative_tolerance.size() == 1){
		rel_tol = std::vector< double >( a.size(), m_relative_tolerance[0] );
		abs_tol = std::vector< double >( a.size(), m_absolute_tolerance[0] );
	} else {
		assert(a.size() == m_absolute_tolerance.size());
		rel_tol = m_relative_tolerance;
		abs_tol = m_absolute_tolerance;
	}
		
	// Check if each elements meets tolerance
	bool withintol = true;
	for( unsigned int i = 0; i < a.size(); ++i){
		withintol = withintol * ( fabs(a[i] - b[i]) < 
			std::max( abs_tol[i], 
				std::max(fabs(a[i]),fabs(b[i])) * rel_tol[i] ));
	}
	return withintol;
}

/// Generates the weights 
/**
 * Generates the weights used in the problem decomposition.
 *
 * @param[in] n_x length of weight vector.
 * @param[in] n_v length of vector of weight vectors.
 * @param[in] fracs if true the sum of a weight vector is always one.
 * @param[in] fracs if true the sum of a weight vector is always one.
 * @param[out] vector of weight vectors
 */
weights_vector_type game_theory::generate_weights(const unsigned int n_x, 
	const unsigned int n_v, const bool fracs, const bool random ) const
{
	// Preform sanity checks
	if ( n_v > n_x ) {
		pagmo_throw(value_error, "The number of weight vector cannot be longer than the number of entries in a weight vector.");
	}

	std::vector< int >  c( n_v, 0);
	std::vector< int >  a( n_x, 0);
	weights_vector_type r( n_v, weights_type( n_x, 0.0 ));
	int shift = 0;

	// Add a shift to prevent bias towards one number if n_x / n_v
	// is not natural
	if( random )
		shift = static_cast< int >( m_drng() * n_v);

	for( unsigned int i = 0; i < n_x; ++i ){
		// Get idx and assign
		int idx = (i + shift) % n_v;
		a[i] = idx;

		// Count the number of occurrences for fraction
		c[idx]++;
	}

	// Shuffle around
	if( random )
		std::random_shuffle(a.begin(), a.end());

	// If no fracs then reset all counters to 1
	if(!fracs)
		c = std::vector< int >( n_v, 1 );

	for( unsigned int i = 0; i < n_x; ++i ){
		r[a[i]][i] = 1.0 / c[a[i]];
	}
	return r;
}

/// Downscale the decomposition
/**
 * 
 */
void game_theory::downscale( ) const
{

	weights_vector_type new_var_weights;
	weights_vector_type new_obj_weights;

	// Option not implemented yet, but allows the selection based 
	// on the current weight -> 0.1 has a 10 percent chance of 
	// being activated.
	bool weighted_random = false;

	for( unsigned int i =  0; i < m_dim; ++i ){

		// Make new weight vector of either activated (1) or 
		// inactivated weights (0), based on threshold
		weights_type new_obj_vector( m_obj_weights[0].size(), 0.0 );
		double threshold = 0.1 / m_obj_weights[0].size();
		int hitcounter = 0;
		for( unsigned int j = 0; j < m_obj_weights[0].size(); ++j ){
			if( m_obj_weights[i][j] > threshold || 
				( static_cast< int >( m_drng() + m_obj_weights[i][j] ) == 1 && weighted_random)){
				new_obj_vector[j] = 1.0;
				hitcounter++;
			}
		}
		// Normalise so sum = 1.0
		for( unsigned int j = 0; j < m_obj_weights[0].size(); ++j ){
			new_obj_vector[j] /= static_cast< double >( hitcounter );
		}


		// Find the matching entry in new_obj_weights (if any)
		bool match_found = false;
		for( unsigned int j =  0; j < new_obj_weights.size(); ++j ){
			if( new_obj_vector == new_obj_weights[j] ){
				// Then add the var_vector to the current
				for( unsigned int k = 0; k < m_var_weights[0].size(); ++k ){
					if( m_var_weights[i][k] == 1.0 ){
						new_var_weights[j][k] = 1.0;
					}
				}
				match_found = true;
				break;
			}
		}
		// Otherwise add new entries to the obj and var weights.
		if( !match_found ){
			new_obj_weights.push_back( new_obj_vector );
			new_var_weights.push_back( m_var_weights[i] );
		}
	}
	m_dim = new_obj_weights.size();
	m_obj_weights = new_obj_weights;
	m_var_weights = new_var_weights;
	m_weight_generation = UNIFORM;
}

/// Evolve implementation.
/**
 * Run the Game Theory algorithm for the number of generations specified in the constructors.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */
void game_theory::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const population::size_type pop_size = pop.size();
	const problem::base::f_size_type prob_objectives = prob.get_f_dimension();
	const problem::base::size_type prob_dimension = prob.get_dimension();

	// Preform sanity checks
	if ( prob_objectives < 2 ) {
		pagmo_throw(value_error, "The problem is not multiobjective, try some other algorithm than Game Theory");
	}

	if ( prob_dimension < 2 ) {
		pagmo_throw(value_error, "The problem has only one decision variable. This is not supported by Game Theory, select a different algorithm.");
	}

	// Number of populations.
	unsigned int min_dim = std::min( prob_objectives, prob_dimension );
	if ( m_dim == 0 ){
		m_dim = prob_dimension;
	}
	if ( m_dim > prob_dimension ) {
		pagmo_throw(value_error, "The dimension of the decomposition can not be greater than the problem dimension.");
	}
	if ( min_dim > m_dim ) {
		pagmo_throw(value_error, "The dimension of the decomposition should be at least the minimum between the problem dimension and the number of objectives.");
	}
	
	// Get out if there is nothing to do.
	if (m_gen == 0) {
		return;
	}

	bool random_weights = false;
	bool adaptive_weights = false;

	problem::decompose::method_type decompose_method = pagmo::problem::decompose::WEIGHTED;

	switch (m_weight_generation)
	{
		case UNIFORM : 
			decompose_method = pagmo::problem::decompose::WEIGHTED;
		        random_weights   = false;
			adaptive_weights = false;
			break;
		case RANDOM :
			// Ignoring previously calculated weight matrices
			m_var_weights.clear();
			m_obj_weights.clear();
			decompose_method = pagmo::problem::decompose::WEIGHTED;
			random_weights   = true;
			adaptive_weights = false;
			break;
		case TCHEBYCHEFF : 
			decompose_method = pagmo::problem::decompose::TCHEBYCHEFF;
			random_weights   = true;
			adaptive_weights = false;
			m_obj_weights = generate_weights( prob_objectives, 1, true, random_weights );
			m_obj_weights = weights_vector_type( m_dim, m_obj_weights[0] );
			break;
		case TCHEBYCHEFF_ADAPTIVE : 
			decompose_method = pagmo::problem::decompose::TCHEBYCHEFF;
			random_weights   = false;
			adaptive_weights = true;
			break;
		case ADAPTIVE : 
			decompose_method = pagmo::problem::decompose::WEIGHTED;
			random_weights   = false;
			adaptive_weights = true;
			break;
	}

	if(m_var_weights.empty()){
		m_var_weights = generate_weights( prob_dimension, m_dim, false, random_weights );
	}

	if( adaptive_weights ){
		// If adaptive every objective should be coupled to
		// every decomposed problem initially.
		m_obj_weights = generate_weights( prob_objectives, 1, true, random_weights );
		m_obj_weights = weights_vector_type( m_dim, m_obj_weights[0] );
	}else if(m_obj_weights.empty()){
		// Generate the objective weights
		weights_vector_type tmp_weights = generate_weights( 
			prob_objectives, min_dim, true, random_weights );
		m_obj_weights = tmp_weights;

		// If m_dim > min_dim weights have to be repeated
		// (i.e. one objective linked to more than one
		// decomposed population).
		int k;
		for( unsigned int i = 0; i < (m_dim - min_dim); ++i ){
			// If m_dim > 2 * min_dim, reset index
			k = i % min_dim;
			// If m_dim > min_dim all potential vector are
			// already in tmp_weights, just reshuffle.
			if( k == 0 && random_weights ){
				std::random_shuffle(tmp_weights.begin(), tmp_weights.end());
			}
			m_obj_weights.push_back( tmp_weights[ k ] );
		}
	}

	// More sanity checks
	if ( m_var_weights.size() != m_dim ) {
		pagmo_throw(value_error, "The vector of variable weights has an incorrect number of entries. Create an entry for each population.");
	}
	if ( m_obj_weights.size() != m_dim ) {
		pagmo_throw(value_error, "The vector of objective weights has an incorrect number of entries. Create an entry for each population.");
	}
	if ( m_var_weights[0].size() != prob_dimension ) {
		pagmo_throw(value_error, "The dimension of the variable weights do not match the problem. The dimension must be equal to the number of decision variables.");
	}
	if ( m_obj_weights[0].size() != prob_objectives ) {
		pagmo_throw(value_error, "The dimension of the objective weights do not match the problem. The dimension must be equal to the number of objectives.");
	}

	// Compute the starting ideal point
	fitness_vector ideal_point = pop.compute_ideal();

	// Create all the decomposed problems 
	std::vector<pagmo::problem::decompose*> prob_vector;
	for( problem::base::f_size_type i = 0; i < m_dim; i++ ) {
		prob_vector.push_back(
			new pagmo::problem::decompose( prob, decompose_method, 
				m_obj_weights[i], ideal_point, true ));
	}

	// Set unconnected topology, only using arch for parallel processing
	topology::unconnected topo;
	
	// Create unconnected archipelago of m_dim islands. Each island solve a different decomposed part of the problem.
	pagmo::archipelago arch(topo);

	// Sets random number generators of the archipelago using the
	// algorithm urng to obtain a deterministic behaviour upon
	// copy.
	arch.set_seeds(m_urng());

	// Assign population to each problem
	for( problem::base::f_size_type i = 0; i < m_dim; i++ ) { 
		
		// Create an empty population for each decomposed
		// problem
		population decomp_pop(*prob_vector[i], 0, m_urng()); 

		// Fill the population from the original and calculate
		// the objective.
		for ( population::size_type j = 0; j<pop.size(); j++ ) {
			decomp_pop.push_back( pop.get_individual(j).cur_x );
		}

		// Put pop in island and add to arch.
		arch.push_back(pagmo::island(*m_solver, decomp_pop, 1.0));
	}

	// Create a last best vector for convergence check
	std::vector< double > last_best_vector( prob_dimension, 0.0 );

	for(int g = 0; g < m_gen; ++g) {

		// Define best decision vector for fixed variables
		std::vector< double > best_vector( prob_dimension, 0.0 );
		for( problem::base::f_size_type i = 0; i < m_dim; i++ ) {
			int best_idx = arch.get_island(i)->get_population().get_best_idx();
			best_vector = sum_of_vec( best_vector,  
				had_of_vec( m_var_weights[i], 
					arch.get_island(i)->get_population().get_individual( best_idx ).cur_x ));
		}

		// Check if Nash equilibrium is reached
		
		if( solution_within_tolerance( best_vector, last_best_vector )){
			break;
		}

		last_best_vector = best_vector;

		// Change the boundaries of each problem, this is
		// equal to fixing!
		for( problem::base::f_size_type i = 0; i < m_dim; i++ ) {
			// Get problem pointer from population from
			// island. Each get creates a clone (except
			// problem_ptr), so that they must be set back
			// later on.
			base_island_ptr isl_i = arch.get_island(i);
			population pop_i = isl_i->get_population();

			// This doesn't work as set_bounds is not a
			// constant function as the bounds itself are
			// not mutable.
			// problem::base &prob_i = pop_i.problem();

			// Need to recast on problem::base_ptr doesn't
			// have all functions of problem::decompose
			boost::shared_ptr< problem::decompose > prob_i = 
				boost::static_pointer_cast< problem::decompose >(
					population_access::get_problem_ptr(pop_i));

			// Inverse of the decision variable weight. 
			weights_type inverse_weight = inv_of_vec( m_var_weights[i] );

			// Calculate modified lower and upper bounds.
			std::vector< double > mod_lb = sum_of_vec(
				had_of_vec( inverse_weight, best_vector ),
				had_of_vec( m_var_weights[i], prob_i->get_lb() ));
			std::vector< double > mod_ub = sum_of_vec(
				had_of_vec( inverse_weight, best_vector ),
				had_of_vec( m_var_weights[i], prob_i->get_ub() ));

			// Change the bounds of the problem.
			prob_i->set_bounds( mod_lb, mod_ub );

			// In case of adaptive weights use the f
			// minmax history to update the weights.
			if( adaptive_weights ){
				std::vector< fitness_vector > f_history = prob_i->get_minmax_history();
				prob_i->reset_minmax_history();
				double obj_sum = 0.0;
				for( unsigned int j = 0; j < prob_objectives; ++j ){
					m_obj_weights[i][j] = f_history[1][j] - f_history[0][j];
					obj_sum += m_obj_weights[i][j];
				}
				for( unsigned int j = 0; j < prob_objectives; ++j ){
					m_obj_weights[i][j] /= obj_sum;
				}
				prob_i->set_weights( m_obj_weights[i] );
			}
			bool reinit_pop_after_bounds_change = true;

			// Change the chromosomes according to bounds.
			for( population::size_type j = 0; j < pop_i.size() && 
				      reinit_pop_after_bounds_change; j++ ) {
				pop_i.set_x(j, sum_of_vec(
						had_of_vec( inverse_weight, best_vector ),
						had_of_vec( m_var_weights[i], pop_i.get_individual(j).cur_x )));
			}

			// Update global ideal point.
			fitness_vector z_i = prob_i->get_ideal_point();
			for( unsigned int j = 0; j < z_i.size(); ++j ){
				if( z_i[j] < ideal_point[j] ){
					ideal_point[j] = z_i[j];
				}
			}
			
			// Update local ideal points. This causes a
			// lag between synchronising of one
			// generation, but it saves pulling everything
			// inside-out again.
			prob_i->set_ideal_point( ideal_point );
			
			// Set population back into island, island
			// back into arch.
			isl_i->set_population( pop_i );
			arch.set_island( i, *isl_i );
		}

		// Evolve entire archipelago once
		arch.evolve_batch(1, m_threads);
	}

	// Calculate best vector for last time and push to population.
	std::vector< double > best_vector( prob_dimension, 0.0 );
	for( problem::base::f_size_type i = 0; i < m_dim; i++ ) {
		best_vector = sum_of_vec( best_vector,  
			had_of_vec( m_var_weights[i], 
				arch.get_island(i)->get_population().champion().x ));
	}
	pop.push_back( best_vector );

        // Fill population will all individuals.
	for (problem::base::f_size_type i=0; i < m_dim; ++i) {
		for (population::size_type j=0; j < pop_size; ++j) {
			pop.push_back(arch.get_island(i)->get_population().get_individual(j).cur_x);
		}
		// Add the number of fevals
		prob.add_fevals( arch.get_island(i)->get_population().problem().get_fevals() );
		prob.add_cevals( arch.get_island(i)->get_population().problem().get_cevals() );
	}

	// Get sorted list from best to worst
	std::vector<population::size_type> best_idx = pop.get_best_idx( pop.size() );

	// Get worst pop.size() - NP
	std::vector<population::size_type> worst_idx( best_idx.begin() + pop_size, best_idx.end() );

	// Sort list in descending order
	std::sort( worst_idx.begin(), worst_idx.end(), std::greater<int>());

	// Remove worst from population (back to front)
	for (population::size_type i = 0; i < worst_idx.size(); ++i) {
		pop.erase( worst_idx[i] );
	}
	
}

/// Algorithm name
std::string game_theory::get_name() const
{
	return m_solver->get_name() + "[Game Theory]";
}

/// Get var_weights.
/**
 * Will return the decision variable weights.
 *
 * @param[out] the vector of weight vectors.
 */
weights_vector_type game_theory::get_var_weights() const {
	return m_var_weights;
}

/// Get obj_weights.
/**
 * Will return the objective variable weights.
 *
 * @param[out] the vector of weight vectors.
 */
weights_vector_type game_theory::get_obj_weights() const {
	return m_obj_weights;
}

/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string game_theory::human_readable_extra() const
{
	std::ostringstream s;
	s << "gen:" << m_gen << ' ';
	s << "threads:" << m_threads << ' ';
	s << "solver:" << m_solver->get_name() << ' ';
	s << "decomposition: WEIGHTED" << ' ';
	s << "weights:";
	switch (m_weight_generation)
	{
		case UNIFORM : 
			s << "UNIFORM" << ' ';
			break;
		case RANDOM : 
			s << "RANDOM" << ' ';
			break;
		case TCHEBYCHEFF : 
			s << "TCHEBYCHEFF" << ' ';
			break;
		case TCHEBYCHEFF_ADAPTIVE : 
			s << "TCHEBYCHEFF_ADAPTIVE" << ' ';
			break;
		case ADAPTIVE : 
			s << "ADAPTIVE" << ' ';
			break;
	}
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::game_theory)
