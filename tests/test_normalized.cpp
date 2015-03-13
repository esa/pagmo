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

// Test code for the normalized meta-problem

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cassert>
#include "../src/pagmo.h"
#include "../src/rng.h"
#include <boost/random/uniform_real.hpp>

using namespace pagmo;

const double EPS = 10e-9;

#define PRINT_VEC(x) do {                             \
	std::cout << "[ " #x " ] = ";                     \
	for(unsigned int iii=0; iii<(x).size(); iii++) {  \
		std::cout<<(x)[iii]<<" ";                     \
	}                                                 \
	std::cout<<std::endl; } while(0);

bool is_eq(fitness_vector f1, fitness_vector f2, double eps)
{
	if (f1.size() != f2.size()) {
		return false;
	}
	for (size_t i = 0; i < f1.size(); i++) {
		if (fabs(f1[i] - f2[i]) > eps) {
			return false;
		}
	}
	return true;
}


// TODO: currently this test is not included as it fails: denormalize does not check
// that the input vector is valid.
/// Test that denormalize accepts ONLY decision_vectors with values in [-1, 1].
/**
 * This test is composed of several checks on the behaviour of denormalize:
 *   MUST NOT throw if input has all values in [-1, 1]
 *   MUST throw if at least one value is -1.0 - EPS
 *   MUST NOT throw if one value is -1.0 + EPS
 *   MUST throw if at least one value is +1.0 + EPS
 *   MUST NOT throw if one value is +1.0 - EPS
 */
int test_normalized_denormalize_outside_bounds()
{
	size_t dim = 20;

	// Create problem and decision vector.
	problem::ackley prob((int)dim);
	problem::normalized norm(prob);
	decision_vector p_test(dim);

	rng_double drng	= rng_generator::get<rng_double>();

	std::cout << "Testing denormalize with invalid bounded input: ";

	// Initialise p_test with random values in [-1, 1].
	for (size_t k = 0; k < dim; ++k) {
		p_test[k] = boost::uniform_real<double>(-1.0,1.0)(drng);
	}

	// Check that we indeed all values are accepted.
	try {
		norm.denormalize(p_test);
	} catch (...) {
		std::cout << "Random [-1, 1] denormalize failed." << std::endl;
	}

	// Begin test procedure:
	for (size_t k = 0; k < dim; ++k) {
		double original_value = p_test[k];
		p_test[k] = -1.0 - EPS;
		try {
			// Should throw.
			norm.denormalize(p_test);
			std::cout << "(-1.0 - EPS) value did not cause exception." << std::endl;
			return 1;
		} catch (const value_error& ex) {
			// Expected exception. Do nothing.
		} catch (...) {
			std::cout << "(-1.0 - EPS) value caused unexpected exception." << std::endl;
			return 1;
		}

		p_test[k] = -1.0 + EPS;
		try {
			norm.denormalize(p_test);
		} catch (...) {
			std::cout << "(-1.0 + EPS) value caused exception." << std::endl;
			return 1;
		}

		p_test[k] = 1.0 + EPS;
		try {
			// Should throw.
			norm.denormalize(p_test);
			std::cout << "(1.0 + EPS) value did not cause exception." << std::endl;
			return 1;
		} catch (const value_error& ex) {
			// Expected exception. Do nothing.
		} catch (...) {
			std::cout << "(1.0 + EPS) value caused unexpected exception." << std::endl;
			return 1;
		}

		p_test[k] = 1.0 - EPS;
		try {
			norm.denormalize(p_test);
		} catch (...) {
			std::cout << "1.0 - EPS value caused exception." << std::endl;
			return 1;
		}

		// Restore original value.
		p_test[k] = original_value;
	}
	std::cout << "pass." << std::endl;
	return 0;
}

/// Test various invariants for denormalize.
/**
 * The invariants are as follows:
 *
 * for any k:
 *   normalized(prob).get_lb[k] == -1
 *   normalized(prob).get_ub[k] == 1
 *
 * for any decision_vector dv:
 *   is_eq(normalized(prob).objfun(dv), prob.objfun(normalized(prob).denormalize(dv)))
 *   is_eq(prob.compute_constraints(dv),
 *         normalized(prob).compute_constraints(normalized(prob).denormalize(dv)))
 */
int test_normalized_invariant(
	const std::vector<problem::base_ptr> &probs)
{
	rng_double drng	= rng_generator::get<rng_double>();

	std::cout << "Start batch testing of normalization invariant" << std::endl;

	for (size_t i = 0; i < probs.size(); ++i) {
		problem::base_ptr cur_prob = probs[i]->clone();
		size_t dim = cur_prob->get_dimension();
		decision_vector p_test(dim);

		std::cout<< std::setw(40) << cur_prob->get_name();

		// Normalize problem.
		pagmo::problem::normalized norm(*(cur_prob));

		if (dim != norm.get_lb().size()) {
			std::cout << " Unexpected bounds vector size" << std::endl;
			std::cout << " Expected " << dim << " got " << norm.get_lb().size();
			return 1;
		}

		// Test that we have indeed normalized: lb[i] = -1 and ub[i] = 1.
		for (size_t k = 0; k < dim; ++k) {
			if ((fabs(norm.get_lb()[k] - (-1.0)) > EPS) ||
				(fabs(norm.get_ub()[k] - ( 1.0)) > EPS)) {
				std::cout << " bounds check failed!" << std::endl;
				PRINT_VEC(norm.get_lb());
				PRINT_VEC(norm.get_ub());
				return 1;
			}
		}

		std::cout << " bounds check pass, ";

		// Generate random (normalized) decision vector: -1 <= p_test[i] <= 1
		for (size_t k = 0; k < dim; ++k) {
			p_test[k] = boost::uniform_real<double>(-1.0,1.0)(drng);
		}

		// Check objfun invariant.
		decision_vector p_denormalized = norm.denormalize(p_test);
		fitness_vector f_expected = cur_prob->objfun(p_denormalized);
		fitness_vector f_actual = norm.objfun(p_test);
		if (!is_eq(f_expected, f_actual, EPS)) {
			std::cout << " objfun check failed!" << std::endl;
			PRINT_VEC(f_expected);
			PRINT_VEC(f_actual);
			return 1;
		}

		std::cout << " objfun invariant pass, ";


		// Check compute_constraints invariant.
		constraint_vector c_expected(cur_prob->get_c_dimension());
		constraint_vector c_actual(norm.get_c_dimension());
		cur_prob->compute_constraints(c_expected, p_denormalized);
		norm.compute_constraints(c_actual, p_test);
		if (!is_eq(c_expected, c_actual, EPS)) {
			std::cout << " compute_constraints check failed!" << std::endl;
			PRINT_VEC(c_expected);
			PRINT_VEC(c_actual);
			return 1;
		}

		std::cout << " compute_constraints invariant pass." << std::endl;
	}

	return 0;
}

/// Ensure that when bounds are (-1, 1), normalization doesn't have any effect.
int test_normalized_11(const std::vector<problem::base_ptr> &probs)
{
	rng_double drng = rng_generator::get<rng_double>();

	std::cout << "Start batch testing with bounds (-1, 1)" << std::endl;

	for (size_t i = 0; i < probs.size(); ++i) {
		problem::base_ptr cur_prob = probs[i]->clone();
		size_t dim = cur_prob->get_dimension();
		decision_vector p_normalized(dim);

		std::cout<< std::setw(40) << cur_prob->get_name();

		// Generate a random decision vector with values in (-1, 1).
		for (size_t k = 0; k < dim; ++k) {
			p_normalized[k] = boost::uniform_real<double>(-1.0,1.0)(drng);
		}

		cur_prob->set_bounds(-1.0, 1.0);

		pagmo::problem::normalized prob_norm(*(cur_prob));

		if (dim != prob_norm.get_lb().size()) {
			std::cout << " Unexpected bounds vector size" << std::endl;
			std::cout << " Expected " << dim << " got " << prob_norm.get_lb().size();
			return 1;
		}

		// Ensure normalization did not change the bounds.
		if (!is_eq(prob_norm.get_lb(), cur_prob->get_lb(), EPS) ||
			!is_eq(prob_norm.get_ub(), cur_prob->get_ub(), EPS)) {
			std::cout << " bounds check failed!" << std::endl;
			PRINT_VEC(prob_norm.get_lb());
			PRINT_VEC(cur_prob->get_lb());
			PRINT_VEC(prob_norm.get_ub());
			PRINT_VEC(cur_prob->get_ub());
			return 1;
		}

		std::cout << " bounds check pass, ";

		// Ensure denormalize() is essentially a no-op.
		decision_vector p_denormalized = prob_norm.denormalize(p_normalized);
		if (!is_eq(p_denormalized, p_normalized, EPS)) {
			std::cout << " denormalize failed!" << std::endl;
			PRINT_VEC(p_denormalized);
			PRINT_VEC(p_normalized);
			return 1;
		}

		std::cout << " denormalize pass." << std::endl;
	}
	return 0;
}


int main()
{
	int dimension = 40;
	std::vector<problem::base_ptr> probs;

	// ZDT and DTLZ are Box-Constrained Continuous Multi-Objective
	for (int i = 1; i <= 6; ++i) {
		probs.push_back(problem::zdt(i, dimension).clone());
	}
	for (int i = 1; i <= 7; ++i) {
		probs.push_back(problem::dtlz(i, dimension).clone());
	}

	// Box-constrained Continuous Single-Objective
	probs.push_back(problem::ackley(dimension).clone());
	probs.push_back(problem::rastrigin(dimension).clone());

	// CEC2006 are Constrained Continuous Single-Objective
	for (int i = 1; i <= 24; ++i) {
		probs.push_back(problem::cec2006(i).clone());
	}

	// CEC2009 problems 1-10 are Constrained Continuous Multi-Objective
	for (int i = 1; i <= 10; ++i) {
		probs.push_back(problem::cec2009(i, dimension, true).clone());
	}

	// Add a few problems with extreme bounds.
	problem::ackley ak(10);
	ak.set_bounds(DBL_MIN, DBL_MAX);
	probs.push_back(ak.clone());
	ak.set_bounds((double)INT_MIN, (double)INT_MAX);
	probs.push_back(ak.clone());
	ak.set_bounds((double)INT_MIN, (double)(INT_MIN + 1));
	probs.push_back(ak.clone());
	ak.set_bounds(1.0 - DBL_EPSILON, 1.0 + DBL_EPSILON);
	probs.push_back(ak.clone());
	ak.set_bounds(0.0 - DBL_EPSILON, 0.0 + DBL_EPSILON);
	probs.push_back(ak.clone());

	return test_normalized_11(probs) ||
		   test_normalized_invariant(probs);
}
