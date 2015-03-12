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

// Test code for the transformation of the optima of meta-problems from the original problems.

#include <iostream>
#include <vector>
#include <cmath>
#include "../src/pagmo.h"

using namespace pagmo;

const double EPS = 10e-9;

bool is_eq(const decision_vector & x1, const decision_vector & x2){
	if (x1.size() != x2.size()) return false;
	for (unsigned int i = 0; i < x1.size(); ++i){
		if (fabs(x1[i] - x2[i]) > EPS) return false;
	}
	return true;
}

// Going to test the code on the Himmelblau function,
// see http://en.wikipedia.org/wiki/Himmelblau%27s_function
// The four minima of the function in the range [-6,6]x[-6,6] are:
// (x0,y0) = (3.0, 2.0)
// (x1,y1) = (-2.805118,3.131312)
// (x2,y2) = (-3.779310,-3.283186)
// (x3,y3) = (3.584428,-1.848126)

// Test Himmelblau function minima in the bounding box [-6,6]:
const std::vector<decision_vector> best_x = {
	{3.0, 2.0},
	{-2.805118,3.131312},
	{-3.779310,-3.283186},
	{3.584428,-1.848126},
};

int test_shifted()
{
	problem::himmelblau prob;
	problem::shifted nonzero_shifted_prob(prob, {100,50});
	// retrieve the shifted minima
	const std::vector<decision_vector> shifted_minima = nonzero_shifted_prob.get_best_x();

	// just check that we have actualy moved (somewhere)
	bool all_equal = true;
	for (unsigned int i = 0; i < best_x.size(); ++i) {
		if (!is_eq(best_x[i],shifted_minima[i])) all_equal = false;
	}
	if (all_equal) return 1;


	// Perform a number of shifts that walk around and get back to the origin in the end.
	// It is possible to add any combination of translations satisfying such a condition.
	const std::vector<std::vector<double>> shifts = {
		{1., 1.},
		{-2., 0.},
		{1., -1.},

		{1., 1.},
		{-2., 0.},
		{1., -1.},
	};

	// test the zero shift
	problem::shifted zero_shifted_prob(prob, {0,0});
	const std::vector<decision_vector> zero_shifted_minima = zero_shifted_prob.get_best_x();
	for (unsigned int i = 0; i < best_x.size(); ++i) {
		if (!is_eq(best_x[i],zero_shifted_minima[i])) return 1;
	}

	std::vector<problem::shifted> problems;
	problems.push_back(zero_shifted_prob);

	for (unsigned int i = 0; i < shifts.size(); ++i)
		problems.push_back(problem::shifted(problems.back(), shifts[i]));

	std::vector<std::vector<decision_vector>> expected_shifted_minima;
	expected_shifted_minima.push_back(zero_shifted_prob.get_best_x());

	for (unsigned int i = 0; i < shifts.size(); ++i) {
		expected_shifted_minima.push_back({
				{expected_shifted_minima.back()[0][0] + shifts[i][0], expected_shifted_minima.back()[0][1] + shifts[i][1]},
				{expected_shifted_minima.back()[1][0] + shifts[i][0], expected_shifted_minima.back()[1][1] + shifts[i][1]},
				{expected_shifted_minima.back()[2][0] + shifts[i][0], expected_shifted_minima.back()[2][1] + shifts[i][1]},
				{expected_shifted_minima.back()[3][0] + shifts[i][0], expected_shifted_minima.back()[3][1] + shifts[i][1]},
		});
	}

	for (unsigned int i = 0; i < problems.size(); ++i) {
		const std::vector<decision_vector> shifted_best_x = problems[i].get_best_x();
		for (unsigned int j = 0; j < shifted_best_x.size(); ++j) {
			if (!is_eq(expected_shifted_minima[i][j],shifted_best_x[j])) return 1;
		}
	}

	// Now check that we have actually returned to the origin back again.
	for (unsigned int i = 0; i < best_x.size(); ++i) {
		if (!is_eq(best_x[i],expected_shifted_minima.back()[i])) return 1;
	}

	return 0;
}

int test_rotated()
{
	// tak the original Himmelblau function
	problem::himmelblau prob;
	// and normalize it
	problem::normalized normalized_prob(prob);

	const std::vector<std::vector<double>> rot_id_mat = {{1.,0.},{0.,1.}}; // identity matrix
	problem::rotated rotated_id_prob(normalized_prob,rot_id_mat);
	// Normalized Himmelblau rotated by the identity matrix and the corresponding minima are:
	const std::vector<decision_vector> rot_id_best_x = rotated_id_prob.get_best_x();

	const std::vector<decision_vector> expected_normalized_minima = {
		{best_x[0][0] / 6, best_x[0][1] / 6},
		{best_x[1][0] / 6, best_x[1][1] / 6},
		{best_x[2][0] / 6, best_x[2][1] / 6},
		{best_x[3][0] / 6, best_x[3][1] / 6},
	};

	for (unsigned int i = 0; i < expected_normalized_minima.size(); ++i) {
		if (!is_eq(expected_normalized_minima[i], rot_id_best_x[i])) return 1;
	}

	const double pi = std::acos(-1);
	const double phi = pi / 2;
	const std::vector<std::vector<double>> rot_mat90 = {{std::cos(phi), -std::sin(phi)},{std::sin(phi), std::cos(phi)}};

	// Take the Himmelblau function and rotate by 90 degrees (plus normalize)
	problem::rotated rotated_prob1(prob,rot_mat90);

	// Rotate again by 90 degrees (plus normalize)
	problem::rotated rotated_prob2(rotated_prob1,rot_mat90);

	// Rotate again by 90 degrees (plus normalize)
	problem::rotated rotated_prob3(rotated_prob2,rot_mat90);

	// Rotate again by 90 degrees (plus normalize). We get back to the original problem (normalized_prob),
	// but there is a twist here, because along the way the bounds were normalized several
	// times (exactly four time). Therefore, all the vectors are scaled by the factor (1 / sqrt(2)),
	// except for the frist translation which is normalized by the factor 1 / 6.
	// See src/problem/rotated.cpp for more details (function configure_new_bounds() in particular).
	problem::rotated rotated_prob4(rotated_prob3,rot_mat90);
	// Previous function rotated by 90 degrees counterclockwise and the minima are (after another normalization)
	// basically the same as in the original problem, but normalized by the factor of (1 / sqrt(2)) ** 3:
	const std::vector<decision_vector> rot4_best_x = rotated_prob4.get_best_x();

	const double sqrt2 = std::sqrt(2);
	// First normalization goes from the bounding box of [-6,6]x[-6,6] to [-1,1]x[-1,1]
	// and every other normalization goes from [-sqrt(2),sqrt(2)]x[-sqrt(2),sqrt(2)] to [-1,1]x[-1,1].
	const double coef = (1 / 6.) * (1 / sqrt2) * (1 / sqrt2) * (1 / sqrt2);

	const std::vector<decision_vector> expected_shrunk_minima = {
		{best_x[0][0] * coef, best_x[0][1] * coef},
		{best_x[1][0] * coef, best_x[1][1] * coef},
		{best_x[2][0] * coef, best_x[2][1] * coef},
		{best_x[3][0] * coef, best_x[3][1] * coef},
	};

	for (unsigned int i = 0; i < expected_shrunk_minima.size(); ++i) {
		if (!is_eq(expected_shrunk_minima[i], rot4_best_x[i])) return 1;
	}

	return 0;
}

int test_normalized()
{
	// tak the original Himmelblau function
	problem::himmelblau prob;
	// and normalize it
	problem::normalized normalized_prob(prob);
	// Normalizing bounds (to [-1,1]x[-1,1]) of the Himmelblau function and the corresponding transformed minima are:
	const std::vector<decision_vector> norm_best_x = normalized_prob.get_best_x();

	const std::vector<decision_vector> expected_normalized_minima = {
		{best_x[0][0] / 6, best_x[0][1] / 6},
		{best_x[1][0] / 6, best_x[1][1] / 6},
		{best_x[2][0] / 6, best_x[2][1] / 6},
		{best_x[3][0] / 6, best_x[3][1] / 6},
	};

	for (unsigned int i = 0; i < norm_best_x.size(); ++i) {
		if (!is_eq(expected_normalized_minima[i], norm_best_x[i])) return 1;
	}

	return 0;
}

int main()
{
	return test_shifted() || test_rotated() || test_normalized();
}
