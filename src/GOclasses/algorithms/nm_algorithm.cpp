/*****************************************************************************
 *   Copyright (C) 2008, 2009 Advanced Concepts Team (European Space Agency) *
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

// 03/02/2009: Initial version by Francesco Biscani.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include "../../exceptions.h"
#include "../basic/individual.h"
#include "../basic/population.h"
#include "../problems/GOproblem.h"
#include "nm_algorithm.h"

nm_algorithm::nm_algorithm(int n_gen, const double &alpha, const double &gamma, const double &rho, const double &sigma):
	m_gen(n_gen),m_alpha(alpha),m_gamma(gamma),m_rho(rho),m_sigma(sigma)
{
	if (n_gen <= 0) {
		pagmo_throw(value_error,"number of generations must be positive");
	}
	if (alpha <= 0) {
		pagmo_throw(value_error,"reflection coefficient must be positive");
	}
	if (gamma <= 1) {
		pagmo_throw(value_error,"expansion coefficient must be greater than one");
	}
	if (rho <= 0 || rho >= 1) {
		pagmo_throw(value_error,"contraction coefficient must be in the ]0,1[ range");
	}
	if (sigma <= 0 || sigma >= 1) {
		pagmo_throw(value_error,"shrink coefficient must be in the ]0,1[ range");
	}
}

struct sorter {
	sorter(const GOProblem &p):problem(p) {}
	bool operator()(const std::pair<std::vector<double>,double> &a, const std::pair<std::vector<double>,double> &b) const {
		return (a.second < b.second);
	}
	const GOProblem &problem;
};

void nm_algorithm::check_bounds(std::vector<double> &v, const GOProblem &p) const
{
	const size_t size = v.size();
	const std::vector<double> &LB = p.getLB(), &UB = p.getUB();
	pagmo_assert(size == LB.size() && size == UB.size());
	for (size_t i = 0; i < size; ++i) {
		if (v[i] < LB[i]) {
			v[i] = LB[i];
		}
		if (v[i] > UB[i]) {
			v[i] = UB[i];
		}
	}
}

std::vector<double> nm_algorithm::center_mass(const simplex &s) const
{
	const size_t size = s.size() - 1;
	std::vector<double> retval(size);
	for (size_t i = 0; i < size; ++i) {
		double tmp = 0;
		for (size_t j = 0; j < size; ++j) {
			pagmo_assert(size == s[j].first.size());
			tmp += s[j].first[i];
		}
		retval[i] = tmp / size;
	}
	return retval;
}

// Compute d + c * (a - b).
std::vector<double> nm_algorithm::sub_mult_add(const std::vector<double> &a, const std::vector<double> &b,
	const double &c, const std::vector<double> &d) const
{
	const size_t size = a.size();
	pagmo_assert(b.size() == size && d.size() == size);
	std::vector<double> retval(size);
	for (size_t i = 0; i < size; ++i) {
		retval[i] = d[i] + c * (a[i] - b[i]);
	}
	return retval;
}

double nm_algorithm::simplex_diameter(const simplex &s) const
{
	const size_t s_size = s.size(), dim = s_size - 1;
	pagmo_assert(s_size > 1);
	double retval = 0;
	// Compute all vertex distances and keep the highest.
	for (size_t i = 0; i < s_size; ++i) {
		for (size_t j = i + 1; j < s_size; ++j) {
			double tmp = 0.;
			for (size_t k = 0; k < dim; ++k) {
				tmp += (s[i].first[k] - s[j].first[k]) * (s[i].first[k] - s[j].first[k]);
			}
			tmp = std::sqrt(tmp);
			if (tmp > retval) {
				retval = tmp;
			}
		}
	}
	return retval;
}

// Shuffle the coordinates of the simplex, excluding the first vertex.
void nm_algorithm::shuffle_simplex(simplex &s, const GOProblem &problem) const
{
	const size_t s_size = s.size();
	pagmo_assert(s_size > 1);
	std::vector<std::vector<double> > simplex_coordinates(s_size - 1, std::vector<double>(s_size - 1));
	for (size_t i = 0; i < s_size - 1; ++i) {
		for (size_t j = 0; j < s_size - 1; ++j) {
			simplex_coordinates[i][j] = s[j + 1].first[i];
		}
	}
	for (size_t i = 0; i < s_size - 1; ++i) {
		std::random_shuffle(simplex_coordinates[i].begin(),simplex_coordinates[i].end());
	}
	for (size_t i = 0; i < s_size - 1; ++i) {
		for (size_t j = 0; j < s_size - 1; ++j) {
			s[j + 1].first[i] = simplex_coordinates[i][j];
			s[j + 1].second = problem.objfun(s[j + 1].first);
		}
	}
}

Population nm_algorithm::evolve(const Population &pop) const
{
	// Preliminary checks and useful variables.
	const GOProblem &problem = pop.problem();
	const size_t pop_size = pop.size(), prob_size = problem.getDimension();
	if (pop_size <= prob_size) {
		pagmo_throw(value_error,"too few individuals in population for Nelder-Mead method");
	}
	if (prob_size < 1) {
		pagmo_throw(value_error,"the dimension of the problem must be at least 1 for Nelder-Mead method");
	}
	Population retval(pop);
	retval.sort();
	// Build a simplex from the sorted input population.
	const size_t simplex_size = prob_size + 1;
	simplex s(simplex_size,std::make_pair(std::vector<double>(prob_size),0));
	for (size_t i = 0; i < simplex_size; ++i) {
		s[i].first = retval[i].getDecisionVector();
		s[i].second = retval[i].getFitness();
	}
	// Perform the NM method for a number of times equal to m_gen.
	for (size_t gen = 0; gen < m_gen; ++gen) {
		const double diameter = simplex_diameter(s);
		// NOTE: hard coded value, is it worth to make it configurable?
		if (diameter < 1E-6) {
			//std::cout << "small diameter: " << diameter << '\n';
			for (size_t i = 1; i < simplex_size; ++i) {
				Individual tmp(problem);
				s[i].first = tmp.getDecisionVector();
				s[i].second = tmp.getFitness();
			}
		}
		// First order the vertices of the simplex according to fitness.
		std::sort(s.begin(),s.end(),sorter(problem));
		// With small probability mix the coordinates between the vertices (excluding the best one).
		if (drng() < .001) {
			shuffle_simplex(s,problem);
		}
		// Compute the center of mass, excluding the worst point.
		const std::vector<double> x0 = center_mass(s);
		// Compute a reflection.
		std::vector<double> xr = sub_mult_add(x0,s.back().first,m_alpha,x0);
		check_bounds(xr,problem);
		// Calculate fitness values. fitness_target is relative to the vertex immediately before the worst one.
		const double fitness_r = problem.objfun(xr), fitness_best = problem.objfun(s[0].first), fitness_target = problem.objfun(s[simplex_size - 2].first);
		if (fitness_r >= fitness_best && fitness_r < fitness_target) {
			s.back().first = xr;
			s.back().second = fitness_r;
		} else if (fitness_r < fitness_best) {
			std::vector<double> xe = sub_mult_add(x0,s.back().first,m_gamma,x0);
			check_bounds(xe,problem);
			const double fitness_e = problem.objfun(xe);
			if (fitness_e < fitness_r) {
				s.back().first = xe;
				s.back().second = fitness_e;
			} else {
				s.back().first = xr;
				s.back().second = fitness_r;
			}
		} else {
			std::vector<double> xc = sub_mult_add(x0,s.back().first,m_rho,s.back().first);
			check_bounds(xc,problem);
			const double fitness_c = problem.objfun(xc), fitness_worst = problem.objfun(s.back().first);
			if (fitness_c <= fitness_worst) {
				s.back().first = xc;
				s.back().second = fitness_c;
			} else {
				for (size_t i = 1; i < simplex_size; ++i) {
					s[i].first = sub_mult_add(s[i].first,s.front().first,m_sigma,s.front().first);
					check_bounds(s[i].first,problem);
					s[i].second = problem.objfun(s[i].first);
				}
			}
		}
	}
	// In retval overwrite the individuals that have evolved (i.e., those involved in the simplex).
	for (size_t i = 0; i < simplex_size; ++i) {
		retval[i] = Individual(s[i].first,retval[i].getVelocity(),s[i].second);
	}
	return retval;
}

void nm_algorithm::log(std::ostream &s) const
{
	s << "NM - generations:" << m_gen << " alpha:" << m_alpha << " gamma:" << m_gamma
		<< " rho:" << m_rho << " sigma:" << m_sigma;
}
