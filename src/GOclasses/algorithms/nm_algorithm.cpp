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
	bool operator()(const std::vector<double> &v1, const std::vector<double> &v2) const {
		return (problem.objfun(v1) < problem.objfun(v2));
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
			pagmo_assert(size == s[j].size());
			tmp += s[j][i];
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
				tmp += (s[i][k] - s[j][k]) * (s[i][k] - s[j][k]);
			}
			tmp = std::sqrt(tmp);
			if (tmp > retval) {
				retval = tmp;
			}
		}
	}
	return retval;
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
	simplex s(simplex_size,std::vector<double>(prob_size));
	for (size_t i = 0; i < simplex_size; ++i) {
		s[i] = retval[i].getDecisionVector();
	}
	// Perform the NM method for a number of times equal to m_gen.
	for (size_t gen = 0; gen < m_gen; ++gen) {
		const double diameter = simplex_diameter(s);
		// NOTE: hard coded value, is it worth to make it configurable?
		if (diameter < 1E-6) {
			//std::cout << "small diameter: " << diameter << '\n';
			for (size_t i = 1; i < simplex_size; ++i) {
				s[i] = Individual(problem).getDecisionVector();
			}
		}
//		if (drng() < .01) {
//			std::cout << "shuffle!\n";
//			for (size_t i = 1; i < simplex_size; ++i) {
//				for (size_t j = 0; j < prob_size; ++j) {
//					ffd
//				}
//			}
//		}
		// First order the vertices of the simplex according to fitness.
		std::sort(s.begin(),s.end(),sorter(problem));
		// Compute the center of mass, excluding the worst point.
		const std::vector<double> x0 = center_mass(s);
		// Compute a reflection.
		std::vector<double> xr = sub_mult_add(x0,s.back(),m_alpha,x0);
		check_bounds(xr,problem);
		// Calculate fitness values. fitness_target is relative to the vertex immediately before the worst one.
		const double fitness_r = problem.objfun(xr), fitness_best = problem.objfun(s[0]), fitness_target = problem.objfun(s[simplex_size - 2]);
		if (fitness_r >= fitness_best && fitness_r < fitness_target) {
			s.back() = xr;
		} else if (fitness_r < fitness_best) {
			std::vector<double> xe = sub_mult_add(x0,s.back(),m_gamma,x0);
			check_bounds(xe,problem);
			const double fitness_e = problem.objfun(xe);
			if (fitness_e < fitness_r) {
				s.back() = xe;
			} else {
				s.back() = xr;
			}
		} else {
			std::vector<double> xc = sub_mult_add(x0,s.back(),m_rho,s.back());
			check_bounds(xc,problem);
			const double fitness_c = problem.objfun(xc), fitness_worst = problem.objfun(s.back());
			if (fitness_c <= fitness_worst) {
				s.back() = xc;
			} else {
				for (size_t i = 1; i < simplex_size; ++i) {
					s[i] = sub_mult_add(s[i],s.front(),m_sigma,s.front());
					check_bounds(s[i],problem);
				}
			}
		}
	}
	// In retval overwrite the individual that have evolved.
	for (size_t i = 0; i < simplex_size; ++i) {
		retval[i] = Individual(s[i],problem.objfun(s[i]));
	}
	return retval;
}
