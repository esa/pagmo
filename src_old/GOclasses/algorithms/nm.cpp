/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
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

// 03/02/2009: Initial version by Francesco Biscani.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include "../../exceptions.h"
#include "../basic/individual.h"
#include "../basic/population.h"
#include "../problems/base.h"
#include "base.h"
#include "nm.h"

namespace pagmo
{
namespace algorithm {

/// Simple constructor.
/**
 * Initialises the algorithm using the following commonly-used values:
 *  - \f$ \alpha = 1 \f$
 *  - \f$ \gamma = 2 \f$
 *  - \f$ \rho = 0.5 \f$
 *  - \f$ \sigma = 0.5 \f$
 * @param[in] n_gen number of generations.
 * @throws value_error if number of generations is not positive.
 */
nm::nm(int n_gen):
		base(),m_gen(n_gen),m_alpha(1),m_gamma(2),m_rho(.5),m_sigma(.5)
{
	if (n_gen <= 0) {
		pagmo_throw(value_error,"number of generations must be positive");
	}
}

/// Full constructor.
/**
 * Allows to specify all the parameters of the algorithm.
 * @param[in] n_gen number of generations.
 * @param[in] alpha reflection coefficient.
 * @param[in] gamma expansion coefficient.
 * @param[in] rho contraction coefficient.
 * @param[in] sigma shrink coefficient.
 * @throws value_error if number of generations is not positive.
 * @throws value_error if alpha is not positive.
 * @throws value_error if gamma is not greater than one.
 * @throws value_error if rho is not in the ]0,1[ range.
 * @throws value_error if sigma is not in the ]0,1[ range.
 */
nm::nm(int n_gen, const double &alpha, const double &gamma, const double &rho, const double &sigma):
		base(),m_gen(n_gen),m_alpha(alpha),m_gamma(gamma),m_rho(rho),m_sigma(sigma)
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

// Functor to sort the components of a simplex.
struct nm::sorter {
	sorter(const problem::base &p):problem(p) {}
	bool operator()(const std::pair<vertex,double> &a, const std::pair<vertex,double> &b) const {
		return (a.second < b.second);
	}
	const problem::base &problem;
};

/// Check vertex coordinates.
/**
 * If the input vertex's coordinates are outside the boundaries of the provided problem,
 * modify the vertex so that it stays inside the boundaries.
 * @param[in] v vertex to be checked.
 * @param[in] p problem that provides upper/lower boundaries.
 */
void nm::check_bounds(vertex &v, const problem::base &p) const
{
	const size_t size = v.size();
	const std::vector<double> &LB = p.get_lb(), &UB = p.get_ub();
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

/// Center of mass of a simplex.
/**
 * @param[in] s input simplex.
 * @param[out] retval vertex representing the center of mass.
 */
nm::vertex nm::center_mass(const simplex &s) const
{
	const size_t size = s.size() - 1;
	vertex retval(size);
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

/// Compute \f$ d + c \cdot \left( a - b \right) \f$.
/**
 * Where \f$ a \f$, \f$ b \f$ and \f$ d \f$ are vertices, and \f$ c \f$ a real number.
 * @throws assertion_error if a, b and d are not of the same size.
 */
nm::vertex nm::sub_mult_add(const vertex &a, const vertex &b,
                            const double &c, const vertex &d) const
{
	const size_t size = a.size();
	pagmo_assert(b.size() == size && d.size() == size);
	vertex retval(size);
	for (size_t i = 0; i < size; ++i) {
		retval[i] = d[i] + c * (a[i] - b[i]);
	}
	return retval;
}

/// Compute the diameter of a simplex.
/**
 * Internally it will compute all the vertex distances and return the highest one.
 * @param[in] s input simplex.
 * @param[out] retval double representing the diameter.
 */
double nm::simplex_diameter(const simplex &s) const
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

/// Shuffle the coordinates of the simplex, excluding the first vertex.
/**
 * @param s input/output simplex.
 * @param[in] problem problem used to compute the fitnesses of the new vertices.
 */
void nm::shuffle_simplex(simplex &s, const problem::base &problem) const
{
	const size_t s_size = s.size();
	pagmo_assert(s_size > 1);
	std::vector<vertex> simplex_coordinates(s_size - 1, vertex(s_size - 1));
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

/// Evolve a population.
/**
 * The individuals' decision vectors are interpreted as simplex vertices, and the Nelder-Mead algorithm is run for
 * a number of iterations equal to nm::m_gen.
 *
 * In case there are not enough individuals in order to build a simplex, a
 * value_error exception will be thrown. If there are more individuals than needed, the best individuals are selected
 * to build the simplex. The remaining individuals will survive unaffected the evolution and will be placed at the
 * end of the output population.
 * @param[in] pop input population.
 * @param[out] retval output population.
 * @throws value_error if the population size is less than or equal to the size of the problem.
 * @throws value_error if the problem size is 0.
 */
population nm::evolve(const population &pop) const
{
	// Preliminary checks and useful variables.
	const problem::base &problem = pop.problem();
	const size_t pop_size = pop.size(), prob_size = problem.getDimension();
	if (pop_size <= prob_size) {
		pagmo_throw(value_error,"too few individuals in population for Nelder-Mead method");
	}
	if (prob_size < 1) {
		pagmo_throw(value_error,"the dimension of the problem must be at least 1 for Nelder-Mead method");
	}
	population retval(pop);
	retval.sort();
	// Build a simplex from the sorted input population.
	const size_t simplex_size = prob_size + 1;
	simplex s(simplex_size,std::make_pair(vertex(prob_size),0));
	for (size_t i = 0; i < simplex_size; ++i) {
		s[i].first = retval[i].get_decision_vector();
		s[i].second = retval[i].get_fitness();
	}
	// Perform the NM method for a number of times equal to m_gen.
	for (size_t gen = 0; gen < m_gen; ++gen) {
		const double diameter = simplex_diameter(s);
		// NOTE: hard coded value, is it worth to make it configurable?
		if (diameter < 1E-6) {
			//std::cout << "small diameter: " << diameter << '\n';
			for (size_t i = 1; i < simplex_size; ++i) {
				individual tmp(problem);
				s[i].first = tmp.get_decision_vector();
				s[i].second = tmp.get_fitness();
			}
		}
		// First order the vertices of the simplex according to fitness.
		std::sort(s.begin(),s.end(),sorter(problem));
		// With small probability mix the coordinates between the vertices (excluding the best one).
		if (drng() < .001) {
			shuffle_simplex(s,problem);
		}
		// Compute the center of mass, excluding the worst point.
		const vertex x0(center_mass(s));
		// Compute a reflection.
		vertex xr = sub_mult_add(x0,s.back().first,m_alpha,x0);
		check_bounds(xr,problem);
		// Calculate fitness values. fitness_target is relative to the vertex immediately before the worst one.
		const double fitness_r = problem.objfun(xr), fitness_best = problem.objfun(s[0].first), fitness_target = problem.objfun(s[simplex_size - 2].first);
		if (fitness_r >= fitness_best && fitness_r < fitness_target) {
			s.back().first = xr;
			s.back().second = fitness_r;
		} else if (fitness_r < fitness_best) {
			vertex xe = sub_mult_add(x0,s.back().first,m_gamma,x0);
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
			vertex xc = sub_mult_add(x0,s.back().first,m_rho,s.back().first);
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
		retval.setIndividual(i, individual(s[i].first, retval[i].get_velocity(), s[i].second));
	}
	return retval;
}

/// Print to stream a description of the algorithm.
void nm::log(std::ostream &s) const
{
	s << "NM - generations:" << m_gen << " alpha:" << m_alpha << " gamma:" << m_gamma
	<< " rho:" << m_rho << " sigma:" << m_sigma;
}

}
}
