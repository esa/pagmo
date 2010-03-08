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

#include <boost/integer_traits.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/random/uniform_real.hpp>
#include <cstddef>
#include <exception>
#include <iterator>
#include <sstream>
#include <string>

#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../types.h"
#include "base.h"
#include "nelder_mead.h"

namespace pagmo
{
namespace algorithm {

/// Constructor.
/**
 * Allows to specify all the parameters of the algorithm.
 *
 * @param[in] n_it number of iterations.
 * @param[in] thresh threshold for simplex rebuilding.
 * @param[in] alpha reflection coefficient.
 * @param[in] gamma expansion coefficient.
 * @param[in] rho contraction coefficient.
 * @param[in] sigma shrink coefficient.
 *
 * @throws value_error if number of iterations is not positive.
 * @throws value_error if alpha is not positive.
 * @throws value_error if gamma is not greater than one.
 * @throws value_error if rho is not in the ]0,1[ range.
 * @throws value_error if sigma is not in the ]0,1[ range.
 */
nelder_mead::nelder_mead(int n_it, const double &thresh, const double &alpha, const double &gamma, const double &rho, const double &sigma):
		base(),m_gen(boost::numeric_cast<std::size_t>(n_it)),m_thresh(thresh),m_alpha(alpha),m_gamma(gamma),m_rho(rho),m_sigma(sigma)
{
	if (n_it <= 0) {
		pagmo_throw(value_error,"number of generations must be positive");
	}
	if (thresh < 0) {
		pagmo_throw(value_error,"simplex rebuild threshold must be nonnegative");
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

/// Clone method.
/**
 * @return algorithm::base_ptr to a copy of this.
 */
base_ptr nelder_mead::clone() const
{
	return base_ptr(new nelder_mead(*this));
}

// Check vertex coordinates.
// If the input vertex's coordinates are outside the boundaries of the provided problem,
// modify the vertex so that it stays inside the boundaries.
bool nelder_mead::check_bounds(decision_vector &x, const problem::base &p)
{
	const decision_vector::size_type size = x.size();
	const decision_vector &lb = p.get_lb(), &ub = p.get_ub();
	pagmo_assert(size == lb.size() && size == ub.size());
	bool retval = true;
	for (decision_vector::size_type i = 0; i < size; ++i) {
		if (x[i] < lb[i]) {
			retval = false;
			x[i] = lb[i];
		}
		if (x[i] > ub[i]) {
			retval = false;
			x[i] = ub[i];
		}
	}
	return retval;
}

// Center of mass of the simplex in the range [begin,end[, considering only the first prob_cont_size
// components of the decision vectors.
void nelder_mead::center_mass(decision_vector &retval, const population::const_iterator &begin,
	const population::const_iterator &end, const problem::base::size_type &prob_cont_size)
{
	pagmo_assert(retval.size() >= prob_cont_size);
	pagmo_assert(std::distance(begin,end) != 0);
	pagmo_assert(prob_cont_size != 0);
	for (population::const_iterator it = begin; it != end; ++it)
	{
		double tmp = 0;
		for (problem::base::size_type j = 0; j < prob_cont_size; ++j) {
			tmp += it->cur_x[j];
		}
		retval[boost::numeric_cast<decision_vector::size_type>(std::distance(begin,it))] =
			tmp / prob_cont_size;
	}
}

// Compute retval = d + c * (a - b), taking into account only the first prob_cont_size elements of the vectors.
void nelder_mead::sub_mult_add(decision_vector &retval, const decision_vector &a, const decision_vector &b,
	const double &c, const decision_vector &d, const problem::base::size_type &prob_cont_size)
{
	const problem::base::size_type size = a.size();
	pagmo_assert(b.size() == size && d.size() == size && retval.size() == size && size >= prob_cont_size);
	for (problem::base::size_type i = 0; i < prob_cont_size; ++i) {
		retval[i] = d[i] + c * (a[i] - b[i]);
	}
}

// Check the ranges in which the components of the vertices of the simplex vary. If they are too small (i.e., smaller than the threshold)
// rebuild the simplex (apart from the best vertex) with random values.
void nelder_mead::simplex_size_check(const population::const_iterator &begin, const population::const_iterator &end,
	population &pop, const problem::base::size_type &prob_cont_size) const
{
	pagmo_assert(prob_cont_size > 0);
	pagmo_assert(std::distance(begin,end) >= 1);
	// Watch, component by component, the min/max value assumed by the vertices
	// of the simplex.
	for (problem::base::size_type i = 0; i < prob_cont_size; ++i) {
		double min = begin->cur_x[i], max = min;
		for (population::const_iterator it = begin + 1; it != end; ++it) {
			if (it->cur_x[i] < min) {
				min = it->cur_x[i];
			}
			if (it->cur_x[i] > max) {
				max = it->cur_x[i];
			}
		}
		pagmo_assert(max >= min);
		if (max - min > m_thresh * (pop.problem().get_ub()[i] - pop.problem().get_lb()[i])) {
			// At least 1 component is greater than the threshold, don't do anything.
			return;
		}
	}
// std::cout << "rebuilding simplex\n";
//std::cout << "best was: " << pop.get_individual(0).cur_x << '\n';
	// All components vary in ranges smaller than the threshold. Rebuild the simplex with random values, apart
	// from the best vertex.
	decision_vector tmp(pop.problem().get_dimension());
	for (population::size_type i = 1; i < boost::numeric_cast<population::size_type>(std::distance(begin,end)); ++i) {
		// Copy over the whole decision vector in order to retain the integer part.
		tmp = pop.get_individual(i).cur_x;
		// Initialise randomly the continuous part of the decision vector.
		for (decision_vector::size_type j = 0; j < prob_cont_size; ++j) {
			tmp[j] = boost::uniform_real<double>(pop.problem().get_lb()[j],pop.problem().get_ub()[j])(m_drng);
		}
		pop.set_x(i,tmp);
	}
//std::cout << "best is: " << pop.get_individual(0).cur_x << '\n';
}

/// Evolve a population.
/**
 * The individuals' decision vectors are interpreted as simplex vertices, and the Nelder-Mead algorithm is run for
 * the number of iterations specified in the constructor.
 *
 * In case there are not enough individuals in order to build a simplex, a
 * value_error exception will be thrown. In case the problem has no continuous part, a value_error exception will be thrown.
 * If there are more individuals than needed, the best individuals are selected
 * to build the simplex. The remaining individuals will survive unaffected the evolution and will be placed at the
 * end of the output population.
 *
 * @param[in,out] pop input/output population.
 *
 * @throws value_error if the population size is less than or equal to the size of the problem.
 * @throws value_error if the problem size is 0.
 * @throws value_error if the problem has no continuous part.
 */
void nelder_mead::evolve(population &pop) const
{
	// Preliminary checks and useful variables.
	const problem::base &problem = pop.problem();
	const population::size_type pop_size = pop.size();
	// We are going to optimise only the continuous part.
	const problem::base::size_type prob_size = problem.get_dimension(), prob_cont_size = prob_size - problem.get_i_dimension();
	const problem::base::f_size_type f_size = problem.get_f_dimension();
	const problem::base::c_size_type c_size = problem.get_c_dimension();
	if (!prob_cont_size) {
		pagmo_throw(value_error,"the problem has no continuous part");
	}
	if (pop_size <= prob_cont_size) {
		pagmo_throw(value_error,"too few individuals in population for Nelder-Mead method");
	}
	if (prob_size < 1) {
		pagmo_throw(value_error,"the dimension of the problem must be at least 1 for Nelder-Mead method");
	}
	if (prob_cont_size == boost::integer_traits<problem::base::size_type>::const_max) {
		pagmo_throw(std::overflow_error,"overflow in Nelder-Mead's simplex size");
	}
	// This is the simplex size.
	const problem::base::size_type simplex_size = prob_cont_size + 1;
	// This will be the center of mass.
	decision_vector x0(prob_size);
	// This will be the vertex used for reflection.
	decision_vector xr(prob_size);
	fitness_vector fr(f_size);
	constraint_vector cr(c_size);
	// This will be the vertex used for expansion.
	decision_vector xe(prob_size);
	fitness_vector fe(f_size);
	constraint_vector ce(c_size);
	// This will be the vertex used for contraction.
	decision_vector xc(prob_size);
	fitness_vector fc(f_size);
	constraint_vector cc(c_size);
	// This will be the vertex used for shrinking.
	decision_vector xs(prob_size);
	// Initial rank of the population.
	pop.rank_current();
	double rescale;
	// Perform the nelder_mead method for a number of times equal to m_gen.
	for (std::size_t gen = 0; gen < m_gen; ++gen) {
		// Iterators to the begin and end of the simplex.
		const population::const_iterator s_begin = pop.begin(),
			s_end = s_begin + boost::numeric_cast<std::iterator_traits<population::const_iterator>::difference_type>(simplex_size);
		pagmo_assert(s_end - s_begin == boost::numeric_cast<std::iterator_traits<population::const_iterator>::difference_type>(simplex_size));
		// Rank the element constituting the simplex.
		pop.rank_current(s_begin,s_end);
		// Check the simplex size.
		simplex_size_check(s_begin,s_end,pop,prob_cont_size);
		// Compute the center of mass, excluding the worst point.
		center_mass(x0,s_begin,s_end - 1,prob_cont_size);
pagmo_assert(check_bounds(x0,problem));
		// Compute a reflection for the worst point.
		rescale = 1;
		do {
			sub_mult_add(xr,x0,(s_end - 1)->cur_x,m_alpha * rescale,x0,prob_cont_size);
			rescale /= 2;
std::cout << "rescale refl: " << rescale << '\n';
if (rescale == 0) {
	std::cout << "null rescale!!!\n";
	check_bounds(xr,problem);
	break;
}
		} while (!check_bounds(xr,problem));
		// Calculate fitness and constraints of the reflected point.
		problem.objfun(fr,xr);
		problem.compute_constraints(cr,xr);
		if (!problem.compare_fc(fr,cr,s_begin->cur_f,s_begin->cur_c) && problem.compare_fc(fr,cr,(s_end - 2)->cur_f,(s_end - 2)->cur_c)) {
			// Reflection generated a point which is not the best but which is better than the second worst: replace the worst.
//std::cout << "Good reflection\n";
			pop.set_x(boost::numeric_cast<population::size_type>(s_end - s_begin - 1),xr);
		} else if (problem.compare_fc(fr,cr,s_begin->cur_f,s_begin->cur_c)) {
			// Reflection generated a point better than the current best. Try to expand.
			rescale = 1;
			do {
				sub_mult_add(xe,x0,(s_end - 1)->cur_x,m_gamma * rescale,x0,prob_cont_size);
				rescale /= 2;
std::cout << "rescale exp: " << rescale << '\n';
if (rescale == 0) {
	std::cout << "null rescale!!!\n";
	check_bounds(xe,problem);
	break;
}
			} while (!check_bounds(xe,problem));
			// Calculate fitness/constraints.
			problem.objfun(fe,xe);
			problem.compute_constraints(ce,xe);
			if (problem.compare_fc(fe,ce,fr,cr)) {
				// Expansion improved again, use the expanded vertex to replace the worst.
				// TODO: use std::distance in these cases?
				pop.set_x(boost::numeric_cast<population::size_type>(s_end - s_begin - 1),xe);
//std::cout << "Best reflection: expansion\n";
			} else {
				// Expansion did not improve, use relfected vertex.
//std::cout << "Best reflection\n";
				pop.set_x(boost::numeric_cast<population::size_type>(s_end - s_begin - 1),xr);
			}
		} else {
			// Reflection did not improve the situation, do the contraction.
			sub_mult_add(xc,x0,(s_end - 1)->cur_x,m_rho,(s_end - 1)->cur_x,prob_cont_size);
			check_bounds(xc,problem);
			// Compute fitness/constraints.
			problem.objfun(fc,xc);
			problem.compute_constraints(cc,xc);
			if (!problem.compare_fc((s_end - 1)->cur_f,(s_end - 1)->cur_c,fc,cc)) {
				// The worst point is no better than the contracted point, replace it.
//std::cout << "Contraction\n";
				pop.set_x(boost::numeric_cast<population::size_type>(s_end - s_begin - 1),xc);
			} else {
//std::cout << "Shrinking\n";
				// The worst point is better than the contracted one. Gotta do the shrinking.
				for (problem::base::size_type i = 1; i < simplex_size; ++i) {
					sub_mult_add(xs,(s_begin + i)->cur_x,s_begin->cur_x,m_sigma,s_begin->cur_x,prob_cont_size);
					check_bounds(xs,problem);
					pop.set_x(boost::numeric_cast<population::size_type>(i),xs);
				}
			}
		}
	}
}

/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string nelder_mead::human_readable_extra() const
{
	std::ostringstream s;
	s << "\tNiter:\t" << m_gen << '\n';
	s << "\talpha:\t" << m_alpha << '\n';
	s << "\tgamma:\t" << m_gamma << '\n';
	s << "\trho:\t" << m_rho << '\n';
	s << "\tsigma:\t" << m_sigma << '\n';
	return s.str();
}

}
}
