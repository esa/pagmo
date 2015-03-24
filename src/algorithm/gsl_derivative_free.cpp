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

#include <algorithm>
#include <boost/numeric/conversion/cast.hpp>
#include <cstddef>
#include <exception>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>

#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../types.h"
#include "base_gsl.h"
#include "gsl_derivative_free.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify all the parameters needed to initialise a GSL minimiser without derivatives, as specified in the GSL documentation.
 * Will fail if max_iter is negative or if at least one of the other parameters is negative.
 *
 * @param[in] max_iter maximum number of iterations the algorithm will be allowed to perform.
 * @param[in] tol desired tolerance in the localisation of the minimum.
 * @param[in] step_size initial step size for the simplex.
 */
gsl_derivative_free::gsl_derivative_free(int max_iter, const double &tol, const double &step_size):
	base_gsl(),m_max_iter(boost::numeric_cast<std::size_t>(max_iter)),m_tol(tol),m_step_size(step_size)
{
	if (step_size <= 0) {
		pagmo_throw(value_error,"step size must be positive");
	}
	if (tol <= 0) {
		pagmo_throw(value_error,"tolerance must be positive");
	}
}

/// Extra information in human-readable format.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string gsl_derivative_free::human_readable_extra() const
{
	std::ostringstream oss;
	oss << "max_iter:" << m_max_iter << ' ';
	oss << "tol:" << m_tol << ' ';
	oss << "step_size:" << m_step_size << ' ';
	return oss.str();
}

/// Evolve method.
/**
 * The best member of the population will be used as starting point for the minimisation process. The algorithm will stop
 * if the size of the simplex falls below the tol parameter, if the maximum number of iterations max_iter is exceeded or if
 * the inner GSL routine call reports an error (which will be logged on std::cout). After the end of the minimisation process,
 * the minimised decision vector will replace the best individual in the population, after being modified to fall within
 * the problem bounds if necessary.
 *
 * @param[in,out] pop population to evolve.
 */
void gsl_derivative_free::evolve(population &pop) const
{
	// Do nothing if the population is empty.
	if (!pop.size()) {
		return;
	}
	// Useful variables.
	const problem::base &problem = pop.problem();
	if (problem.get_f_dimension() != 1) {
		pagmo_throw(value_error,"this algorithm does not support multi-objective optimisation");
	}
	if (problem.get_c_dimension()) {
		pagmo_throw(value_error,"this algorithm does not support constrained optimisation");
	}
	const problem::base::size_type cont_size = problem.get_dimension() - problem.get_i_dimension();
	if (!cont_size) {
		pagmo_throw(value_error,"the problem has no continuous part");
	}
	// Extract the best individual.
	const population::size_type best_ind_idx = pop.get_best_idx();
	const population::individual_type &best_ind = pop.get_individual(best_ind_idx);
	// GSL wrapper parameters structure.
	objfun_wrapper_params params;
	params.p = &problem;
	// Integer part of the temporay decision vector must be filled with the integer part of the best individual,
	// which will not be optimised.
	params.x.resize(problem.get_dimension());
	std::copy(best_ind.cur_x.begin() + cont_size, best_ind.cur_x.end(), params.x.begin() + cont_size);
	params.f.resize(1);
	// GSL function structure.
	gsl_multimin_function gsl_func;
	// Number of function components.
	gsl_func.n = boost::numeric_cast<std::size_t>(cont_size);
	gsl_func.f = &objfun_wrapper;
	gsl_func.params = (void *)&params;
	// Mimimizer.
	gsl_multimin_fminimizer *s = 0;
	// Starting point and step sizes.
	gsl_vector *x = 0, *ss = 0;
	// Here we start the allocations.
	// Recast as size_t here, in order to avoid potential overflows later.
	const std::size_t s_cont_size = boost::numeric_cast<std::size_t>(cont_size);
	// Allocate and check the allocation results.
	x = gsl_vector_alloc(s_cont_size);
	ss = gsl_vector_alloc(s_cont_size);
	const gsl_multimin_fminimizer_type *minimiser = get_gsl_minimiser_ptr();
	pagmo_assert(minimiser);
	s = gsl_multimin_fminimizer_alloc(minimiser,s_cont_size);
	// Check the allocations.
	check_allocs(x,ss,s);
	// Starting point comes from the best individual.
	for (std::size_t i = 0; i < s_cont_size; ++i) {
		gsl_vector_set(x,i,best_ind.cur_x[i]);
	}
	// Set initial step sizes.
	gsl_vector_set_all(ss,m_step_size);
	// Init the solver.
	gsl_multimin_fminimizer_set(s,&gsl_func,x,ss);
	// Iterate.
	std::size_t iter = 0;
	int status;
	double size;
	try {
		do
		{
			status = gsl_multimin_fminimizer_iterate(s);
			++iter;
			if (status) {
				break;
			}
			size = gsl_multimin_fminimizer_size(s);
			status = gsl_multimin_test_size(size, m_tol);
			if (m_screen_output) {
				if (!((iter-1)%20)) {
					std::cout << std::endl << std::left << std::setw(20) << 
					"Iter." << std::setw(20) << 
					"Best " << std::setw(20) <<
					"Size "<< std::endl; 
				}
			std::cout << std::left << std::setprecision(14) << std::setw(20) << 
			 iter << std::setw(20) << 
			 gsl_multimin_fminimizer_minimum(s) << std::setw(20) << 
			 size << std::endl;
			}
		} while (status == GSL_CONTINUE && iter < m_max_iter);
	} catch (const std::exception &e) {
		// Cleanup and re-throw.
		cleanup(x,ss,s);
		throw e;
	} catch (...) {
		// Cleanup and throw.
		cleanup(x,ss,s);
		pagmo_throw(std::runtime_error,"unknown exception caught in gsl_derivative_free::evolve");
	}
	// Free up resources.
	cleanup(x,ss,s);
	// Check the generated individual and change it to respect the bounds as necessary.
	for (problem::base::size_type i = 0; i < cont_size; ++i) {
		if (params.x[i] < problem.get_lb()[i]) {
			params.x[i] = problem.get_lb()[i];
		}
		if (params.x[i] > problem.get_ub()[i]) {
			params.x[i] = problem.get_ub()[i];
		}
	}
	// Replace the best individual.
	pop.set_x(best_ind_idx,params.x);
}

// Check GSL allocations done during evolve().
void gsl_derivative_free::check_allocs(gsl_vector *x, gsl_vector *ss, gsl_multimin_fminimizer *s)
{
	// If at least one allocation failed, cleanup and throw a memory error.
	if (!x || !ss || !s) {
		cleanup(x,ss,s);
		throw std::bad_alloc();
	}
}

// Cleanup allocations done during evolve.
void gsl_derivative_free::cleanup(gsl_vector *x, gsl_vector *ss, gsl_multimin_fminimizer *s)
{
	if (x) {
		gsl_vector_free(x);
	}
	if (ss) {
		gsl_vector_free(ss);
	}
	if (s) {
		gsl_multimin_fminimizer_free(s);
	}
}

}}
