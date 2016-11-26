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
#include <gsl/gsl_deriv.h>
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
#include "gsl_gradient.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify all the parameters needed to initialise a GSL minimiser with derivatives, as specified in the GSL documentation.
 * Will fail if max_iter is negative or if at least one of the other parameters is negative.
 *
 * @param[in] max_iter maximum number of iterations allowed.
 * @param[in] grad_tol tolerance when testing the norm of the gradient as stopping criterion.
 * @param[in] numdiff_step_size step size for the numerical computation of the gradient.
 * @param[in] tol accuracy of the line minimisation.
 * @param[in] step_size size of the first trial step.
 */
gsl_gradient::gsl_gradient(int max_iter, const double &grad_tol, const double &numdiff_step_size, const double &step_size, const double &tol):
	base_gsl(),
	m_max_iter(boost::numeric_cast<std::size_t>(max_iter)),m_grad_tol(grad_tol),m_numdiff_step_size(numdiff_step_size),
	m_step_size(step_size),m_tol(tol)
{
	if (step_size <= 0) {
		pagmo_throw(value_error,"step size must be positive");
	}
	if (tol <= 0) {
		pagmo_throw(value_error,"tolerance must be positive");
	}
	if (numdiff_step_size <= 0) {
		pagmo_throw(value_error,"step size for numerical differentiation must be positive");
	}
	if (grad_tol <= 0) {
		pagmo_throw(value_error,"gradient tolerance must be positive");
	}
}

// Wrapper for the numerical differentiation of the objective function of a problem via GSL.
double gsl_gradient::objfun_numdiff_wrapper(double x, void *params)
{
	objfun_numdiff_wrapper_params *pars = (objfun_numdiff_wrapper_params *)params;
	pars->x[pars->coord] = x;
	pars->prob->objfun(pars->f,pars->x);
	return pars->f[0];
}

// Write into retval the gradient of the continuous part of the objective function of prob calculated in input.
void gsl_gradient::objfun_numdiff_central(gsl_vector *retval, const problem::base &prob, const decision_vector &input, const double &step_size)
{
	if (input.size() != prob.get_dimension()) {
		pagmo_throw(value_error,"invalid input vector dimension in numerical differentiation of the objective function");
	}
	if (prob.get_f_dimension() != 1) {
		pagmo_throw(value_error,"numerical differentiation of the objective function cannot work on multi-objective problems");
	}
	// Size of the continuous part of the problem.
	const problem::base::size_type cont_size = prob.get_dimension() - prob.get_i_dimension();
	// Structure to pass data to the wrapper.
	objfun_numdiff_wrapper_params pars;
	pars.x = input;
	pars.f.resize(1);
	pars.prob = &prob;
	// GSL function.
	gsl_function F;
	F.function = &objfun_numdiff_wrapper;
	F.params = (void *)&pars;
	double result, abserr;
	// Numerical differentiation component by component.
	for (problem::base::size_type i = 0; i < cont_size; ++i) {
		pars.coord = i;
		gsl_deriv_central(&F,input[i],step_size,&result,&abserr);
		gsl_vector_set(retval,i,result);
	}
}

// Objective function's derivative wrapper.
void gsl_gradient::d_objfun_wrapper(const gsl_vector *v, void *params, gsl_vector *df)
{
	objfun_wrapper_params *par = (objfun_wrapper_params *)params;
	// Size of the continuous part of the problem.
	const problem::base::size_type cont_size = par->p->get_dimension() - par->p->get_i_dimension();
	// Fill up the continuous part of temporary storage with the contents of v.
	for (problem::base::size_type i = 0; i < cont_size; ++i) {
		par->x[i] = gsl_vector_get(v,i);
	}
	// Calculate the gradient.
	objfun_numdiff_central(df,*par->p,par->x,par->step_size);
}

// Simmultaneous function/derivative computation wrapper for the objective function.
void gsl_gradient::fd_objfun_wrapper(const gsl_vector *v, void *params, double *retval, gsl_vector *df)
{
	*retval = objfun_wrapper(v,params);
	d_objfun_wrapper(v,params,df);
}

/// Evolve method.
/**
 * The best member of the population will be used as starting point for the minimisation process. The algorithm will stop
 * if the gradient falls below the grad_tol parameter, if the maximum number of iterations max_iter is exceeded or if
 * the inner GSL routine call reports an error (which will be logged on std::cout). After the end of the minimisation process,
 * the minimised decision vector will replace the best individual in the population, after being modified to fall within
 * the problem bounds if necessary.
 *
 * @param[in,out] pop population to evolve.
 */
void gsl_gradient::evolve(population &pop) const
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
	params.step_size = m_numdiff_step_size;
	// GSL function structure.
	gsl_multimin_function_fdf gsl_func;
	gsl_func.n = boost::numeric_cast<std::size_t>(cont_size);
	gsl_func.f = &objfun_wrapper;
	gsl_func.df = &d_objfun_wrapper;
	gsl_func.fdf = &fd_objfun_wrapper;
	gsl_func.params = (void *)&params;
	// Minimiser.
	gsl_multimin_fdfminimizer *s = 0;
	// This will be the starting point.
	gsl_vector *x = 0;
	// Here we start the allocations.
	// Recast as size_t here, in order to avoid potential overflows later.
	const std::size_t s_cont_size = boost::numeric_cast<std::size_t>(cont_size);
	// Allocate and check the allocation results.
	x = gsl_vector_alloc(s_cont_size);
	const gsl_multimin_fdfminimizer_type *minimiser = get_gsl_minimiser_ptr();
	pagmo_assert(minimiser);
	s = gsl_multimin_fdfminimizer_alloc(minimiser,s_cont_size);
	// Check the allocations.
	check_allocs(x,s);
	// Fill in the starting point (from the best individual).
	for (std::size_t i = 0; i < s_cont_size; ++i) {
		gsl_vector_set(x,i,best_ind.cur_x[i]);
	}
	// Init the solver.
	gsl_multimin_fdfminimizer_set(s,&gsl_func,x,m_step_size,m_tol);
	// Iterate.
	std::size_t iter = 0;
	int status;
	try {
		do
		{
			++iter;
			status = gsl_multimin_fdfminimizer_iterate(s);
			if (status) {
				break;
			}
			status = gsl_multimin_test_gradient(s->gradient,m_grad_tol);
		} while (status == GSL_CONTINUE && iter < m_max_iter);
	} catch (const std::exception &e) {
		// Cleanup and re-throw.
		cleanup(x,s);
		throw e;
	} catch (...) {
		// Cleanup and throw.
		cleanup(x,s);
		pagmo_throw(std::runtime_error,"unknown exception caught in gsl_gradient::evolve");
	}
	// Free up resources.
	cleanup(x,s);
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
void gsl_gradient::check_allocs(gsl_vector *x, gsl_multimin_fdfminimizer *s)
{
	// If at least one allocation failed, cleanup and throw a memory error.
	if (!x || !s) {
		cleanup(x,s);
		throw std::bad_alloc();
	}
}

// Cleanup allocations done during evolve.
void gsl_gradient::cleanup(gsl_vector *x, gsl_multimin_fdfminimizer *s)
{
	if (x) {
		gsl_vector_free(x);
	}
	if (s) {
		gsl_multimin_fdfminimizer_free(s);
	}
}

/// Extra information in human-readable format.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string gsl_gradient::human_readable_extra() const
{
	std::ostringstream oss;
	oss << "max_iter:" << m_max_iter << ' ';
	oss << "step_size:" << m_step_size << ' ';
	oss << "tol:" << m_tol << ' ';
	oss << "grad_step_size:" << m_numdiff_step_size << ' ';
	oss << "grad_tol:" << m_grad_tol << ' ';


	return oss.str();
}

}}
