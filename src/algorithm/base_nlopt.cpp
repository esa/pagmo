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
#include <nlopt.hpp>
#include <sstream>
#include <string>
#include <vector>

#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../types.h"
#include "base.h"
#include "base_nlopt.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Will build a wrapper for the NLopt algorithm algo, specifying the maximum number of iterations allowed and the tolerance (e.g., the relative
 * improvement of the objective function under which the algorithm will stop iterating). The user must also specify whether the chosen NLopt algorithm supports
 * or not constrained optimisation.
 *
 * @param[in] algo NLopt algorithm (e.g., nlopt::LN_COBYLA).
 * @param[in] constrained true if the algorithm supports nonlinear constraints, false otherwise.
 * @param[in] only_ineq true if the algorithm only support inequality constraints
 * @param[in] max_iter stop-criteria (number of iterations)
 * @param[in] ftol stop-criteria (absolute on the obj-fun)
 * @param[in] xtol stop-criteria (absolute on the chromosome)
 */
base_nlopt::base_nlopt(nlopt::algorithm algo, bool constrained, bool only_ineq, int max_iter, const double &ftol, const double &xtol):base(),
	m_algo(algo),m_constrained(constrained),m_only_ineq(only_ineq),m_max_iter(boost::numeric_cast<std::size_t>(max_iter)),m_ftol(ftol),m_xtol(xtol)
{
	if ( (ftol <= 0) || (xtol <= 0) ) {
		pagmo_throw(value_error,"tolerances must be positive");
	}
	//Dummy init for m_opt
	nlopt::opt opt(algo,1);
	m_opt = opt;
}

std::string base_nlopt::human_readable_extra() const
{
	std::ostringstream oss;
	oss << "max_iter:" << m_max_iter << ' ';
	oss << "ftol:" << m_ftol << " ";
	oss << "xtol:" << m_xtol;
	return oss.str();
}

// Objective function wrapper.
double base_nlopt::objfun_wrapper(const std::vector<double> &x, std::vector<double> &grad, void* data)
{
	nlopt_wrapper_data *d = (nlopt_wrapper_data *)data;
	pagmo_assert(d->f.size() == 1);

	// Compute the gradient by central diffs if necessary TODO: also ipopt has this code (or similar)
	// It should be moved elswhere in PaGMO. Plus be aware that here a chromsome outside the bounds
	// can be created, thus invaidating its compatibility with the problem (exception will be thrown)

	if (!grad.empty()) {
		std::copy(x.begin(),x.end(),d->dx.begin());
		double central_diff;
		const double h0=1e-8;
		double h;
		double mem;
		for (size_t i =0; i < d->dx.size(); ++i)
		{
			h = h0 * std::max(1.,fabs(d->dx[i]));
			mem = d->dx[i];
			d->dx[i] += h;
			d->prob->objfun(d->f,d->dx);
			central_diff = d->f[0];
			d->dx[i] -= 2*h;
			d->prob->objfun(d->f,d->dx);
			central_diff = (central_diff-d->f[0]) / 2 / h;
			grad[i] = central_diff;
			d->dx[i] = mem;
		}
	}

	// Calculate the objective function.
	d->prob->objfun(d->f,x);
	// Return the fitness.
	return (d->f)[0];
}


// Constraint function wrapper.
double base_nlopt::constraints_wrapper(const std::vector<double> &x, std::vector<double> &grad, void* data)
{
	nlopt_wrapper_data *d = (nlopt_wrapper_data *)data;
	pagmo_assert(d->c.size() == d->prob->get_c_dimension());

	// Compute the gradient by central diffs (if necessary). TODO: also ipopt has this code (or similar)
	// It should be moved elswhere in PaGMO. Plus be aware that here a chromsome outside the bounds
	// can be created, thus invaidating its compatibility with the problem (exception will be thrown)

	if (!grad.empty()) {
		std::copy(x.begin(),x.end(),d->dx.begin());
		double central_diff;
		const double h0=1e-8;
		double h;
		double mem;
		for (size_t i =0; i < d->dx.size(); ++i)
		{
			h = h0 * std::max(1.,fabs(d->dx[i]));
			mem = d->dx[i];
			d->dx[i] += h;
			d->prob->compute_constraints(d->c,d->dx);
			central_diff = d->c[d->c_comp];
			d->dx[i] -= 2*h;
			d->prob->compute_constraints(d->c,d->dx);
			central_diff = (central_diff-d->c[d->c_comp]) / 2 / h;
			grad[i] = central_diff;
			d->dx[i] = mem;
		}
	}

	// Calculate the constraints.
	d->prob->compute_constraints(d->c,x);
	// Return the constraints component.
	return (d->c)[d->c_comp];
}

// Evolve method.
void base_nlopt::evolve(population &pop) const
{
	// Useful variables.
	const problem::base &problem = pop.problem();
	if (problem.get_f_dimension() != 1) {
		pagmo_throw(value_error,"this algorithm does not support multi-objective optimisation");
	}
	const problem::base::c_size_type c_size = problem.get_c_dimension();
	const problem::base::c_size_type ec_size = problem.get_c_dimension() - problem.get_ic_dimension();
	if (c_size && !m_constrained) {
		pagmo_throw(value_error,"this algorithm does not support constraints");
	}
	if (ec_size && m_only_ineq) {
		pagmo_throw(value_error,"this algorithm does not support equality constraints");
	}
	const problem::base::size_type cont_size = problem.get_dimension() - problem.get_i_dimension();
	if (!cont_size) {
		pagmo_throw(value_error,"the problem has no continuous part");
	}
	// Do nothing if the population is empty.
	if (!pop.size()) {
		return;
	}
	// Extract the best individual and set the inital point
	const population::size_type best_ind_idx = pop.get_best_idx();
	const population::individual_type &best_ind = pop.get_individual(best_ind_idx);

	
	// Structure to pass data to the objective function wrapper.
	nlopt_wrapper_data data_objfun;

	data_objfun.prob = &problem;
	data_objfun.x.resize(problem.get_dimension());
	data_objfun.dx.resize(problem.get_dimension());
	data_objfun.f.resize(1);
	
	// Structure to pass data to the constraint function wrapper.
	std::vector<nlopt_wrapper_data> data_constrfun(boost::numeric_cast<std::vector<nlopt_wrapper_data>::size_type>(c_size));
	for (problem::base::c_size_type i = 0; i < c_size; ++i) {
		data_constrfun[i].prob = &problem;
		data_constrfun[i].x.resize(problem.get_dimension());
		data_constrfun[i].dx.resize(problem.get_dimension());
		data_constrfun[i].c.resize(problem.get_c_dimension());
		data_constrfun[i].c_comp = i;
	}

	// Main NLopt call.
	nlopt::opt opt(m_algo, problem.get_dimension());
	m_opt = opt;
	// Sets local optimizer for aug_lag methods, do nothing otherwise
	set_local(problem.get_dimension());
	m_opt.set_lower_bounds(problem.get_lb());
	m_opt.set_upper_bounds(problem.get_ub());
	m_opt.set_min_objective(objfun_wrapper, &data_objfun);
	for (problem::base::c_size_type i =0; i<ec_size; ++i) {
		m_opt.add_equality_constraint(constraints_wrapper, &data_constrfun[i], problem.get_c_tol().at(i));
	}
	for (problem::base::c_size_type i =ec_size; i<c_size; ++i) {
		m_opt.add_inequality_constraint(constraints_wrapper, &data_constrfun[i], problem.get_c_tol().at(i));
	}

	m_opt.set_ftol_abs(m_ftol);
	m_opt.set_xtol_abs(m_xtol);
	m_opt.set_maxeval(m_max_iter);

	//nlopt::result result;
	double dummy;
	decision_vector x0(best_ind.cur_x);
	m_opt.optimize(x0, dummy);
	pop.set_x(best_ind_idx,x0);
}

void base_nlopt::set_local(size_t d) const
{
	(void)d;
}
}}
