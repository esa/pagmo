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

#include <algorithm>
#include <boost/numeric/conversion/cast.hpp>
#include <cstddef>
#include <nlopt.h>
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
 * @param[in] algo NLopt algorithm (e.g., NLOPT_LN_COBYLA).
 * @param[in] constrained true if the algorithm supports nonlinear constraints, false otherwise.
 * @param[in] max_iter maximum number of iterations.
 * @param[in] tol optimality tolerance.
 */
base_nlopt::base_nlopt(nlopt_algorithm algo, bool constrained, int max_iter, const double &tol):base(),
	m_algo(algo),m_constrained(constrained),m_max_iter(boost::numeric_cast<std::size_t>(max_iter)),m_tol(tol),m_last_status(0)
{
	if (tol <= 0) {
		pagmo_throw(value_error,"tolerance must be positive");
	}
}

std::string base_nlopt::human_readable_extra() const
{
	std::ostringstream oss;
	oss << "\tmax_iter:\t\t" << m_max_iter << '\n';
	oss << "\ttol:\t\t\t" << m_tol << '\n';
	oss << "\tlast status:\t\t" << m_last_status << '\n';
	return oss.str();
}

// Objective function wrapper.
double base_nlopt::objfun_wrapper(int n, const double *x, double *, void *data)
{
	nlopt_wrapper_data *d = (nlopt_wrapper_data *)data;
	const problem::base::size_type cont_size = boost::numeric_cast<problem::base::size_type>(n);
	pagmo_assert(cont_size == d->prob->get_dimension() - d->prob->get_i_dimension());
	pagmo_assert(d->f->size() == 1);
	// Copy over the decision vector in the data structure.
	std::copy(x,x + cont_size,d->x->begin());
	// Calculate the objective function.
	d->prob->objfun(*d->f,*d->x);
	// Return the fitness.
	return (*d->f)[0];
}

// Constraint function wrapper.
double base_nlopt::constraints_wrapper(int n, const double *x, double *, void *datum)
{
	nlopt_wrapper_data *d = (nlopt_wrapper_data *)datum;
	const problem::base::size_type cont_size = boost::numeric_cast<problem::base::size_type>(n);
	pagmo_assert(cont_size == d->prob->get_dimension() - d->prob->get_i_dimension());
	pagmo_assert(d->c->size() == d->prob->get_c_dimension());
	// Copy over the decision vector in the data structure.
	std::copy(x,x + cont_size,d->x->begin());
	// Calculate the constraints.
	d->prob->compute_constraints(*d->c,*d->x);
	// Return the constraints component.
	return (*d->c)[d->c_comp];
}

// Evolve method.
void base_nlopt::evolve(population &pop) const
{
	// Useful variables.
	const problem::base &problem = pop.problem();
	if (problem.get_f_dimension() != 1) {
		pagmo_throw(value_error,"this algorithm does not support multi-objective optimisation");
	}
	const problem::base::c_size_type c_size = problem.get_ic_dimension();
	if (problem.get_c_dimension() != c_size) {
		pagmo_throw(value_error,"this algorithm does not support equality constraints");
	}
	if (problem.get_c_dimension() && !m_constrained) {
		pagmo_throw(value_error,"this algorithm does not support constraints");
	}
	const problem::base::size_type cont_size = problem.get_dimension() - problem.get_i_dimension();
	if (!cont_size) {
		pagmo_throw(value_error,"the problem has no continuous part");
	}
	// Do nothing if the population is empty.
	if (!pop.size()) {
		return;
	}
	// Extract the best individual.
	const population::size_type best_ind_idx = pop.get_best_idx();
	const population::individual_type &best_ind = pop.get_individual(best_ind_idx);
	// Structure to pass data to the objective function.
	nlopt_wrapper_data data_objfun;
	decision_vector objfun_x(problem.get_dimension());
	// Copy over the integer part.
	std::copy(best_ind.cur_x.begin() + cont_size, best_ind.cur_x.end(), objfun_x.begin() + cont_size);
	fitness_vector objfun_f(1);
	data_objfun.prob = &problem;
	data_objfun.x = &objfun_x;
	data_objfun.f = &objfun_f;
	// Structure to pass data to the constraint functions.
	std::vector<nlopt_wrapper_data> data_constrfun(boost::numeric_cast<std::vector<nlopt_wrapper_data>::size_type>(c_size));
	decision_vector constrfun_x(problem.get_dimension());
	constraint_vector constrfun_c(c_size);
	for (problem::base::c_size_type i = 0; i < c_size; ++i) {
		data_constrfun[i].prob = &problem;
		data_constrfun[i].x = &constrfun_x;
		data_constrfun[i].c = &constrfun_c;
		data_constrfun[i].c_comp = i;
	}
	// Main NLopt call.
	double retval;
	decision_vector x(best_ind.cur_x);
	m_last_status = nlopt_minimize_constrained(
		m_algo,
		boost::numeric_cast<int>(cont_size),
		objfun_wrapper,
		(void *)&data_objfun,
		boost::numeric_cast<int>(c_size),
		constraints_wrapper,
		(void *)&data_constrfun[0],
		sizeof(nlopt_wrapper_data),
		&problem.get_lb()[0],
		&problem.get_ub()[0],
		&x[0],
		&retval,
		-HUGE_VAL,
		m_tol,0,0,NULL,boost::numeric_cast<int>(m_max_iter),0
	);
	pop.set_x(best_ind_idx,x);
}

}}
