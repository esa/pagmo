/*****************************************************************************
 *   Copyright (C) 2015 The PaGMO development team,                          *
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


#include "worhp.h"

namespace pagmo { namespace algorithm {

/// Constructor
 /**
 * Constructs a WORHP algorithm
 */
worhp::worhp() {
	CHECK_WORHP_VERSION

	int status;
	params.initialised = false;
	ReadParams(&status, const_cast<char*>("param.xml"), &params);
	assert(status == OK);
	params.MatrixCC = false;
}


/// Evolve implementation.
/**
 * Run the WORHP algorithm.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */
void worhp::evolve(pagmo::population& pop) const {
	if (pop.size() == 0) {
		return;
	}

	const auto& prob = pop.problem();

	if (prob.get_i_dimension() != 0) {
		pagmo_throw(value_error,
		            "The problem has an integer part and WORHP is not suitable to solve it.");
	}

	if (prob.get_f_dimension() != 1) {
		pagmo_throw(value_error,
		            "The problem is not single objective and WORHP is not suitable to solve it");
	}

	OptVar opt;
	opt.initialised = false;
	opt.n = prob.get_dimension(); // number of variables
	opt.m = prob.get_c_dimension(); // number of constraints
	auto n_eq = prob.get_c_dimension() - prob.get_ic_dimension(); // number of equality constraints

	Workspace workspace;
	workspace.initialised = false;

	// specify nonzeros of derivative matrixes
	workspace.DF.nnz = opt.n; // dense
	workspace.DG.nnz = opt.n * opt.m; // dense
	workspace.HM.nnz = opt.n;

	Control control;
	control.initialised = false;

	WorhpInit(&opt, &workspace, &params, &control);
	assert(control.status == FirstCall);

	params.UserDF = false;
	params.UserDG = false;
	params.UserHM = false;

	const auto best_idx = pop.get_best_idx();
	pagmo::decision_vector x = pop.get_individual(best_idx).cur_x;

	for (int i = 0; i < opt.n; ++i) {
		opt.X[i] = x[i];
		opt.Lambda[i] = 0;
		opt.XL[i] = prob.get_lb()[i];
		opt.XU[i] = prob.get_ub()[i];
	}

	// Equality constraints
	for (int i = 0; i < n_eq; ++i) {
		opt.Mu[i] = 0;
		opt.GL[i] = 0;
		opt.GU[i] = 0;
	}

	// Inequality constraints
	for (int i = n_eq; i < opt.m; ++i) {
		opt.Mu[i] = 0;
		opt.GL[i] = -params.Infty;
		opt.GU[i] = 0;
	}

	// Define HM as a diagonal LT-CS-matrix, but only if needed by WORHP
	if (workspace.HM.NeedStructure) {
		for(int i = 0; i < workspace.HM.nnz; ++i) 
		{
			workspace.HM.row[i] = i + 1;
			workspace.HM.col[i] = i + 1;
		}
	}

	while (control.status < TerminateSuccess && control.status > TerminateError) {
		if (GetUserAction(&control, callWorhp)) {
			Worhp(&opt, &workspace, &params, &control);
		}

		if (GetUserAction(&control, iterOutput)) {
			//if (m_screen_output) {
				IterationOutput(&opt, &workspace, &params, &control);
			//}
			DoneUserAction(&control, iterOutput);
		}

		if (GetUserAction(&control, evalF)) {
			for (int i = 0; i < opt.n; ++i) {
				x[i] = opt.X[i];
			}
			auto f = prob.objfun(x);
			opt.F = workspace.ScaleObj * f[0];
			DoneUserAction(&control, evalF);
		}

		if (GetUserAction(&control, evalG)) {
			for (int i = 0; i < opt.n; ++i) {
				x[i] = opt.X[i];
			}
			auto g = prob.compute_constraints(x);
			for (int i = 0; i < opt.m; ++i) {
				opt.G[i] = g[i];
			}
			DoneUserAction(&control, evalG);
		}

		if (GetUserAction(&control, fidif)) {
			WorhpFidif(&opt, &workspace, &params, &control);
		}
	}

	StatusMsg(&opt, &workspace, &params, &control);

	for (int i = 0; i < opt.n; ++i) {
		x[i] = opt.X[i];
	}
	pop.set_x(best_idx, x);
}

/// Clone method.
base_ptr worhp::clone() const
{
	return pagmo::algorithm::base_ptr(new worhp(*this));
}

/// Algorithm name
std::string worhp::get_name() const
{
	return "WORHP";
}





}} // namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::worhp)
