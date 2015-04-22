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

// Dummy print function used to suppress all screen output
void no_screen_output (int, const char []) {}
void default_output(int mode, const char *message)
{
	if (mode & WORHP_PRINT_MESSAGE) {
		printf(" %s\n",message);
	}
	if (mode & WORHP_PRINT_WARNING) {
		printf(" %s\n",message);
	}
	if (mode & WORHP_PRINT_ERROR) {
		fprintf(stderr," %s\n",message);
	}
}

/// Constructor
 /**
 * Constructs a WORHP algorithm
 */
worhp::worhp(const int iter, const double feas, const double opt, const bool screen_output)
{
	// We construct the map between parameters and integers used to set and get them
	define_param_map();

	// We set the screen output member from algorithm::base
	set_screen_output(screen_output);

	// We deactivate WORHP keyboard handler as it does introduce funny problems
	setenv("WORHP_DISABLE_KEYBOARD_HANDLER", "1", 0);

	// We deal with screen output (this is buggy in lworhp 1.8.0, hopefully future releases can fix the problem
	// and we will be able to restore the screen output upon request
	if (m_screen_output) {
		SetWorhpPrint(default_output);
	} else {
		SetWorhpPrint(no_screen_output);
	}

	// We read the algorithm parameters from the xml file. If this is not found
	// we set default values and ignore the issue.
	int status;
	ReadParams(&status, const_cast<char*>("param.xml"), &m_params);
	m_params.MatrixCC = false; // Not sure what this does exactly

	// We set some of the parameters exposed in the constructor
	set_param("TolFeas", feas);
	set_param("TolOpti", opt);
	set_param("MaxIter", iter);
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

	// We check the screen output again as the user may have changed it
	if (m_screen_output) {
		SetWorhpPrint(default_output);
	} else {
		SetWorhpPrint(no_screen_output);
	}

	OptVar opt;
	Control control;
	Params params;
	Workspace workspace;
	WorhpPreInit(&opt, &workspace, &params, &control);

	opt.n = prob.get_dimension(); // number of variables
	opt.m = prob.get_c_dimension(); // number of constraints
	auto n_eq = prob.get_c_dimension() - prob.get_ic_dimension(); // number of equality constraints

	// specify nonzeros of derivative matrixes
	workspace.DF.nnz = opt.n; // dense
	workspace.DG.nnz = opt.n * opt.m; // dense
	workspace.HM.nnz = opt.n;

	WorhpInit(&opt, &workspace, &params, &control);
	assert(control.status == FirstCall);
	params = m_params;

	// Specify a derivative free case
	params.UserDF = false;
	params.UserDG = false;
	params.UserHM = false;
	params.UserHMstructure = false;

	// Initialization of variables
	const auto best_idx = pop.get_best_idx();
	pagmo::decision_vector x = pop.get_individual(best_idx).cur_x;

	for (int i = 0; i < opt.n; ++i) {
		opt.X[i] = x[i];
		opt.Lambda[i] = 0;
		opt.XL[i] = prob.get_lb()[i];
		opt.XU[i] = prob.get_ub()[i];
	}

	// Equality constraints
	for (auto i = 0u; i < n_eq; ++i) {
		opt.Mu[i] = 0;
		opt.GL[i] = 0;
		opt.GU[i] = 0;
	}

	// Inequality constraints
	for (auto i = n_eq; i < unsigned(opt.m); ++i) {
		opt.Mu[i] = 0;
		opt.GL[i] = -params.Infty;
		opt.GU[i] = 0;
	}

	while (control.status < TerminateSuccess && control.status > TerminateError) {
		if (GetUserAction(&control, callWorhp)) {
			Worhp(&opt, &workspace, &params, &control);
		}

		if (GetUserAction(&control, iterOutput))
		{
			IterationOutput(&opt, &workspace, &params, &control);
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

	for (int i = 0; i < opt.n; ++i) {
		x[i] = opt.X[i];
	}
	pop.set_x(best_idx, x);

	StatusMsg(&opt, &workspace, &params, &control);
	WorhpFree(&opt, &workspace, &params, &control);
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

/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string worhp::human_readable_extra() const
{
	std::ostringstream s;
	s << "MaxIter:" << get_param("MaxIter") << " TolFeas:"<<get_param("TolFeas")<< " TolOpti:"<< get_param("TolOpti") << std::endl;
	return s.str();
}

}} // namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::worhp)
