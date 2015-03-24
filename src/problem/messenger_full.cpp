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

#include <string>
#include <keplerian_toolbox/epoch.h>

#include "messenger_full.h"
#include "../AstroToolbox/mga_dsm.h"

namespace pagmo { namespace problem {

const int messenger_full::sequence[7] = {3, 2, 2, 1, 1, 1, 1};

/// Problem Constructor
/**
 * @see problem::base constructors.
 */
messenger_full::messenger_full():base(26),problem(orbit_insertion,sequence,7,0,0,0,0.704,2440 + 200)
{
	// Set bounds.
	const double lb[26] = {1900, 3,    0, 0, 100, 100, 100, 100, 100, 100, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 1.1, 1.1, 1.05, 1.05, 1.05, -M_PI, -M_PI, -M_PI, -M_PI, -M_PI};
	const double ub[26] = {2200, 4.05, 1, 1, 500, 500, 500, 500, 500, 550, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99,   6,   6,    6,    6,    6,  M_PI,  M_PI,  M_PI,  M_PI,  M_PI};
	set_bounds(lb,lb+26,ub,ub+26);

	// Set sequence

}

/// Clone method.
base_ptr messenger_full::clone() const
{
	return base_ptr(new messenger_full(*this));
}

/// Implementation of the objective function.
void messenger_full::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	MGA_DSM(x, problem,f[0]);
}

/// Outputs a stream with the trajectory data
/**
 * While the chromosome contains all necessary information to describe a trajectory, mission analysits
 * often require a different set of data to evaluate its use. This method outputs a stream with
 * information on the trajectory that is otherwise 'hidden' in the chromosome
 *
 * \param[in] x chromosome representing the trajectory in the optimization process
 * \returns an std::string with launch dates, DV magnitues and other information on the trajectory
 */

std::string messenger_full::pretty(const std::vector<double> &x) const
{
	double obj = 0;
	MGA_DSM(x, problem, obj);
	std::ostringstream s;
	s.precision(15);
	s << std::scientific;
	const size_t seq_size = (x.size() + 2) / 4;
	pagmo_assert((x.size() + 2) % 4 == 0 && seq_size >= 2);
	pagmo_assert(problem.sequence.size() == seq_size);
	s << "Flyby sequence:\t\t\t";
	for (size_t i = 0; i < seq_size; ++i) {
		s << problem.sequence[i];
	}
	s << '\n';
	s << "Departure epoch (mjd2000):\t" << x[0] << '\n';
	s << "Departure epoch:\t\t" << ::kep_toolbox::epoch(x[0],::kep_toolbox::epoch::MJD2000) << '\n';
	s << "Vinf polar components:\t\t";
	for (size_t i = 0; i < 3; ++i) {
		s << x[i + 1] << ' ';
	}
	s << '\n';
	double totaltime = 0;
	for (size_t i = 0; i < seq_size - 1; ++i) {
		s << "Leg time of flight:\t\t" << x[i + 4] << '\n';
		totaltime += x[i + 4];
	}
	s << "Total time of flight:\t\t" << totaltime << '\n';
	for (size_t i = 0; i < seq_size - 2; ++i) {
		s << "Flyby radius:\t\t\t" << x[i + 2 * (seq_size + 1)] << '\n';
	}
	totaltime=x[0];
	for (size_t i = 0; i < seq_size - 2; ++i) {
	totaltime += x[i + 4];
		s << "Flyby date:\t\t\t" << ::kep_toolbox::epoch(totaltime,::kep_toolbox::epoch::MJD2000) << '\n';
	}
	for (size_t i = 0; i < seq_size - 2; ++i) {
		s << "Vinf at flyby:\t\t\t" << std::sqrt(problem.vrelin_vec[i]) << '\n';
	}
	for (size_t i = 0; i < seq_size - 1; ++i) {
		s << "dsm" << i+1 << ":\t\t\t\t" << problem.DV[i+1] << '\n';
	}
	s << "Final DV:\t\t\t" << problem.DV.back() << '\n';
	return s.str();
}

/// Implementation of the sparsity structure.
void messenger_full::set_sparsity(int &lenG, std::vector<int> &iGfun, std::vector<int> &jGvar) const
{
	lenG=26;
	iGfun.resize(26);
	jGvar.resize(26);
	for (int i = 0; i<lenG; ++i)
	{
		iGfun[i] = 0;
		jGvar[i] = i;
	}
}

std::string messenger_full::get_name() const
{
	return "Messenger full";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::messenger_full)
