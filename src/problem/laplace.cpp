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

#include "laplace.h"
#include "../AstroToolbox/mga_dsm.h"
#include "../AstroToolbox/misc4Tandem.h"
#include "../keplerian_toolbox/epoch.h"

namespace pagmo { namespace problem {

/// Constructor
/**
* Instantiates one of the possible TandEM problems
* \param[in] problemid This is an integer number from 1 to 24 encoding the fly-by sequence to be used (default is EVEES). Check http://www.esa.int/gsp/ACT/inf/op/globopt/TandEM.htm for more information
* \param[in] tof_ (in years) This is a number setting the constraint on the total time of flight (10 from the GTOP database). If -1 (default) an unconstrained problem is instantiated
*/
laplace::laplace(const std::vector<int> &seq):base(-10000.0,10000.0,4*seq.size() - 2), problem(0)
{
	if (seq.size() < 2)
	{
		pagmo_throw(value_error,"fly-by sequence must contain at least two planets!! Earth and Jupiter");
	}
	if (seq[0]!=3)
	{
		pagmo_throw(value_error,"Starting planet must be the Erath (3) for the Laplace mission");
	}
	if (seq.back()!=5)
	{
		pagmo_throw(value_error,"Final planet must be the Jupiter (5) for the Laplace mission");
	}
	//Set the sequence as defined by the user
	problem.reset(new mgadsmproblem(orbit_insertion,&seq[0],seq.size(),0,0,0,.97,4 * 71492.0));
	// Set the bounds.
	set_lb(0,5475.0);
	set_ub(0,9132.0);
	set_lb(1,0.1);
	set_ub(1,3.5);
	set_lb(2,0.0);
	set_ub(2,1.0);
	set_lb(3,0.0);
	set_ub(3,1.0);
	for (size_t i = 0; i < seq.size() - 1; ++i)
	{
		set_lb(i+4,20);
		set_ub(i+4,2500);
	}
	for (size_t i = 0; i < seq.size() - 1; ++i)
	{
		set_lb(i+3+seq.size(),0.01);
		set_ub(i+3+seq.size(),0.99);
	}
	for (size_t i = 0; i < seq.size() - 2; ++i)
	{
		set_lb(i+2*(seq.size()+ 1),1.05);
		set_ub(i+2*(seq.size()+ 1),100);
	}
	for (size_t i = 0; i < seq.size() - 2; ++i)
	{
		set_lb(i+3*seq.size(),-M_PI);
		set_ub(i+3*seq.size(),M_PI);
	}
}

laplace::laplace(const laplace &other):base(other),problem(new mgadsmproblem(*(other.problem))) {}

/// Clone method.
base_ptr laplace::clone() const
{
	return base_ptr(new laplace(*this));
}

/// Implementation of the objective function.
void laplace::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	f[0] = 0;
	MGA_DSM(x, *problem, f[0]);
	const size_t sequence_size = (x.size() + 2) / 4;
	double totaltime = 0;
	for (size_t i = 0; i < sequence_size - 1 ; ++i) {
		totaltime += x[i + 4];
	}
	const double delta = totaltime - 8 * 365.25;
	// Penalise trajectory longer than 8 years by 200 meters/s per month.
	f[0] = std::max<double>(f[0],f[0] + 0.2 / 30 * delta);
}

/// Prints the chromosome in a pretty way!!!
/**
 Outputs to screen the trajectory data in a humar readable format.
 */

 std::string laplace::pretty(const std::vector<double> &x) const
{
	double obj = 0;
	MGA_DSM(x, *problem, obj);
	std::ostringstream s;
	s.precision(15);
	s << std::scientific;
	const size_t seq_size = (x.size() + 2) / 4;
	pagmo_assert((x.size() + 2) % 4 == 0 && seq_size >= 2);
	pagmo_assert(problem->sequence.size() == seq_size);
	s << "Flyby sequence:        ";
	for (size_t i = 0; i < seq_size; ++i) {
		s << problem->sequence[i];
	}
	s << '\n';
	s << "Departure epoch (mjd2000):     " << x[0] << '\n';
	s << "Departure epoch:     " << ::kep_toolbox::epoch(x[0],::kep_toolbox::epoch::MJD2000) << '\n';
	s << "Vinf polar components: ";
	for (size_t i = 0; i < 3; ++i) {
		s << x[i + 1] << ' ';
	}
	s << '\n';
	double totaltime = 0;
	for (size_t i = 0; i < seq_size - 1; ++i) {
		s << "Leg time of flight:    " << x[i + 4] << '\n';
		totaltime += x[i + 4];
	}
	s << "Total time of flight:  " << totaltime << '\n';
	for (size_t i = 0; i < seq_size - 2; ++i) {
		s << "Flyby radius:          " << x[i + 2 * (seq_size + 1)] << '\n';
	}
	for (size_t i = 0; i < seq_size - 2; ++i) {
		s << "Vinf at flyby:         " << std::sqrt(problem->vrelin_vec[i]) << '\n';
	}
	for (size_t i = 0; i < seq_size - 1; ++i) {
		s << "dsm" << i+1 << ":         " << problem->DV[i+1] << '\n';
	}
	return s.str();
}

/// Implementation of the sparsity structure.
/**
 * This is necessary and cannot be left to the automatic algorithm implemented in problem::base
 * as the numerical difficulties introduced by the objective function definition through a logarithm
 * makes automated detection unreliable (e.g. also SNOPT algorithm fails).
 * CLearly, as the problem is box constrained no sarsity is present.
 */

void laplace::set_sparsity(int &lenG, std::vector<int> &iGfun, std::vector<int> &jGvar) const
{
	lenG=(4*problem->sequence.size()-2);
	iGfun.resize(lenG);
	jGvar.resize(lenG);
	for (int i = 0; i<lenG; ++i)
	{
		iGfun[i] = 0;
		jGvar[i] = i;
	}
}

}} //namespaces
