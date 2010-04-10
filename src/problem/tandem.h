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

#ifndef PAGMO_PROBLEM_TANDEM_H
#define PAGMO_PROBLEM_TANDEM_H

#include "../config.h"
#include "../types.h"
#include "base.h"
#include "../AstroToolbox/mga_dsm.h"
#include "../AstroToolbox/misc4Tandem.h"

namespace pagmo{ namespace problem {

/// TandEM problem
/**
 * This interplanetary trajectory problem has 25 different instances, depending on the fly-by sequence adopted.
 * The mission data are taken from a joint ESA/NASA working group that performed, together with industrial partners
 * a preliminary trajectory design during the first months of 2009.
 * Please refer to http://www.esa.int/gsp/ACT/inf/op/globopt/TandEM.htm to select the proper instance. A default
 * choice is 6, which corrsponds to an EVEES fly-by sequence. The possibility of adding one additional constraint
 * on the time-of-flight is given. If such a constraint is specified the components x[4]-x[7] are no longer time of flights
 * but they represent the percentage of the remaining time of flight to be used in one leg.
 *
 * The problem is transcribed as an MGA-DSM problem allowing one chemical manouvre per trajectory leg.
 * The objective function is defined as \f$ f = -\log(m_f)\f$ where \f$m_f\f$ is the final spacecraft mass.
 * The problem is also part of the Global Trajectory Optimization database (GTOP)
 *
 * tandem is a box constrained single objective, continuous optimization problem of dimension 18.
 *
 * @see http://www.esa.int/gsp/ACT/inf/op/globopt/TandEM.htm
 * @author Dario Izzo (dario.izzo@esa.int)
 */
class __PAGMO_VISIBLE tandem: public base
{
	public:
		tandem(const int problemid = 6, const double tof_ = -1);
		base_ptr clone() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		void set_sparsity(int &, std::vector<int> &, std::vector<int> &) const;
	private:
		mgadsmproblem problem;
		const double tof;
		mutable vector<double> copy_of_x;

};

}}

#endif // PAGMO_PROBLEM_TANDEM_H
