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

#ifndef PAGMO_PROBLEM_MESSENGER_FULL_H
#define PAGMO_PROBLEM_MESSENGER_FULL_H

#include "../config.h"
#include "../types.h"
#include "base.h"
#include "../AstroToolbox/mga_dsm.h"
#include "boost/array.hpp"

namespace pagmo{ namespace problem {

/// Messenger problem (dufficult version, e.g. messenger_full)
/**
 * This is a rather complex interplanetary trajectory problem that is related to the mission
 * Messenger (NASA). Fly-by sequence, spacecraft mass, objective and other parameter are set as to
 * resemble the Messenger spacecraft mission profile. The problem is transcribed as an MGA-DSM
 * problem allowing one chemical manouvre per trajectory leg.
 * Fly-by sequence is Earth-Venus-Venus-Mercury-Mercury-Mercury-Mercury.
 * The problem is also part of the Global Trajectory Optimization database (GTOP)
 *
 * messenger_full is a box constrained single objective, continuous optimization problem of dimension 26.
 *
 * @see http://www.esa.int/gsp/ACT/inf/op/globopt/MessengerFull.html
 * @see http://www.nasa.gov/mission_pages/messenger/timeline/index.html
 * @author Dario Izzo (dario.izzo@esa.int)
 */
class __PAGMO_VISIBLE messenger_full: public base
{
	public:
		messenger_full();
		base_ptr clone() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
	private:
		mgadsmproblem problem;

};

}}

#endif // PAGMO_PROBLEM_MESSENGER_FULL_H
