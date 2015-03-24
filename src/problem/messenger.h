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

#ifndef PAGMO_PROBLEM_MESSENGER_H
#define PAGMO_PROBLEM_MESSENGER_H

#include <string>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "../AstroToolbox/mga_dsm.h"
#include "base.h"

namespace pagmo{ namespace problem {

/// Messenger MGA-DSM Problem (reduced version)
/**
 *
 * An interplanetary trajectory transfer transcribed as an MGA-DSM problem allowing one chemical manouvre per trajectory leg.
 * The objective function is defined as total \f$\Delta V\f$ (km/sec) for a Mercury randevouz and
 * the considered fly-by sequence is Earth-Earth-Venus-Venus-Mercury and corresponds to the first
 * part of the real Messenger fly-by sequence.
 *
 * messenger is a box constrained single objective, continuous optimization problem of dimension 18.
 * The problem is also part of the Global Trajectory Optimization database (GTOP)

 * @see http://www.esa.int/gsp/ACT/inf/op/globopt/MessengerFull.html
 * @see http://www.nasa.gov/mission_pages/messenger/timeline/index.html
 * @author Dario Izzo (dario.izzo@esa.int)
 */
class __PAGMO_VISIBLE messenger: public base
{
	public:
		messenger();
		base_ptr clone() const;
		std::string get_name() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		void set_sparsity(int &, std::vector<int> &, std::vector<int> &) const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & problem;
		} 
		static const int sequence[5];
		mgadsmproblem problem;

};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::problem::messenger)

#endif // PAGMO_PROBLEM_MESSENGER_H
