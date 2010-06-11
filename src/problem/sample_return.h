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

#ifndef PAGMO_SAMPLE_RETURN_H
#define PAGMO_SAMPLE_RETURN_H

#include <string>

#include "../config.h"
#include "../types.h"
#include "base.h"
#include "../keplerian_toolbox/keplerian_toolbox.h"
#include "../AstroToolbox/mga_dsm.h"


namespace pagmo{ namespace problem {

/// Human mission to asteroids.
 /**
 * \image html orion.jpg "Artist visualization of the ORION concpet."
 * \image latex orion.jpg "Artist visualization of the ORION concpet." width=5cm
 * The interplanetary trajectory for a human mission to asteroids is here transcribed
 * into a global optimization problem. Box-constrained single-objective, continuous. The mission profile chosen is that
 * of an Earth-DSM-Asteroid, waiting time on the asteroid, Asteroid-DSM-Earth. The arrival hyperbolic
 * velocity is considered to be smaller than a threshold. The cuncurrent minimization of the launch DV, and
 * all the DV given by the on-board propulsion system allows the flexibility of not having to choose a launcher
 */
class __PAGMO_VISIBLE sample_return: public base
{
	public:
		sample_return(::kep_toolbox::planet asteroid, const double &Tmax = 600);
		base_ptr clone() const;
		std::string pretty(const std::vector<double> &x) const;
		std::string get_name() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		void set_sparsity(int &, std::vector<int> &, std::vector<int> &) const;
	private:
		const ::kep_toolbox::planet	m_target;
		mutable mgadsmproblem		m_leg1;
		mutable mgadsmproblem		m_leg2;
		mutable std::vector<double>	x_leg1;
		mutable std::vector<double>	x_leg2;
		const double			m_Tmax;

};

}}

#endif // PAGMO_SAMPLE_RETURN_H
