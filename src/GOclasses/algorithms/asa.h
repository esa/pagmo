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

// 16/05/2008: Initial version by Dario Izzo.

#ifndef PAGMO_ALGORITHM_ASA_H
#define PAGMO_ALGORITHM_ASA_H

#include <iostream>
#include <string>

#include "../../config.h"
#include "../basic/population.h"
#include "base.h"

namespace pagmo
{
namespace algorithm {

/// Simulated Annealing with Adaptive Neighbourhood solver (ASA or SA-AN)
/**
 * The implementation comes from the paper in http://www.csit.fsu.edu/~navon/5420a/corana.pdf by Corana et al.
 * It is, essentially, a simulated annealing where the new points are generated perturbing each component separately
 * and the neighbourhood size of each component is controlled by the number of point accepted or rejected according
 * to the Metropolis criteria
 */
class __PAGMO_VISIBLE asa: public base
{
	public:
		/// Constructor.
		/**
		 * Creates the SA-AN algorithm using a minimal amount of parameters. While it does not give full
		 * control over the algorithm, it assumes reasonable default values for parameters difficult to set
		 * for non expert users.
		 * \param[in] niter Maximum number of function evaluation allowed
		 * \param[in] Ts Starting temperature
		 * \param[in] Tf Final temperature
		 */
		asa(int niter, const double &Ts, const double &Tf);

		/// Constructor.
		/**
		 * Creates the SA-AN algorithm specifying all the algorithm parameters. This constructor should be used only
		 * by users familiar with Corana SA-AN algorithm
		 * \param[in] niter Maximum number of function evaluation allowed
		 * \param[in] Ts Starting temperature
		 * \param[in] Tf Final temperature
		 * \param[in] niterT The number of temperature adjustments will be proportional to this number
		 * \param[in] niterR The number of range adjustments will be proportional to this number
		 * \param[in] startRange Initial range (equal for all components of the decision vector)
		 */
		asa(int niter, const double &Ts, const double &Tf, const int niterT, const int niterR, const double startRange);

		/// Algorithm code.
		/**
		* This method contains the actual code of the SA-AN algorithm. It performs the algorithm starting from
		* the first individual of the population passed as input
		* \param[in] pop Population to be evolved (only the first individual will change)
		*/
		virtual population evolve(const population & pop) const;
		virtual asa *clone() const {
			return new asa(*this);
		}
		virtual std::string id_object() const {
			return id_name();
		}
	private:
		virtual void log(std::ostream &) const;
		size_t niterTot;
		size_t niterTemp;
		size_t niterRange;
		double Ts;
		double Tf;
		double StartStep;
};

}
} // Close namespaces.

#endif
