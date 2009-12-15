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

// 05/2009: initial version by Dario Izzo and Francesco Biscani.

#ifndef PAGMO_PROBLEM_EARTH_MARS_LT_H
#define PAGMO_PROBLEM_EARTH_MARS_LT_H

#include <iostream>
#include <vector>

#include "base.h"

namespace pagmo
{
namespace problem {

/// The Earth-Mars Low-thrust problem,Sims-Flanagan transcription
/**
 * In the paper form 1999 by Sims-Flanagan "Preliminary Design of Low-Thrust Interplanetary Missions"
 * a new direct transcription method was presented. The method, alternative to Hermite-Simpson type of
 * transcription, was noted by the authors to have really good convergence properties and speed.
 * The same method can be transported into a global optimisation framework with only minor modifications
 * giving birth to what we call here an impulsive transcription method.
 * earth_mars_lt creates a global optimisation problem that is an impulsive transcription of the OCP
 * describing a low-thrust trajectory from an Earth launch to a Mars randezvous
 */
class __PAGMO_VISIBLE earth_mars_lt: public base
{
	public:
		/// Constructor.
		/**
		* This instantiate a "earth_mars_lt" problem. This is an impulsive transcription of the low-thrust
		* trajectory optimisation from Earth Launch to Mars Randezvous.
		*
		* \param[in] segments number of segments the trajectory is divided into
		* \param[in] mass     spacecraft launch mass in kg
		* \param[in] thrust   maximumm thrust achievable in N
		* \param[in] Isp      thrusters specific impulse in seconds
		*/
		earth_mars_lt(int, const double &, const double &, const double &);
		virtual earth_mars_lt *clone() const {
			return new earth_mars_lt(*this);
		}
		virtual std::string id_object() const {
			return "hippo's problem";
		}
		void human_readable(const std::vector<double> &) const;
	private:
		virtual double objfun_(const std::vector<double> &) const;
		double main_objfun(const std::vector<double> &) const;
		void state_mismatch(const std::vector<double> &, double *, double *, double *, double *) const;
		static void ruv2cart(double *, const double *);
		static void earth_eph(const double &, double *, double *);
		static void mars_eph(const double &, double *, double *);
		static void kick(double *, const double *);
		static void punch(double *, const double *, const double &, const double &);
		static void back_punch(double *, const double *, const double &, const double &);
		static void propagate(double *, double *, const double &);
		static void back_propagate(double *, double *, const double &);
	private:
		int 	n;
		double 	M;
		double 	thrust;
		double	Isp;
};

}
}

#endif
