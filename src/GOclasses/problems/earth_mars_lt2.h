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

#ifndef PAGMO_PROBLEM_EARTH_MARS_LT2_H
#define PAGMO_PROBLEM_EARTH_MARS_LT2_H

#include <iostream>
#include <vector>

#include "base.h"

namespace pagmo
{
namespace problem {

/// The Earth-Mars Low-thrust problem, a novel transcription
/**
 * In the paper form 1999 by Sims-Flanagan "Preliminary Design of Low-Thrust Interplanetary Missions"
 * a new direct transcription method was presented. The method, alternative to Hermite-Simpson type of
 * transcription, was noted by the authors to have really good convergence properties and speed.
 * We do not think the authors realised that their method coul be significantly generalised and extended to Optimal
 * Control Problems in general. Their transcription allows for the dynamic to be explicitly accounted for in the objective function
 * evaluation and not to be considered as a constraint to be solved by the NLP algorithm choosen (as it is in the case of Hermite-Simpson
 * transcriptions). In this problem we implement a generalization of the Sims-Flanagan method that does not use impulsive velocity changes,
 * but rather integrates for each segment the equation of motion using a continuous thrust having fixed inertial direction to be optimised.
 * This creates a problem whose solution is actually a feasible trajectory for the considered spacecraft and needs no further feasibility
 * correction. The problem of an Earth Launch, Mars randevouz is chosen as in the case of earth_mars_lt2. The two problems are indeed
 * comparable both in terms of computational speed and objective function value.
 */
class __PAGMO_VISIBLE earth_mars_lt2: public base
{
	public:
		/// Constructor
		/**
		* This instantiate a "earth_mars_lt2" problem. This is a transcription of the low-thrust
		* trajectory optimisation from Earth Launch to Mars Randezvous using a novel transcription method.
		* \param[in] segments number of segments the trajectory is divided into
		* \param[in] mass     spacecraft launch mass in kg
		* \param[in] thrust   maximumm thrust achievable in N
		* \param[in] Isp      thrusters specific impulse in seconds
		*/
		earth_mars_lt2(int segments, const double & mass, const double & thrust, const double & Isp);
		virtual earth_mars_lt2 *clone() const {
			return new earth_mars_lt2(*this);
		}
		virtual std::string id_object() const {
			return "hippo's problem";
		}

		/// Decision Vector log in human-readable format
		/**
		 * This function uses the standard output device to print the decision vector
		 * in a human-readable format
		 *
		 * \param[in] x decision vector (or chromosome)
		 */
		void human_readable(const std::vector<double> & x) const;
		std::vector<std::vector<double> > visualize(const std::vector<double> &x, const int N) const;

	private:
		virtual double objfun_(const std::vector<double> &) const;
		double main_objfun(const std::vector<double> &) const;
		void state_mismatch(const std::vector<double> &, double *, double *, double *, double *) const;

		///To follow is a number of functions that are declared as private static. They should be somewhere else (i.e. a toolbox)
		static void ruv2cart(double *, const double *);
		static void earth_eph(const double &, double *, double *);
		static void mars_eph(const double &, double *, double *);
		static void kick(double *, const double *);
		void propagate(double *r, double *v, const double &t, double* fixed_thrust) const;
		void back_propagate(double *r, double *v, const double &t, double* fixed_thrust) const;

	private:
		int 	n;		//Number of segments
		double 	M;		//Spacecraft mass
		double 	thrust;		//Spacecraft Thrust
		double	Isp;		//Spacecraft Specific Impulse

		//The following are needed for numerical integartor
		const double eps;		//Accuracy of the integrator (both absolute and relative)
};

}
}

#endif

