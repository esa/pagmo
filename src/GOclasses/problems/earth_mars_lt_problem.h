/*****************************************************************************
 *   Copyright (C) 2008, 2009 Advanced Concepts Team (European Space Agency) *
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

#ifndef EARTH_MARS_LT_H
#define EARTH_MARS_LT_H

#include <iostream>
#include <vector>

#include "GOproblem.h"

class __PAGMO_VISIBLE earth_mars_lt_problem: public GOProblem {
	public:
		earth_mars_lt_problem(int, const double &, const double &, const double &);
		virtual earth_mars_lt_problem *clone() const {return new earth_mars_lt_problem(*this);}
		virtual std::string id_object() const {return "hippo's problem";}
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

#endif
