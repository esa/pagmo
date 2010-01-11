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

// 16/05/08 Created by Dario Izzo.

#include "../../exceptions.h"
#include "../problems/base.h"
#include "individual.h"

namespace pagmo
{

individual::individual(const problem::base &problem)
		:x(problem.get_lb().size()),
		v(problem.get_lb().size())
{
	static_rng_double drng;

	// Store local references.
	const std::vector<double> &LB = problem.get_lb(), &UB = problem.get_ub();
	const size_t size = LB.size();

	// Fill a new random chromosome and velocity vector.
	for (size_t i = 0; i < size; ++i) {
		x[i] = LB[i] + drng() * (UB[i] - LB[i]);
		v[i] = drng() * (UB[i] - LB[i]);
	}

	// Evaluation of fitness.
	fitness = problem.objfun(x);
}

individual::individual(const problem::base &problem, const std::vector<double> &x_, const std::vector<double> &v_)
		:x(x_),
		v(v_)
{
	if (x.size() != v.size()) {
		pagmo_throw(value_error,"while constructing individual, size mismatch between decision vector and velocity vector");
	}
	check(problem);
	fitness = problem.objfun(x);
}

individual::individual(const problem::base &problem, const std::vector<double> &x_)
		:x(x_),
		v(x_.size())
{
	check(problem);
	fitness = problem.objfun(x);
}

individual &individual::operator=(const individual &i)
{
	if (this != &i) {
		if (i.get_decision_vector().size() != x.size()) {
			pagmo_throw(value_error,"individuals are incompatible");
		}
		x = i.x;
		v = i.v;
		fitness = i.fitness;
	}
	return *this;
}

void individual::check(const problem::base &p) const
{
	const size_t size = x.size();
	if (size != p.getDimension()) {
		pagmo_throw(value_error,"mismatch between individual size and problem size");
	}
	for (size_t i = 0; i < size; ++i) {
		if (x[i] > p.get_ub()[i] || x[i] < p.get_lb()[i]) {
			pagmo_throw(value_error,"individual's decision vector is incompatible with the boundaries of the problem");
		}
	}
}

int individual::compare_by_fitness(const individual& ind1, const individual& ind2)
{
	return ind1.fitness < ind2.fitness;
}

std::ostream &operator<<(std::ostream &s, const individual &ind)
{
	s << std::scientific;
	s.precision(15);
	for (size_t i = 0; i < ind.x.size(); ++i) {
		s << ind.x[i];
		if (i != ind.x.size() - 1) {
			s << ',';
		}
	}
	return s;
}

}
