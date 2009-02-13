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

// 09/06/09 Created by Francesco Biscani.

#include <iostream>
#include <vector>

#include "../../atomic_counters/atomic_counters.h"
#include "../../exceptions.h"
#include "GOproblem.h"

PaGMO::atomic_counter_size_t GOProblem::m_objfun_counter(0);

double GOProblem::objfun(const std::vector<double> &v) const
{
	if (v.size() != getDimension()) {
		pagmo_throw(value_error,"mismatch between the size of the decision vector and the size of the problem when calling the objective function");
	}
	const double retval = objfun_(v);
	++m_objfun_counter;
	return retval;
}

size_t GOProblem::objfun_calls()
{
	return (size_t)(m_objfun_counter);
}

std::ostream &operator<<(std::ostream &s, const GOProblem &p) {
	s << "Problem type: " << p.id_name() << '\n';
	s << "Dimension: " << p.getDimension() << '\n';
	const size_t size = p.getDimension();
	s << "Lower bounds:\n";
	for (size_t i = 0; i < size; ++i) {
		s << p.LB[i];
		if (i < size - 1) {
			s << ' ';
		}
	}
	s << '\n';
	s << "Upper bounds:\n";
	for (size_t i = 0; i < size; ++i) {
		s << p.UB[i];
		if (i < size - 1) {
			s << ' ';
		}
	}
	s << '\n';
	return s;
}
