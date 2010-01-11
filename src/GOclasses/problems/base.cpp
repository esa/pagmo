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

// 09/06/09 Created by Francesco Biscani.

#include <iostream>
#include <string>
#include <vector>

#include "../../atomic_counters/atomic_counters.h"
#include "../../exceptions.h"
#include "base.h"

namespace pagmo
{
namespace problem {

atomic_counter_size_t base::m_objfun_counter(0);

base::base(int n)
{
	if (n <= 0) {
		pagmo_throw(value_error,"size of problem must be positive");
	}
	const size_t size = (size_t)(n);
	LB.resize(size);
	UB.resize(size);
	for (size_t i = 0; i < size; ++i) {
		LB[i] = 0;
		UB[i] = 1;
	}
}

base::base(const size_t &d, const double *l, const double *u):LB(l, l + d),UB(u, u + d)
{
	check_boundaries();
}

// Constructor with vectors initialisers
base::base(const std::vector<double>& l, const std::vector<double>& u):LB(l),UB(u)
{
	if (l.size() != u.size()) {
		pagmo_throw(value_error,"size mismatch in base problem constructor from vectors");
	}
	check_boundaries();
}

void base::check_boundaries() const
{
	const size_t size = LB.size();
	for (size_t i = 0; i < size; ++i) {
		if (LB[i] >= UB[i]) {
			pagmo_throw(value_error,"inconsistent boundaries in base problem");
		}
	}
}

const std::vector<double> &base::get_lb() const
{
	return LB;
}

const std::vector<double> &base::get_ub() const
{
	return UB;
}

void base::set_lb(const std::vector<double> &lb)
{
	const size_t size = lb.size();
	if (size != getDimension()) {
		pagmo_throw(value_error,"size of lower bounds vector is incompatible with problem size");
	}
	for (size_t i = 0; i < size; ++i) {
		if (lb[i] >= UB[i]) {
			pagmo_throw(value_error,"lower bound is not less than upper bound");
		}
	}
	LB = lb;
}

void base::set_ub(const std::vector<double> &ub)
{
	const size_t size = ub.size();
	if (size != getDimension()) {
		pagmo_throw(value_error,"size of upper bounds vector is incompatible with problem size");
	}
	for (size_t i = 0; i < size; ++i) {
		if (ub[i] <= LB[i]) {
			pagmo_throw(value_error,"upper bound is not greater than lower bound");
		}
	}
	UB = ub;
}

size_t base::getDimension() const
{
	return LB.size();
}

std::string base::id_name() const
{
	return typeid(*this).name();
}

void base::set_lb(int n, const double &value)
{
	if (n < 0 || (size_t)n >= LB.size()) {
		pagmo_throw(index_error,"invalid index for lower bound");
	}
	if (UB[n] <= value) {
		pagmo_throw(value_error,"invalid value for lower bound");
	}
	LB[n] = value;
}

void base::set_ub(int n, const double &value)
{
	if (n < 0 || (size_t)n >= UB.size()) {
		pagmo_throw(index_error,"invalid index for upper bound");
	}
	if (LB[n] >= value) {
		pagmo_throw(value_error,"invalid value for upper bound");
	}
	UB[n] = value;
}

double base::objfun(const std::vector<double> &v) const
{
	if (v.size() != getDimension()) {
		pagmo_throw(value_error,"mismatch between the size of the decision vector and the size of the problem when calling the objective function");
	}
	const double retval = objfun_(v);
	// Actually do the increment only if we have fast incrementing capabilities in m_objfun_counter.
	if (m_objfun_counter.is_increment_fast) {
		++m_objfun_counter;
	}
	return retval;
}

bool base::operator==(const base &p) const
{
	const size_t size = getDimension();
	if (typeid(*this) != typeid(p) || size != p.getDimension()) {
		return false;
	}
	for (size_t i = 0; i < size; ++i) {
		if (LB[i] != p.LB[i] || UB[i] != p.UB[i]) {
			return false;
		}
	}
	return true;
}

bool base::operator!=(const base &p) const
{
	return !(*this == p);
}

std::ostream &base::print(std::ostream &s) const
{
	s << "Problem type: " << id_name() << '\n';
	s << "Dimension: " << getDimension() << '\n';
	const size_t size = getDimension();
	s << "Lower bounds:\n[";
	for (size_t i = 0; i < size; ++i) {
		s << LB[i];
		if (i < size - 1) {
			s << ',';
		}
	}
	s << "]\n";
	s << "Upper bounds:\n[";
	for (size_t i = 0; i < size; ++i) {
		s << UB[i];
		if (i < size - 1) {
			s << ' ';
		}
	}
	s << "]\n";
	return s;
}

std::ostream &operator<<(std::ostream &s, const base &p)
{
	return p.print(s);
}

size_t objfun_calls()
{
	if (!base::m_objfun_counter.is_increment_fast) {
		pagmo_throw(not_implemented_error,"fast atomic counters are not available in this version of PaGMO");
	}
	return (base::m_objfun_counter).get_value();
}

void reset_objfun_calls()
{
	base::m_objfun_counter = atomic_counter_size_t();
}

}
}
