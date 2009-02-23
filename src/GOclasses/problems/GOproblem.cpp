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
#include <string>
#include <vector>

#include "../../atomic_counters/atomic_counters.h"
#include "../../exceptions.h"
#include "GOproblem.h"

PaGMO::atomic_counter_size_t GOProblem::m_objfun_counter(0);

GOProblem::GOProblem(int n)
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

GOProblem::GOProblem(const size_t &d, const double *l, const double *u):LB(l, l + d),UB(u, u + d)
{
	check_boundaries();
}

// Constructor with vectors initialisers
GOProblem::GOProblem(const std::vector<double>& l, const std::vector<double>& u):LB(l),UB(u)
{
	if (l.size() != u.size()) {
		pagmo_throw(value_error,"size mismatch in base problem constructor from vectors");
	}
	check_boundaries();
}

void GOProblem::check_boundaries() const
{
	const size_t size = LB.size();
	for (size_t i = 0; i < size; ++i) {
		if (LB[i] >= UB[i]) {
			pagmo_throw(value_error,"inconsistent boundaries in base problem");
		}
	}
}

const std::vector<double> &GOProblem::getLB() const
{
	return LB;
}

const std::vector<double> &GOProblem::getUB() const
{
	return UB;
}

void GOProblem::setLB(const std::vector<double> &lb)
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

void GOProblem::setUB(const std::vector<double> &ub)
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

size_t GOProblem::getDimension() const
{
	return LB.size();
}

std::string GOProblem::id_name() const
{
	return typeid(*this).name();
}

void GOProblem::set_lb(int n, const double &value) {
	if (n < 0 || (size_t)n >= LB.size()) {
		pagmo_throw(index_error,"invalid index for lower bound");
	}
	if (UB[n] <= value) {
		pagmo_throw(value_error,"invalid value for lower bound");
	}
	LB[n] = value;
}

void GOProblem::set_ub(int n, const double &value) {
	if (n < 0 || (size_t)n >= UB.size()) {
		pagmo_throw(index_error,"invalid index for upper bound");
	}
	if (LB[n] >= value) {
		pagmo_throw(value_error,"invalid value for upper bound");
	}
	UB[n] = value;
}

double GOProblem::objfun(const std::vector<double> &v) const
{
	if (v.size() != getDimension()) {
		pagmo_throw(value_error,"mismatch between the size of the decision vector and the size of the problem when calling the objective function");
	}
	const double retval = objfun_(v);
	++m_objfun_counter;
	return retval;
}

bool GOProblem::operator==(const GOProblem &p) const
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

bool GOProblem::operator!=(const GOProblem &p) const
{
	return !(*this == p);
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

size_t objfun_calls()
{
	return (size_t)(GOProblem::m_objfun_counter);
}
