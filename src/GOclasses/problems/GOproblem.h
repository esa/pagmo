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

// 04/06/08 Created by Dario Izzo.

#ifndef PAGMO_GOPROBLEM_H
#define PAGMO_GOPROBLEM_H

#include <iostream>
#include <string>
#include <typeinfo>
#include <vector>

#include "../../../config.h"
#include "../../atomic_counters/atomic_counters.h"
#include "../../exceptions.h"

class __PAGMO_VISIBLE GOProblem {
	friend std::ostream &operator<<(std::ostream &, const GOProblem &);
public:
	// Constructor from size; construct problem of size n, with lower bounds to zero and upper bounds to one
	GOProblem(int n) {
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
	// Virtual destructor - required because the class contains a pure virtual member function
	virtual ~GOProblem() {}
	// Bounds getters and setters via reference
	const std::vector<double> &getLB() const {return LB;}
	const std::vector<double> &getUB() const {return UB;}
	// Dimension getter
	size_t getDimension() const {return LB.size();}
	// The objective function - must be implemented in subclasses
	double objfun(const std::vector<double> &) const;
	virtual GOProblem *clone() const = 0;
	std::string id_name() const {return typeid(*this).name();}
	virtual void pre_evolution() const {}
	virtual void post_evolution() const {}
	void set_lb(int n, const double &value) {
		if (n < 0 || (size_t)n >= LB.size()) {
			pagmo_throw(index_error,"invalid index for lower bound");
		}
		if (UB[n] <= value) {
			pagmo_throw(value_error,"invalid value for lower bound");
		}
		LB[n] = value;
	}
	void set_ub(int n, const double &value) {
		if (n < 0 || (size_t)n >= UB.size()) {
			pagmo_throw(index_error,"invalid index for upper bound");
		}
		if (LB[n] >= value) {
			pagmo_throw(value_error,"invalid value for upper bound");
		}
		UB[n] = value;
	}
	static size_t objfun_calls();
protected:
	virtual double objfun_(const std::vector<double> &) const = 0;
	// Constructor with array bounds initialisers
	GOProblem(const size_t &d, const double *l, const double *u):LB(l, l + d),UB(u, u + d) {
		check_boundaries();
	}
	// Constructor with vectors initialisers
	GOProblem(const std::vector<double>& l, const std::vector<double>& u):LB(l),UB(u) {
		if (l.size() != u.size()) {
			pagmo_throw(value_error,"size mismatch in base problem constructor from vectors");
		}
		check_boundaries();
	}
	// These need to be protected and cannot be private and const as some problems need to deifne the LB and UB at run time (i.e. LJ)
	std::vector<double> LB;
	std::vector<double> UB;
private:
	void check_boundaries() const {
		const size_t size = LB.size();
		for (size_t i = 0; i < size; ++i) {
			if (LB[i] >= UB[i]) {
				pagmo_throw(value_error,"inconsistent boundaries in base problem");
			}
		}
	}
	static PaGMO::atomic_counter_size_t m_objfun_counter;
};

std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const GOProblem &);

#endif
