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

class __PAGMO_VISIBLE GOProblem {
	friend std::ostream &operator<<(std::ostream &, const GOProblem &);
	friend size_t objfun_calls();
public:
	// Virtual destructor - required because the class contains a pure virtual member function
	virtual ~GOProblem() {}
	// Bounds getters and setters via reference
	const std::vector<double> &getLB() const;
	const std::vector<double> &getUB() const;
	void setLB(const std::vector<double> &);
	void setUB(const std::vector<double> &);
	// Dimension getter
	size_t getDimension() const;
	double objfun(const std::vector<double> &) const;
	virtual GOProblem *clone() const = 0;
	std::string id_name() const;
	virtual void pre_evolution() const {}
	virtual void post_evolution() const {}
	void set_lb(int, const double &);
	void set_ub(int, const double &);
	bool operator==(const GOProblem &) const;
	bool operator!=(const GOProblem &) const;
protected:
	// Print function: called by operator<<, can be re-implemented in sublcasses.
	virtual std::ostream &print(std::ostream &) const;
	// The objective function - must be implemented in subclasses
	virtual double objfun_(const std::vector<double> &) const = 0;
	// Constructor from size; construct problem of size n, with lower bounds to zero and upper bounds to one
	GOProblem(int);
	// Constructor with array bounds initialisers
	GOProblem(const size_t &, const double *, const double *);
	// Constructor with vectors initialisers
	GOProblem(const std::vector<double> &, const std::vector<double> &);
	// These need to be protected and cannot be private and const as some problems need to deifne the LB and UB at run time (i.e. LJ)
	std::vector<double> LB;
	std::vector<double> UB;
private:
	void check_boundaries() const;
	static PaGMO::atomic_counter_size_t m_objfun_counter;
};

std::ostream __PAGMO_VISIBLE_FUNC &operator<<(std::ostream &, const GOProblem &);

size_t __PAGMO_VISIBLE_FUNC objfun_calls();

#endif
