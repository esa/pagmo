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

// 16/05/08 Created by Dario Izzo.

#ifndef PAGMO_INDIVIDUAL_H
#define PAGMO_INDIVIDUAL_H

#include <iostream>
#include <vector>

#include "../../config.h"
#include "GOproblem.h"
#include "rng.h"

class __PAGMO_VISIBLE Individual {

public:
	///Constructs an Individual with x randomly placed within a GOProblem bounds and a random velocity of maximum magnitude (UB-LB).
	Individual(const GOProblem &);
	///Constructs an Individual with x and v, fitness will be calculated from the problem.
	Individual(const GOProblem &, const std::vector<double> &, const std::vector<double> &);
	///Constructs an Individual with x and an empty velocity of size x.size(). Fitness will be calculated from the problem.
	Individual(const GOProblem &, const std::vector<double> &);
	Individual &operator=(const Individual &);
	///Returns the Individual fitness.
	double getFitness() const {return fitness;}
	///Returns the Individual chromosome (position).
	const std::vector<double> &getDecisionVector() const {return x;}
	///Returns the Individual velocity.
	const std::vector<double> &getVelocity() const {return v;}
private:
	///Operator << for the Individual.
	friend std::ostream &operator<<(std::ostream &, const Individual &);
	///Individual chromosome (position).
	std::vector<double>			x;
	///Individual velocity.
	std::vector<double>			v;
	///Individual fitness.
	double						fitness;
};

inline std::ostream &operator<<(std::ostream &s, const Individual &ind) {
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

#endif
