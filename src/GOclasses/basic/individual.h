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

#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <iostream>
#include <vector>

#include "../../config.h"
#include "GOproblem.h"
#include "rng.h"

class __PAGMO_VISIBLE Individual{

public:
	Individual(const std::vector<double> &, const std::vector<double> &, const double &);
	Individual(const GOProblem &);
	double getFitness() const {return fitness;}
	const std::vector<double> &getDecisionVector() const {return x;}
	const std::vector<double> &getVelocity() const {return v;}
	friend std::ostream &operator<<(std::ostream &, const Individual &);
private:
	void init(const GOProblem &);
	std::vector<double> x;  //this is the "chromosome" or "decision vector"
	std::vector<double> v;  //this is the "velocity" or "heading" of each individual
	double fitness;
};

inline std::ostream &operator<<(std::ostream &s, const Individual &ind) {
	for (size_t i = 0; i < ind.x.size(); ++i) {
		s << ind.x[i] << " ";
	}
	return s;
}

#endif
