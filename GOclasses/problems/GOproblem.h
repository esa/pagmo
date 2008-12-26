/*
 *  GOProblem.h
 *  SeGMO
 *
 *  Created by Dario Izzo on 6/4/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef GOPROBLEM_H
#define GOPROBLEM_H

#include <vector>

class GOProblem {
public:
	//Constructor with array bounds initialisers
	GOProblem(const size_t &d, const double *l, const double *u):dimension(d),LB(l,l+d),UB(u,u+d) {};
	//Virtual destructor - required because the class contains a pure virtual member function
	virtual ~GOProblem() {};

	// Bounds getters and setters via reference
	const std::vector<double> &getLB() const {return LB;}
	const std::vector<double> &getUB() const {return UB;}

	// Dimension getter
	size_t getDimension() const {return dimension;}

	// The objective function - must be implemented in subclasses
	virtual double objfun(const std::vector<double> &) = 0;

private:
	const size_t dimension;
	const std::vector<double> LB;
	const std::vector<double> UB;

};

#endif
