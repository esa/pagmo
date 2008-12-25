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
	//Legacy constructor - when using it, remember than all members must be properly initialised afterwards
	GOProblem() {};
	//Constructor with vector bounds initialisers
	GOProblem(int, const std::vector<double> &, const std::vector <double> &);
	//Constructor with array bounds initialisers
	GOProblem(int, const double *, const double *);
	//Virtual destructor - required because the class contains a pure virtual member function
	virtual ~GOProblem() {};

	//Bounds getters and setters (via copy)
	void getBounds(std::vector<double> &, std::vector <double> &) const;
	void setBounds(const std::vector<double> &, const std::vector <double> &);
	void setBounds(const double *, const double *);  //this requires setDimension to be called first

	//Bounds getters and setters via reference
	//read only mode
	const std::vector<double> &getLB() const {return LB;}
	const std::vector<double> &getUB() const {return UB;}

	//Dimension getter and setter
	int getDimension() const {return dimension;}
	void setDimension(int D) {dimension = D;}

	//The objective function - must be implemented in subclasses
	virtual double objfun(const std::vector<double> &) = 0;

private:
	int dimension;
	std::vector<double> LB;
	std::vector<double> UB;

};

#endif
