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
	GOProblem(int dimension, const std::vector<double>& lower, const std::vector <double>& upper);
	//Constructor with array bounds initialisers
	GOProblem(int dimension, const double lower[], const double upper[]);
	//Virtual destructor - required because the class contains a pure virtual member function
	virtual ~GOProblem() { };

	//Bounds getters and setters (via copy)
	void getBounds(std::vector<double>& lower, std::vector <double>& upper) const;
	void setBounds(const std::vector<double>& lower, const std::vector <double>& upper);
	void setBounds(const double lower[], const double upper[]);  //this requires setDimension to be called first

	//Bounds geeters and setters via reference
	//read only mode
	const std::vector<double>& getLB() const { return LB; }
	const std::vector<double>& getUB() const { return UB; }
	//rewrite mode
	std::vector<double>& getLB() { return LB; }
	std::vector<double>& getUB() { return UB; }

	//Dimension getter and setter
    int getDimension() const { return dimension; };
	void setDimension(const int D) { dimension = D; };

	//The objective function - must be implemented in subclasses
	virtual double objfun(const std::vector<double>&) = 0;

private:
	int dimension;
	std::vector <double> LB;
	std::vector <double> UB;

};	//end class GOProblems

#endif
