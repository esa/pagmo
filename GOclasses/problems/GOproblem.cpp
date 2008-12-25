/*
 *  GOProblem.cpp
 *  SeGMO
 *
 *  Created by Dario Izzo on 6/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "GOproblem.h"
#include <math.h>
#include <vector>

GOProblem::GOProblem(int dimension, const std::vector<double> &lower, const std::vector <double> &upper) {
	setDimension(dimension);
	setBounds(lower, upper);
}

GOProblem::GOProblem(int dimension, const double *lower, const double *upper) {
	setDimension(dimension);
	setBounds(lower, upper);
}

void GOProblem::getBounds(std::vector<double> &lower, std::vector <double> &upper) const{
	lower = LB;
	upper = UB;
}

void GOProblem::setBounds(const std::vector<double> &lower, const std::vector <double> &upper){
	LB = lower;
	UB = upper;
}

void GOProblem::setBounds(const double *lower, const double *upper){
	LB.resize(dimension);
	UB.resize(dimension);
	for (int i = 0; i < dimension; i++) {
		LB[i] = lower[i];
		UB[i] = upper[i];
	}
}
