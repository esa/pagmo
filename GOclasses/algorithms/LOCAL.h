/*
 *  LOCAL.h
 *  SeGMO
 *
 *  Created by Dario Izzo on 6/10/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef LOCAL_H
#define LOCAL_H

#include "individual.h"
#include <vector>
#include <math.h>


//The following class defines a simple local optimisation algorithm gradient free
class SIMPLELOCALalgorithm{
public:

Individual evolve(Individual x0, double(*objfun)(const std::vector<double>&), std::vector<double> LB, std::vector<double> UB );

void initSIMPLELOCAL(double rangeInit,
					 double reduxCoeffInit,
					 double minrangeInit);

private:
	 double range;
	 double reduxCoeff;
	 double minrange;
};

#endif
