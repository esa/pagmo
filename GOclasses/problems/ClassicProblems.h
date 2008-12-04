/*
 *  GOProblem.h
 *  SeGMO
 *
 *  Created by Dario Izzo on 6/4/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CLASSICPROBLEMS_H
#define CLASSICPROBLEMS_H

#include <vector>
#include "classicobjfuns.h"
#include "GOproblem.h"


//***********************************************************************************
//Classical problems
//***********************************************************************************


class TestProb : public GOProblem {
public:
	TestProb(int dim);
	virtual double objfun(const std::vector<double>& x) { return testfunction(x); }
};	//end class testfunctionProb

class rastriginProb : public GOProblem{
public:
	rastriginProb(int dim);
	virtual double objfun(const std::vector<double>& x) { return rastrigin(x); }
};	//end class rastriginProb

class schwefelProb : public GOProblem{
public:
	schwefelProb(int dim);
	virtual double objfun(const std::vector<double>& x) { return schwefel(x); }
};	//end class schwefelProb

class ackleyProb : public GOProblem{
public:
	ackleyProb(int dim);
	virtual double objfun(const std::vector<double>& x) { return ackley(x); }
};	//end class ackleyProb

class rosenbrockProb : public GOProblem{
public:
	rosenbrockProb(int dim);
	virtual double objfun(const std::vector<double>& x) { return rosenbrock(x); }
};	//end class rosenbrockProb

class lennardjonesProb : public GOProblem{
public:
	lennardjonesProb(int atoms);
	virtual double objfun(const std::vector<double>& x) { return lennardjones(x); };
};	//end class lennardjonesProb

class levyProb : public GOProblem{
public:
	levyProb(int dim);
	virtual double objfun(const std::vector<double>& x) { return levy(x); };
};	//end class levyProb


#endif
