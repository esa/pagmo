/*
 *  GOProblem.h
 *  SeGMO
 *
 *  Created by Dario Izzo on 6/4/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef TRAJECTORYPROBLEMS_H
#define TRAJECTORYPROBLEMS_H

#include <vector>
#include "trajobjfuns.h"
#include "GOproblem.h"
#include "mga_dsm.h"
#include "misc4Tandem.h"

//***********************************************************************************
//Trajectory problems MGA
//***********************************************************************************

class cassini1Prob : public GOProblem {
public:
	cassini1Prob();
	virtual double objfun(const std::vector<double>&);
};	//end class cassini1Prob

class gtoc1Prob : public GOProblem {
public:
	gtoc1Prob();
	virtual double objfun(const std::vector<double>&);
};	//end class gtoc1Prob



//***********************************************************************************
//Trajectory problems MGA-1DSM
//***********************************************************************************

class messengerProb : public GOProblem {
public:
	messengerProb();
	virtual ~messengerProb();
	virtual double objfun(const std::vector<double>&);
	
private:
	mgadsmproblem mgadsm;
};	//end class messengerProb

class messengerfullProb : public GOProblem {
public:

	messengerfullProb();
	virtual ~messengerfullProb();
	virtual double objfun(const std::vector<double>&);
	
private:
	mgadsmproblem mgadsm;
	static const double lb[26];
	static const double ub[26];
};	//end class messengerfullProb

class tandemProb : public GOProblem {
public:
	tandemProb();
	virtual ~tandemProb();
	virtual double objfun(const std::vector<double>&);
private:
	mgadsmproblem mgadsm;
};	//end class tandemProb


class cassini2Prob : public GOProblem {
public:
	cassini2Prob();
	virtual ~cassini2Prob();
	virtual double objfun(const std::vector<double>&);
private:
	mgadsmproblem mgadsm;
};	//end class cassini2Prob

class rosettaProb : public GOProblem {
public:
	rosettaProb();
	virtual ~rosettaProb();
	virtual double objfun(const std::vector<double>&);
private:
	mgadsmproblem mgadsm;
};	//end class rosettaProb

class sagasProb : public GOProblem {
public:
	sagasProb();
	virtual ~sagasProb();
	virtual double objfun(const std::vector<double>&);
private:
	mgadsmproblem mgadsm;
};	//end class sagasProb

#endif
