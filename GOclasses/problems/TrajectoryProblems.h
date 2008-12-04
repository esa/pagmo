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
//Trajectory problems
//***********************************************************************************

class messengerProb : public GOProblem {
public:
	messengerProb();
	virtual ~messengerProb();
	virtual double objfun(const std::vector<double>&);
	
private:
	mgadsmproblem mgadsm;
	double Delta_V[6]; //Dummy array, required for calling MGA_DSM
};	//end class messengerProb

class tandemProb : public GOProblem {
public:
	tandemProb();
	virtual ~tandemProb();
	virtual double objfun(const std::vector<double>&);
private:
	mgadsmproblem mgadsm;
	double Delta_V[6]; //Dummy array, required for calling MGA_DSM
};	//end class tandemProb

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

class cassini2Prob : public GOProblem {
public:
	cassini2Prob();
	virtual ~cassini2Prob();
	virtual double objfun(const std::vector<double>&);
private:
	mgadsmproblem mgadsm;
	double Delta_V[7]; //Dummy array, required for calling MGA_DSM
};	//end class cassini1Prob

class rosettaProb : public GOProblem {
public:
	rosettaProb();
	virtual ~rosettaProb();
	virtual double objfun(const std::vector<double>&);
private:
	double Delta_V[7];
	mgadsmproblem mgadsm;
};	//end class rosettaProb

class sagasProb : public GOProblem {
public:
	sagasProb();
	virtual ~sagasProb();
	virtual double objfun(const std::vector<double>&);
private:
	double Delta_V[4];
	mgadsmproblem mgadsm;
};	//end class sagasProb

#endif
