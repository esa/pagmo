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

#include "GOproblem.h"
#include "mga_dsm.h"
#include "misc4Tandem.h"
#include "trajobjfuns.h"

//***********************************************************************************
//Trajectory problems MGA
//***********************************************************************************

class cassini1Prob : public GOProblem {
public:
	cassini1Prob();
	virtual double objfun(const std::vector<double>&) const;
	virtual cassini1Prob *clone() const {return new cassini1Prob(*this);}
private:
	static const double lb[6];
	static const double ub[6];
};	//end class cassini1Prob

class gtoc1Prob : public GOProblem {
public:
	gtoc1Prob();
	virtual double objfun(const std::vector<double>&) const ;
	virtual gtoc1Prob *clone() const {return new gtoc1Prob(*this);}
private:
	static const double lb[8];
	static const double ub[8];
};	//end class gtoc1Prob



//***********************************************************************************
//Trajectory problems MGA-1DSM
//***********************************************************************************

class messengerProb : public GOProblem {
public:
	messengerProb();
	virtual ~messengerProb() {};
	virtual double objfun(const std::vector<double>&) const;
	virtual messengerProb *clone() const {return new messengerProb(*this);}
private:
	mgadsmproblem mgadsm;
	static const double lb[18];
	static const double ub[18];
	static const int sequence[5];
};	//end class messengerProb

class messengerfullProb : public GOProblem {
public:

	messengerfullProb();
	virtual ~messengerfullProb() {};
	virtual double objfun(const std::vector<double>&) const;
	virtual messengerfullProb *clone() const {return new messengerfullProb(*this);}
private:
	mgadsmproblem mgadsm;
	static const double lb[26];
	static const double ub[26];
	static const int sequence[7];
};	//end class messengerfullProb

class tandemProb : public GOProblem {
public:
	tandemProb();
	virtual ~tandemProb() {};
	virtual double objfun(const std::vector<double>&) const;
	virtual tandemProb *clone() const {return new tandemProb(*this);}
private:
	mgadsmproblem mgadsm;
	static const double lb[18];
	static const double ub[18];
	static const int sequence[5];
};	//end class tandemProb


class cassini2Prob : public GOProblem {
public:
	cassini2Prob();
	virtual ~cassini2Prob() {};
	virtual double objfun(const std::vector<double>&) const;
	virtual cassini2Prob *clone() const {return new cassini2Prob(*this);}
private:
	mgadsmproblem mgadsm;
	static const double lb[22];
	static const double ub[22];
	static const int sequence[6];
};	//end class cassini2Prob

class rosettaProb : public GOProblem {
public:
	rosettaProb();
	virtual ~rosettaProb() {};
	virtual double objfun(const std::vector<double>&) const;
	virtual rosettaProb *clone() const {return new rosettaProb(*this);}
private:
	mgadsmproblem mgadsm;
	static const double lb[22];
	static const double ub[22];
	static const int sequence[6];
};	//end class rosettaProb

class sagasProb : public GOProblem {
public:
	sagasProb();
	virtual ~sagasProb() {};
	virtual double objfun(const std::vector<double>&) const;
	virtual sagasProb *clone() const {return new sagasProb(*this);}
private:
	mgadsmproblem mgadsm;
	static const double lb[12];
	static const double ub[12];
	static const int sequence[3];
};	//end class sagasProb

#endif
