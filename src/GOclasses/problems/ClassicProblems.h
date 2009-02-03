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

// 04/06/08 Created by Dario Izzo.

#ifndef PAGMO_CLASSICPROBLEMS_H
#define PAGMO_CLASSICPROBLEMS_H

#include <vector>

#include "../../config.h"
#include "classicobjfuns.h"
#include "GOproblem.h"


//***********************************************************************************
//Classical problems
//***********************************************************************************


class __PAGMO_VISIBLE TestProb : public GOProblem {
public:
	TestProb(int dim);
	virtual double objfun(const std::vector<double>& x) const { return testfunction(x); }
	virtual TestProb *clone() const {return new TestProb(*this);}
};	//end class testfunctionProb

class __PAGMO_VISIBLE rastriginProb : public GOProblem{
public:
	rastriginProb(int dim);
	virtual double objfun(const std::vector<double>& x) const { return rastrigin(x); }
	virtual rastriginProb *clone() const {return new rastriginProb(*this);}
};	//end class rastriginProb

class __PAGMO_VISIBLE schwefelProb : public GOProblem{
public:
	schwefelProb(int dim);
	virtual double objfun(const std::vector<double>& x) const { return schwefel(x); }
	virtual schwefelProb *clone() const {return new schwefelProb(*this);}
};	//end class schwefelProb

class __PAGMO_VISIBLE ackleyProb : public GOProblem{
public:
	ackleyProb(int dim);
	virtual double objfun(const std::vector<double>& x) const { return ackley(x); }
	virtual ackleyProb *clone() const {return new ackleyProb(*this);}
};	//end class ackleyProb

class __PAGMO_VISIBLE rosenbrockProb : public GOProblem{
public:
	rosenbrockProb(int dim);
	virtual double objfun(const std::vector<double>& x) const { return rosenbrock(x); }
	virtual rosenbrockProb* clone() const {return new rosenbrockProb(*this);}
};	//end class rosenbrockProb

class __PAGMO_VISIBLE lennardjonesProb : public GOProblem{
public:
	lennardjonesProb(int atoms);
	virtual double objfun(const std::vector<double>& x) const { return lennardjones(x); }
	virtual lennardjonesProb *clone() const {return new lennardjonesProb(*this);}
};	//end class lennardjonesProb

class __PAGMO_VISIBLE levyProb : public GOProblem{
public:
	levyProb(int dim);
	virtual double objfun(const std::vector<double>& x) const { return levy(x); };
	virtual levyProb *clone() const {return new levyProb(*this);}
};	//end class levyProb


#endif
