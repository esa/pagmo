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

#include "../../../config.h"
#include "../../Functions/objfuns/classicobjfuns.h"
#include "GOproblem.h"


//***********************************************************************************
//Classical problems
//***********************************************************************************

/// A test problem used to check/debug algorithms
/**
 * The objective function here is the sum of the components of the decision vector. Each component is limited to be
 * between 0 and 1. Trivially the best solution is 0 when all decision vector components are 0
 */
class __PAGMO_VISIBLE testProb : public GOProblem {
public:
        /// Constructor
       /**
        * It instantiate a test problem
        * \param[in] dim problem dimensions
        */
        testProb(int dim);
        virtual testProb *clone() const {return new testProb(*this);}
	virtual std::string id_object() const;
private:
	virtual double objfun_(const std::vector<double>& x) const { return testfunction(x); }
};	//end class testfunctionProb

/// The Rastrigin Problem
/**
 * The objective function here is the Rastrigin function of arbitrary dimension. Bounds are set to -5.12,5.12
 */
class __PAGMO_VISIBLE rastriginProb : public GOProblem{
public:
       /// Constructor
       /**
        * It instantiate a Rastrigin problem
        * \param[in] dim problem dimensions
        */
	rastriginProb(int dim);
	virtual rastriginProb *clone() const {return new rastriginProb(*this);}
	virtual std::string id_object() const;
private:
	virtual double objfun_(const std::vector<double>& x) const { return rastrigin(x); }
};	//end class rastriginProb

/// The Schwefel Problem
/**
 * The objective function here is the Schwefel function of arbitrary dimension. Bounds are set to -500,500
 */
class __PAGMO_VISIBLE schwefelProb : public GOProblem{
public:
        /// Constructor
       /**
        * It instantiate a Schwefel problem
        * \param[in] dim problem dimensions
        */
	schwefelProb(int dim);
	virtual schwefelProb *clone() const {return new schwefelProb(*this);}
	virtual std::string id_object() const;
private:
	virtual double objfun_(const std::vector<double>& x) const { return schwefel(x); }
};	//end class schwefelProb

/// The Ackley Problem
/**
 * The objective function here is the Ackley function of arbitrary dimension. Bounds are set to -15,30
 */
class __PAGMO_VISIBLE ackleyProb : public GOProblem{
public:
        /// Constructor
       /**
        * It instantiate a Ackley problem
        * \param[in] dim problem dimensions
        */
	ackleyProb(int dim);
	virtual ackleyProb *clone() const {return new ackleyProb(*this);}
	virtual std::string id_object() const;
private:
	virtual double objfun_(const std::vector<double>& x) const { return ackley(x); }
};	//end class ackleyProb

/// The Rosenbrock Problem
/**
 * \image html RosenbrockFunction.PNG
 *
 * The objective function here is the Rosenbrock function of arbitrary dimension. Bounds are set to -5,10
 */
class __PAGMO_VISIBLE rosenbrockProb : public GOProblem{
public:
        /// Constructor
       /**
        * It instantiate a Rosenbrock problem
        * \param[in] dim problem dimensions
        */
	rosenbrockProb(int dim);
	virtual rosenbrockProb* clone() const {return new rosenbrockProb(*this);}
	virtual std::string id_object() const;
private:
	virtual double objfun_(const std::vector<double>& x) const { return rosenbrock(x); }
};	//end class rosenbrockProb

/// The Lennard-Jones Problem
/**
 * The objective function here is the Lennard-Jones potential in an ensemble of atoms.
 * Bounds on atoms positions are set in a hypercube of width 6; The online database http://physchem.ox.ac.uk/~doye/jon/structures/LJ/tables.150.html
 * contains all the best solution found up to large atoms numbers.
 */
class __PAGMO_VISIBLE lennardjonesProb : public GOProblem{
public:
       /// Constructor
       /**
        * It instantiate a Lennard-Jones problem
        * \param[in] atoms Number of atoms. The problem dimension will correspondingly be 3*atoms - 6
        */
	lennardjonesProb(int atoms);
	virtual lennardjonesProb *clone() const {return new lennardjonesProb(*this);}
	virtual std::string id_object() const;
private:
	virtual double objfun_(const std::vector<double>& x) const { return lennardjones(x); }
};	//end class lennardjonesProb


/// The Levy Problem
/**
 * The objective function here is the Levy function of arbitrary dimension. Bounds are set to -10,10
 */
class __PAGMO_VISIBLE levyProb : public GOProblem{
public:
       /// Constructor
       /**
        * It instantiate a Levy problem
        * \param[in] dim problem dimensions
        */
	levyProb(int dim);
	virtual std::string id_object() const;
private:
	virtual levyProb *clone() const {return new levyProb(*this);}
	virtual double objfun_(const std::vector<double>& x) const { return levy(x); }
};	//end class levyProb

/// The Griewank Problem
/**
 * The objective function here is the Griewank function of arbitrary dimension. Bounds are set to -600,600
 */
class __PAGMO_VISIBLE griewankProb : public GOProblem{
public:
       /// Constructor
       /**
        * It instantiate a Griewank problem
        * \param[in] dim problem dimensions
        */
        griewankProb(int dim);
        virtual std::string id_object() const;
private:
        virtual griewankProb *clone() const {return new griewankProb(*this);}
        virtual double objfun_(const std::vector<double>& x) const { return griewank(x); }
};	//end class griewankProb


#endif
