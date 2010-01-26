/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
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

#ifndef PAGMO_PROBLEM_CLASSIC_H
#define PAGMO_PROBLEM_CLASSIC_H

#include <string>
#include <vector>

#include "../../config.h"
#include "base.h"
#include "classicobjfuns.h"

//***********************************************************************************
//Classical problems
//***********************************************************************************

namespace pagmo
{
namespace problem {

/// A test problem used to check/debug algorithms
/**
 * The objective function here is the sum of the components of the decision vector. Each component is limited to be
 * between 0 and 1. Trivially the best solution is 0 when all decision vector components are 0
 */
class __PAGMO_VISIBLE test : public base
{
	public:
		/// Constructor
		/**
		 * It instantiate a test problem. This is essentially the minimization of the vector norm defined
		* as the sum of the absolute values of its components. The problem
		* can be used to test algorithm performances and debug the algorithm code.
		 * \param[in] dim problem dimensions
		 */
		test(int dim);
		virtual test *clone() const {
			return new test(*this);
		}
		virtual std::string id_object() const;
	private:
		virtual double objfun_(const std::vector<double>& x) const {
			return testfunction(x);
		}
};	//end class testfunction

/// The Rastrigin Problem
/**
 * The objective function here is the Rastrigin function of arbitrary dimension. Bounds are set to -5.12,5.12
 */
class __PAGMO_VISIBLE rastrigin : public base
{
	public:
		/// Constructor
		/**
		 * It instantiate a Rastrigin problem
		 * \param[in] dim problem dimensions
		 */
		rastrigin(int dim);
		virtual rastrigin *clone() const {
			return new rastrigin(*this);
		}
		virtual std::string id_object() const;
	private:
		virtual double objfun_(const std::vector<double>& x) const {
			return rastriginf(x);
		}
};	//end class rastrigin

/// The Schwefel Problem
/**
 * The objective function here is the Schwefel function of arbitrary dimension. Bounds are set to -500,500
 */
class __PAGMO_VISIBLE schwefel : public base
{
	public:
		/// Constructor
		/**
		 * It instantiate a Schwefel problem
		 * \param[in] dim problem dimensions
		 */
		schwefel(int dim);
		virtual schwefel *clone() const {
			return new schwefel(*this);
		}
		virtual std::string id_object() const;
	private:
		virtual double objfun_(const std::vector<double>& x) const {
			return schwefelf(x);
		}
};	//end class schwefel

/// The Ackley Problem
/**
 * The objective function here is the Ackley function of arbitrary dimension. Bounds are set to -15,30
 */
class __PAGMO_VISIBLE ackley : public base
{
	public:
		/// Constructor
		/**
		 * It instantiate a Ackley problem
		 * \param[in] dim problem dimensions
		 */
		ackley(int dim);
		virtual ackley *clone() const {
			return new ackley(*this);
		}
		virtual std::string id_object() const;
	private:
		virtual double objfun_(const std::vector<double>& x) const {
			return ackleyf(x);
		}
};	//end class ackley

/// The Rosenbrock Problem
/**
 * \image html RosenbrockFunction.PNG
 *
 * The objective function here is the Rosenbrock function of arbitrary dimension. Bounds are set to -5,10
 */
class __PAGMO_VISIBLE rosenbrock : public base
{
	public:
		/// Constructor
		/**
		 * It instantiate a Rosenbrock problem
		 * \param[in] dim problem dimensions
		 */
		rosenbrock(int dim);
		virtual rosenbrock* clone() const {
			return new rosenbrock(*this);
		}
		virtual std::string id_object() const;
	private:
		virtual double objfun_(const std::vector<double>& x) const {
			return rosenbrockf(x);
		}
};	//end class rosenbrock

/// The Lennard-Jones Problem
/**
 * The objective function here is the Lennard-Jones potential in an ensemble of atoms.
 * Bounds on atoms positions are set in a hypercube of width 6; The online database http://physchem.ox.ac.uk/~doye/jon/structures/LJ/tables.150.html
 * contains all the best solution found up to large atoms numbers.
 */
class __PAGMO_VISIBLE lennardjones : public base
{
	public:
		/// Constructor
		/**
		 * It instantiate a Lennard-Jones problem
		 * \param[in] atoms Number of atoms. The problem dimension will correspondingly be 3*atoms - 6
		 */
		lennardjones(int atoms);
		virtual lennardjones *clone() const {
			return new lennardjones(*this);
		}
		virtual std::string id_object() const;
	private:
		virtual double objfun_(const std::vector<double>& x) const {
			return lennardjonesf(x);
		}
};	//end class lennardjones


/// The Levy Problem
/**
 * The objective function here is the Levy function of arbitrary dimension. Bounds are set to -10,10
 */
class __PAGMO_VISIBLE levy : public base
{
	public:
		/// Constructor
		/**
		 * It instantiate a Levy problem
		 * \param[in] dim problem dimensions
		 */
		levy(int dim);
		virtual std::string id_object() const;
	private:
		virtual levy *clone() const {
			return new levy(*this);
		}
		virtual double objfun_(const std::vector<double>& x) const {
			return levyf(x);
		}
};	//end class levy

/// The Griewank Problem
/**
 * The objective function here is the Griewank function of arbitrary dimension. Bounds are set to -600,600
 */
class __PAGMO_VISIBLE griewank : public base
{
	public:
		/// Constructor
		/**
		 * It instantiate a Griewank problem
		 * \param[in] dim problem dimensions
		 */
		griewank(int dim);
		virtual std::string id_object() const;
	private:
		virtual griewank *clone() const {
			return new griewank(*this);
		}
		virtual double objfun_(const std::vector<double>& x) const {
			return griewankf(x);
		}
};	//end class griewank

}
}

#endif
