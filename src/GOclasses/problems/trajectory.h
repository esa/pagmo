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

#ifndef PAGMO_PROBLEM_TRAJECTORY_H
#define PAGMO_PROBLEM_TRAJECTORY_H

#include <boost/scoped_ptr.hpp>
#include <iostream>
#include <string>
#include <vector>

#include "../../AstroToolbox/mga_dsm.h"
#include "../../AstroToolbox/misc4Tandem.h"
#include "../../config.h"
#include "../problems/base.h"
#include "trajobjfuns.h"

namespace pagmo
{
namespace problem {

//***********************************************************************************
//Trajectory problems MGA
//***********************************************************************************
/// Cassini Multiple Gravity Assist interplanetary trajectory problem (from the GTOP database)
/**
  * This is a rather simple six dimensional MGA problem that is related to the Cassini spacecraft trajectory design problem. See
  * (http://www.esa.int/gsp/ACT/inf/op/globopt/evvejs.htm) for further informations
 *
*/
class __PAGMO_VISIBLE cassini1 : public base
{
	public:
		/// Constructor
		/**
		* It instantiate the cassini1 problem as defined in the ESA GTOP database (http://www.esa.int/gsp/ACT/inf/op/globopt/evvejs.htm)
		*/
		cassini1();
		virtual cassini1 *clone() const {
			return new cassini1(*this);
		}
		virtual std::string id_object() const {
			return id_name();
		}
	private:
		virtual double objfun_(const std::vector<double>&) const;
		static const double lb[6];
		static const double ub[6];
};	//end class cassini1

/// GTOC1 Multiple Gravity Assist interplanetary trajectory problem (from the GTOP database)
/**
 * The problem is part of the European Space Agency GTOP database (http://www.esa.int/gsp/ACT/inf/op/globopt/evvejs.htm)
 * This problem draws inspiration from the first edition of the Global Trajectory Optimisation Competition (GTOC1)
 * It is an 8 dimensional MGA problem with a rather long fly-by sequence including mainly Earth and Venus. The final target
 * is the asteroid TW229. The objective of the mission is to maximise the change in sami-major axis of the
 * asteroid orbit following an anaelastic impact of the spacecraft with the asteroid.
 */
class gtoc1 : public base
{
	public:
		/// Constructor
		/**
		* It instantiate the gtoc problem as defined in the ESA GTOP database (http://www.esa.int/gsp/ACT/inf/op/globopt/evevejsa.htm)
		*/
		gtoc1();
		virtual gtoc1 *clone() const {
			return new gtoc1(*this);
		}
		virtual std::string id_object() const {
			return id_name();
		}
	private:
		virtual double objfun_(const std::vector<double>&) const ;
		static const double lb[8];
		static const double ub[8];
};	//end class gtoc1



//***********************************************************************************
//Trajectory problems MGA-1DSM
//***********************************************************************************

class __PAGMO_VISIBLE messenger : public base
{
	public:
		messenger();
		virtual ~messenger() {}
		virtual messenger *clone() const {
			return new messenger(*this);
		}
		virtual std::string id_object() const {
			return id_name();
		}
	private:
		virtual double objfun_(const std::vector<double>&) const;
		mgadsmproblem mgadsm;
		static const double lb[18];
		static const double ub[18];
		static const int sequence[5];
};	//end class messenger

class __PAGMO_VISIBLE messengerfull : public base
{
	public:

		messengerfull();
		virtual ~messengerfull() {};
		virtual messengerfull *clone() const {
			return new messengerfull(*this);
		}
		virtual std::string id_object() const {
			return id_name();
		}
	private:
		virtual double objfun_(const std::vector<double>&) const;
		mgadsmproblem mgadsm;
		static const double lb[26];
		static const double ub[26];
		static const int sequence[7];
};	//end class messengerfull

/// Unconstrained TandEM trajectory problem (from the GTOP database)
/**
 * This interplanetary trajectory problem has 25 different instances, depending on the fly-by sequence adopted.
 * Please refer to http://www.esa.int/gsp/ACT/inf/op/globopt/TandEM.htm to select the proper instance. A classical
 * choice is 6, which corrsponds to an EVEES fly-by sequence. The possibility of adding one additional constraint
 * on the time-of-flight is given. If such a constraint is specified the components x[4]-x[7] are no longer time of flights
 * but they represent the percentage of the remaining time of flight to be used in one leg.
 */
class __PAGMO_VISIBLE tandem : public base
{
	public:
		/// Constructor
		/**
		* It instantiates an unconstrained TandEM problem.
		* \param[in] problemid This is an integer number from 1 to 24 encoding the fly-by sequence to be used. Please Check
		* http://www.esa.int/gsp/ACT/inf/op/globopt/TandEM.htm for more information
		*/
		tandem(const int problemid);

		/// Constructor
		/**
		* It instantiates a TandEM problem with a time constraint on the maximum time-of-flight.
		* \param[in] problemid This is an integer number from 1 to 24 encoding the fly-by sequence to be used. Please Check
		* http://www.esa.int/gsp/ACT/inf/op/globopt/TandEM.htm for more information
		* \param[in] tof (in years) This is a number setting the constraint on the total time of flight (10 from the GTOP database)
		*/
		tandem(const int problemid, const double tof);
		virtual ~tandem() {}
		virtual tandem *clone() const {
			return new tandem(*this);
		}
		virtual std::string id_object() const {
			return id_name();
		}
	private:
		virtual std::ostream &print(std::ostream &) const;
		virtual double objfun_(const std::vector<double>&) const;
		mgadsmproblem mgadsm;
		static const double lbunc[18];
		static const double ubunc[18];
		static const double lbcon[18];
		static const double ubcon[18];
		static const int Data[24][5];
		const double tof;
		mutable vector<double> copy_of_x;
};	//end class tandem


class cassini2 : public base
{
	public:
		cassini2();
		virtual ~cassini2() {}
		virtual cassini2 *clone() const {
			return new cassini2(*this);
		}
		virtual std::string id_object() const {
			return id_name();
		}
	private:
		virtual double objfun_(const std::vector<double>&) const;
		mgadsmproblem mgadsm;
		static const double lb[22];
		static const double ub[22];
		static const int sequence[6];
};	//end class cassini2

class rosetta : public base
{
	public:
		rosetta();
		virtual ~rosetta() {}
		virtual rosetta *clone() const {
			return new rosetta(*this);
		}
		virtual std::string id_object() const {
			return id_name();
		}
	private:
		virtual double objfun_(const std::vector<double>&) const;
		mgadsmproblem mgadsm;
		static const double lb[22];
		static const double ub[22];
		static const int sequence[6];
};	//end class rosetta

class sagas : public base
{
	public:
		sagas();
		virtual ~sagas() {}
		virtual sagas *clone() const {
			return new sagas(*this);
		}
		virtual std::string id_object() const {
			return id_name();
		}
	private:
		virtual double objfun_(const std::vector<double>&) const;
		mgadsmproblem mgadsm;
		static const double lb[12];
		static const double ub[12];
		static const int sequence[3];
};	//end class sagas

class __PAGMO_VISIBLE laplace : public base
{
	public:
		laplace(const std::vector<int> &);
		laplace(const laplace &);
		virtual ~laplace() {}
		virtual laplace *clone() const {
			return new laplace(*this);
		}
		std::string solution(const std::vector<double> &) const;
		virtual std::string id_object() const {
			return id_name();
		}
	private:
		virtual std::ostream &print(std::ostream &) const;
		void operator=(const laplace &) {}
		virtual double objfun_(const std::vector<double>&) const;
		boost::scoped_ptr<mgadsmproblem> mgadsm;
};

}
}

#endif
