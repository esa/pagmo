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

#ifndef TRAJECTORYPROBLEMS_H
#define TRAJECTORYPROBLEMS_H

#include <boost/scoped_ptr.hpp>
#include <iostream>
#include <string>
#include <vector>

#include "../../config.h"
#include "GOproblem.h"
#include "mga_dsm.h"
#include "misc4Tandem.h"
#include "trajobjfuns.h"

//***********************************************************************************
//Trajectory problems MGA
//***********************************************************************************
/// Cassini Multiple Gravity Assist interplanetary trajectory problem
/**
 * The problem is part of the European Space Agency GTOP database (http://www.esa.int/gsp/ACT/inf/op/globopt/evvejs.htm)
 * This is a rather simple six dimensional MGA problem that is related to the Cassini spacecraft trajectory design problem.
 * The objective of this mission is to reach Saturn and to be captured by its gravity into
 * an orbit having pericenter radius r_p=108950 km, and eccentricity e=0.98.
 * The planetary fly-by sequence considered is Earth, Venus, Venus, Earth, Jupiter, Saturn (as the one used by Cassini spacecraft).
 * As objective function we use the total deltaV accumulated during the mission, including the
 * launch deltaV and the various deltaV one needs to give at the planets and upon arrival to perform the
 * final orbit injection
*/
class __PAGMO_VISIBLE cassini1Prob : public GOProblem {
public:
       /// Constructor
       /**
        * It instantiate the cassini1 problem as defined in the ESA GTOP database (http://www.esa.int/gsp/ACT/inf/op/globopt/evvejs.htm)
        */
	cassini1Prob();
	virtual cassini1Prob *clone() const {return new cassini1Prob(*this);}
	virtual std::string id_object() const {return id_name(); }
private:
	virtual double objfun_(const std::vector<double>&) const;
	static const double lb[6];
	static const double ub[6];
};	//end class cassini1Prob

/// GTOC1 Multiple Gravity Assist interplanetary trajectory problem
/**
 * The problem is part of the European Space Agency GTOP database (http://www.esa.int/gsp/ACT/inf/op/globopt/evvejs.htm)
 * This problem draws inspiration from the first edition of the Global Trajectory Optimisation Competition (GTOC1)
 * It is an 8 dimensional MGA problem with a rather long fly-by sequence including mainly Earth and Venus. The final target
 * is the asteroid TW229. The objective of the mission is to maximise the change in sami-major axis of the
 * asteroid orbit following an anaelastic impact of the spacecraft with the asteroid.
 */
class gtoc1Prob : public GOProblem {
public:
       /// Constructor
       /**
        * It instantiate the gtoc problem as defined in the ESA GTOP database (http://www.esa.int/gsp/ACT/inf/op/globopt/evevejsa.htm)
        */
	gtoc1Prob();
	virtual gtoc1Prob *clone() const {return new gtoc1Prob(*this);}
	virtual std::string id_object() const {return id_name(); }
private:
	virtual double objfun_(const std::vector<double>&) const ;
	static const double lb[8];
	static const double ub[8];
};	//end class gtoc1Prob



//***********************************************************************************
//Trajectory problems MGA-1DSM
//***********************************************************************************

class __PAGMO_VISIBLE messengerProb : public GOProblem {
public:
	messengerProb();
        virtual ~messengerProb() {}
	virtual messengerProb *clone() const {return new messengerProb(*this);}
	virtual std::string id_object() const {return id_name(); }
private:
	virtual double objfun_(const std::vector<double>&) const;
	mgadsmproblem mgadsm;
	static const double lb[18];
	static const double ub[18];
	static const int sequence[5];
};	//end class messengerProb

class __PAGMO_VISIBLE messengerfullProb : public GOProblem {
public:

	messengerfullProb();
	virtual ~messengerfullProb() {};
	virtual messengerfullProb *clone() const {return new messengerfullProb(*this);}
	virtual std::string id_object() const {return id_name(); }
private:
	virtual double objfun_(const std::vector<double>&) const;
	mgadsmproblem mgadsm;
	static const double lb[26];
	static const double ub[26];
	static const int sequence[7];
};	//end class messengerfullProb

class tandemProb : public GOProblem {
public:
	tandemProb();
        virtual ~tandemProb() {}
	virtual tandemProb *clone() const {return new tandemProb(*this);}
	virtual std::string id_object() const {return id_name(); }
private:
	virtual double objfun_(const std::vector<double>&) const;
	mgadsmproblem mgadsm;
	static const double lb[18];
	static const double ub[18];
	static const int sequence[5];
};	//end class tandemProb


class cassini2Prob : public GOProblem {
public:
	cassini2Prob();
        virtual ~cassini2Prob() {}
	virtual cassini2Prob *clone() const {return new cassini2Prob(*this);}
	virtual std::string id_object() const {return id_name(); }
private:
	virtual double objfun_(const std::vector<double>&) const;
	mgadsmproblem mgadsm;
	static const double lb[22];
	static const double ub[22];
	static const int sequence[6];
};	//end class cassini2Prob

class rosettaProb : public GOProblem {
public:
	rosettaProb();
	virtual ~rosettaProb() {};
	virtual rosettaProb *clone() const {return new rosettaProb(*this);}
	virtual std::string id_object() const {return id_name(); }
private:
	virtual double objfun_(const std::vector<double>&) const;
	mgadsmproblem mgadsm;
	static const double lb[22];
	static const double ub[22];
	static const int sequence[6];
};	//end class rosettaProb

class sagasProb : public GOProblem {
public:
	sagasProb();
	virtual ~sagasProb() {};
	virtual sagasProb *clone() const {return new sagasProb(*this);}
	virtual std::string id_object() const {return id_name(); }
private:
	virtual double objfun_(const std::vector<double>&) const;
	mgadsmproblem mgadsm;
	static const double lb[12];
	static const double ub[12];
	static const int sequence[3];
};	//end class sagasProb

class __PAGMO_VISIBLE laplaceProb : public GOProblem {
	public:
		laplaceProb(const std::vector<int> &);
		laplaceProb(const laplaceProb &);
		virtual ~laplaceProb() {};
		virtual laplaceProb *clone() const {return new laplaceProb(*this);}
		std::string solution(const std::vector<double> &) const;
		virtual std::string id_object() const {return id_name(); }
	private:
		virtual std::ostream &print(std::ostream &) const;
		void operator=(const laplaceProb &) {};
		virtual double objfun_(const std::vector<double>&) const;
		boost::scoped_ptr<mgadsmproblem> mgadsm;
};

#endif
