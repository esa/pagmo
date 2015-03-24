/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
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

#ifndef MISSION_H
#define MISSION_H

#include <vector>
#include "Pl_Eph_An.h"
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/split_member.hpp>

// problem types
#define orbit_insertion          0 // Tandem
#define total_DV_orbit_insertion 1 // Cassini 1
#define rndv                     2 // Rosetta
#define total_DV_rndv            3 // Cassini 2 and Messenger
#define asteroid_impact          4 // gtoc1
#define time2AUs                 5 // SAGAS 

struct customobject
{
	customobject()
	{
		for (int i = 0; i < 6; ++i) {
			keplerian[i] = 0;
		}
		epoch = 0;
		mu = 0;
	}
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int){
		ar & keplerian;
		ar & epoch;
		ar & mu;
	}
	double keplerian[6];
	double epoch;
	double mu;
};


struct mgaproblem {
	mgaproblem():type(0),e(0),rp(0),Isp(0),mass(0),DVlaunch(0) {}
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int){
		ar & type;
		ar & sequence;
		ar & rev_flag;
		ar & e;
		ar & rp;
		ar & asteroid;
		ar & Isp;
		ar & mass;
		ar & DVlaunch;
	}
	int type;							//problem type
	std::vector<int> sequence;				//fly-by sequence (ex: 3,2,3,3,5,is Earth-Venus-Earth-Earth-Jupiter)
	std::vector<int> rev_flag;				//vector of flags for clockwise legs
	double e;							//insertion e (only in case total_DV_orbit_insertion)
	double rp;							//insertion rp in km (only in case total_DV_orbit_insertion)
	customobject asteroid;				//asteroid data (in case fly-by sequence has a final number = 10)
	double Isp;
	double mass;
	double DVlaunch;
};

int MGA( 
		 //INPUTS
		 std::vector<double>,
		 mgaproblem, 
		
		 //OUTPUTS
		 std::vector <double>&, std::vector<double>&, double&); 

#endif
