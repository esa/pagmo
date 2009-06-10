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

#ifndef MISSION_H
#define MISSION_H

#include <vector>
#include "Pl_Eph_An.h"

using namespace std;

// problem types
#define orbit_insertion          0 // Tandem
#define total_DV_orbit_insertion 1 // Cassini 1
#define rndv                     2 // Rosetta
#define total_DV_rndv            3 // Cassini 2 and Messenger
#define asteroid_impact          4 // gtoc1
#define time2AUs                 5 // SAGAS 

struct customobject
{
	double keplerian[6];
	double epoch;
	double mu;
};


struct mgaproblem {
	int type;							//problem type
	vector<int> sequence;				//fly-by sequence (ex: 3,2,3,3,5,is Earth-Venus-Earth-Earth-Jupiter)
	vector<int> rev_flag;				//vector of flags for clockwise legs
	double e;							//insertion e (only in case total_DV_orbit_insertion)
	double rp;							//insertion rp in km (only in case total_DV_orbit_insertion)
	customobject asteroid;				//asteroid data (in case fly-by sequence has a final number = 10)
	double Isp;
	double mass;
	double DVlaunch;
};

int MGA( 
		 //INPUTS
		 vector<double>,
		 mgaproblem, 
		
		 //OUTPUTS
		 vector <double>&, vector<double>&, double&); 

#endif
