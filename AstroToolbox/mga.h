// ------------------------------------------------------------------------ //
// This source file is part of the 'ESA Advanced Concepts Team's			//
// Space Mechanics Toolbox' software.                                       //
//                                                                          //
// The source files are for research use only,                              //
// and are distributed WITHOUT ANY WARRANTY. Use them on your own risk.     //
//                                                                          //
// Copyright (c) 2004-2007 European Space Agency                            //
// ------------------------------------------------------------------------ //

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
