// ------------------------------------------------------------------------ //
// This source file is part of the 'ESA Advanced Concepts Team's			//
// Space Mechanics Toolbox' software.                                       //
//                                                                          //
// The source files are for research use only,                              //
// and are distributed WITHOUT ANY WARRANTY. Use them on your own risk.     //
//                                                                          //
// Copyright (c) 2004-2007 European Space Agency                            //
// ------------------------------------------------------------------------ //

#ifndef MGA_DSM_H
#define MGA_DSM_H

#include <vector>

#include "mga.h"

struct mgadsmproblem {
	int type;							//problem type
	std::vector<int> sequence;				//fly-by sequence (ex: 3,2,3,3,5,is Earth-Venus-Earth-Earth-Jupiter)
	double e;							//insertion e (only in case total_DV_orbit_insertion)
	double rp;							//insertion rp in km (only in case total_DV_orbit_insertion)
	customobject asteroid;			//asteroid data (in case fly-by sequence has a final number = 10)
	double AUdist;						//Distance to reach in AUs (only in case of time2AUs)
	double DVtotal;						//Total DV allowed in km/s (only in case of time2AUs)
	double DVonboard;					//Total DV on the spacecraft in km/s (only in case of time2AUs)

	//Pre-allocated memory, in order to remove allocation of heap space in MGA_DSM calls
	std::vector<double*> r;// = std::vector<double*>(n);
	std::vector<double*> v;// = std::vector<double*>(n);
	std::vector<double> DV;// = std::vector<double>(n+1);
};


int MGA_DSM(
			/* INPUT values: */
			std::vector<double> x ,	// it is the decision vector
			mgadsmproblem,  // contains the problem specific data

			/* OUTPUT values: */
			double& J,    // J output
			double* DVVec// DVVec
			);

#endif
