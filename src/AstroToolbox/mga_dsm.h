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

class mgadsmproblem {
public:
	mgadsmproblem():size(0) {};
	mgadsmproblem(int t, const int *seq, const size_t &size_, const double &AUdist_, const double &DVtotal_, const double &DVonboard_,
		const double &e_, const double &rp_):
		size(size_),type(t),sequence(seq,seq + size),e(e_),rp(rp_),AUdist(AUdist_),DVtotal(DVtotal_),
		DVonboard(DVonboard_),r(size),v(size),DV(size + 1) {
		for (size_t i = 0; i < size; ++i) {
			r[i] = new double[3];
			v[i] = new double[3];
		}
	}
	mgadsmproblem(const mgadsmproblem &m):size(m.size),type(m.type),sequence(m.sequence),e(m.e),rp(m.rp),asteroid(m.asteroid),
		AUdist(m.AUdist),DVtotal(m.DVtotal),DVonboard(m.DVonboard),r(m.size),v(m.size),DV(m.DV) {
		for (size_t i = 0; i < size; ++i) {
			r[i] = new double[3];
			v[i] = new double[3];
			for (size_t j = 0; j < 3; ++j) {
				r[i][j] = m.r[i][j];
				v[i][j] = m.v[i][j];
			}
		}
	}
	~mgadsmproblem() {
		for (size_t i = 0; i < size; ++i) {
			delete[] r[i];
			delete[] v[i];
		}
	}
	const size_t size;
	int type;							//problem type
	std::vector<int> sequence;				//fly-by sequence (ex: 3,2,3,3,5,is Earth-Venus-Earth-Earth-Jupiter)
	double e;							//insertion e (only in case total_DV_orbit_insertion)
	double rp;							//insertion rp in km (only in case total_DV_orbit_insertion)
	customobject asteroid;			//asteroid data (in case fly-by sequence has a final number = 10)
	double AUdist;						//Distance to reach in AUs (only in case of time2AUs)
	double DVtotal;						//Total DV allowed in km/s (only in case of time2AUs)
	double DVonboard;					//Total DV on the spacecraft in km/s (only in case of time2AUs)

	//Pre-allocated memory, in order to remove allocation of heap space in MGA_DSM calls
	mutable std::vector<double*> r;// = std::vector<double*>(n);
	mutable std::vector<double*> v;// = std::vector<double*>(n);
	mutable std::vector<double> DV;// = std::vector<double>(n+1);
};


int MGA_DSM(
			/* INPUT values: */
			std::vector<double> x ,	// it is the decision vector
			const mgadsmproblem &mgadsm,  // contains the problem specific data

			/* OUTPUT values: */
			double &J    // J output
			);

#endif
