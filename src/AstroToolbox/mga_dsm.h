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

#ifndef MGA_DSM_H
#define MGA_DSM_H

#include <vector>
#include "mga.h"

class mgadsmproblem {
public:
	mgadsmproblem():size(0),type(0),e(0),rp(0),AUdist(0),DVtotal(0),DVonboard(0) {};
	mgadsmproblem(int t, const int *seq, const size_t &size_, const double &AUdist_, const double &DVtotal_, const double &DVonboard_,
		const double &e_, const double &rp_):
		size(size_),type(t),sequence(seq,seq + size),e(e_),rp(rp_),AUdist(AUdist_),DVtotal(DVtotal_),
		DVonboard(DVonboard_),r(size),v(size),DV(size + 1),vrelin_vec(size-2) {
		for (size_t i = 0; i < size; ++i) {
			r[i] = new double[3];
			v[i] = new double[3];
			for (size_t j = 0; j < 3; ++j) {
				r[i][j] = 0;
				v[i][j] = 0;
			}
		}
	}
	mgadsmproblem(const mgadsmproblem &m):size(m.size),type(m.type),sequence(m.sequence),e(m.e),rp(m.rp),asteroid(m.asteroid),
		AUdist(m.AUdist),DVtotal(m.DVtotal),DVonboard(m.DVonboard),r(m.size),v(m.size),DV(m.DV),vrelin_vec(m.vrelin_vec) {
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
	int type;						//problem type
	std::vector<int> sequence;				//fly-by sequence (ex: 3,2,3,3,5,is Earth-Venus-Earth-Earth-Jupiter)
	double e;						//insertion e (only in case total_DV_orbit_insertion)
	double rp;						//insertion rp in km (only in case total_DV_orbit_insertion)
	customobject asteroid;					//asteroid data (in case fly-by sequence has a final number = 10)
	double AUdist;						//Distance to reach in AUs (only in case of time2AUs)
	double DVtotal;						//Total DV allowed in km/s (only in case of time2AUs)
	double DVonboard;					//Total DV on the spacecraft in km/s (only in case of time2AUs)

	//Pre-allocated memory, in order to remove allocation of heap space in MGA_DSM calls
	mutable std::vector<double*> r;		// = std::vector<double*>(n);
	mutable std::vector<double*> v;		// = std::vector<double*>(n);
	mutable std::vector<double> DV;		// = std::vector<double>(n+1);
	mutable std::vector<double> vrelin_vec;	// = std::vector<double>(n+1);
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int version) {
		ar & const_cast<size_t &>(size);
		ar & type;
		ar & sequence;
		ar & e;
		ar & rp;
		ar & asteroid;
		ar & AUdist;
		ar & DVtotal;
		ar & DVonboard;
		ar & DV;
		ar & vrelin_vec;
		boost::serialization::split_member(ar, *this, version);	// spliting the serialization into save/load to handle "r" and "v" vectors of pointers
	}
	template<class Archive>
	void save(Archive & ar, const unsigned int) const {
		std::vector<std::vector<double> > r_values(size,std::vector<double>(3));
		std::vector<std::vector<double> > v_values(size,std::vector<double>(3));
		for (int i = 0; i < (int)size; ++i) {
			for(int j = 0; j < 3; ++j){
				r_values[(size_t)i][j] = r[i][j];
			}
		}
		for (int i = 0; i < (int)size; ++i) {		
			for(int j = 0; j < 3; ++j){
				v_values[(size_t)i][j] = v[i][j];
			}
		}
		ar << static_cast<const std::vector<std::vector<double> > &>(r_values);
		ar << static_cast<const std::vector<std::vector<double> > &>(v_values);
	}
	template<class Archive>
	void load(Archive & ar, const unsigned int) {
		std::vector<std::vector<double> > r_values;
		std::vector<std::vector<double> > v_values;
		ar >> r_values;
		ar >> v_values;

		for (int i = 0; i < (int)size; ++i) {
			for(int j = 0; j < 3; ++j){
				r[(size_t)i][j] = r_values[i][j];
			}
		}
		for (int i = 0; i < (int)size; ++i) {
			for(int j = 0; j < 3; ++j){
				v[(size_t)i][j] = v_values[i][j];
			}
		}
	}
};


int MGA_DSM(
			/* INPUT values: */
			const std::vector<double> &x ,	// it is the decision vector
			const mgadsmproblem &mgadsm,  // contains the problem specific data

			/* OUTPUT values: */
			double &J    // J output
			);

#endif
