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

// 11/06/08 Created by Dario Izzo.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "../../Functions/objfuns/trajobjfuns.h"
#include "../../constants.h"
#include "../../exceptions.h"
#include "base.h"
#include "trajectory.h"

//***********************************************************************************
//Trajectory problems
//***********************************************************************************

namespace pagmo
{
namespace problem {

const double messengerfull::lb[26] = {1900, 3,    0, 0, 100, 100, 100, 100, 100, 100, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 1.1, 1.1, 1.05, 1.05, 1.05, -M_PI, -M_PI, -M_PI, -M_PI, -M_PI};
const double messengerfull::ub[26] = {2200, 4.05, 1, 1, 500, 500, 500, 500, 500, 550, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99,   6,   6,    6,    6,    6,  M_PI,  M_PI,  M_PI,  M_PI,  M_PI};
const int messengerfull::sequence[7] = {3, 2, 2, 1, 1, 1, 1};

messengerfull::messengerfull():base(26,lb,ub),mgadsm(orbit_insertion,sequence,7,0,0,0,0.704,2440 + 200) {}

double messengerfull::objfun_(const std::vector<double>& x) const
{
	double obj = 0.0;
	MGA_DSM(x, mgadsm,
	        obj);
	return obj;
}


const double messenger::lb[18] = {1000, 1, 0, 0, 200, 30,  30,  30,  0.01, 0.01, 0.01, 0.01, 1.1, 1.1, 1.1, -M_PI, -M_PI, -M_PI};
const double messenger::ub[18] = {4000, 5, 1, 1, 400, 400, 400, 400, 0.99, 0.99, 0.99, 0.99, 6,   6,   6,    M_PI,  M_PI,  M_PI};
const int messenger::sequence[5] = {3, 3, 2, 2, 1};

messenger::messenger():base(18,lb,ub),mgadsm(total_DV_rndv,sequence,5,0,0,0,0,0) {}

double messenger::objfun_(const std::vector<double>& x) const
{
	double obj = 0.0;
	MGA_DSM(x, mgadsm,
	        obj);
	return obj;
}

const double tandem::lbunc[18] = {5475, 2.5, 0, 0, 20   , 20  ,  20 , 20  , 0.01, 0.01, 0.01, 0.01, 1.05, 1.05, 1.05, -M_PI, -M_PI, -M_PI};
const double tandem::ubunc[18] = {9132, 4.9, 1, 1, 2500 , 2500, 2500, 2500, 0.99, 0.99, 0.99, 0.99,    10,    10,    10,  M_PI,  M_PI,  M_PI};

const double tandem::lbcon[18] = {5475, 2.5, 0, 0, 0, 0, 0, 0, 0.01, 0.01, 0.01, 0.01, 1.05, 1.05, 1.05, -M_PI, -M_PI, -M_PI};
const double tandem::ubcon[18] = {9132, 4.9, 1, 1, 1, 1, 1, 1, 0.99, 0.99, 0.99, 0.99,    10,    10,    10,  M_PI,  M_PI,  M_PI};


const int tandem::Data[24][5] = {
	{3,2,2,2,6},
	{3,2,2,3,6},
	{3,2,2,4,6},
	{3,2,2,5,6},
	{3,2,3,2,6},
	{3,2,3,3,6},
	{3,2,3,4,6},
	{3,2,3,5,6},
	{3,2,4,2,6},
	{3,2,4,3,6},
	{3,2,4,4,6},
	{3,2,4,5,6},
	{3,3,2,2,6},
	{3,3,2,3,6},
	{3,3,2,4,6},
	{3,3,2,5,6},
	{3,3,3,2,6},
	{3,3,3,3,6},
	{3,3,3,4,6},
	{3,3,3,5,6},
	{3,3,4,2,6},
	{3,3,4,3,6},
	{3,3,4,4,6},
	{3,3,4,5,6}
};

const int sequence[5] = {1,1,1,1,1};

tandem::tandem(const int probid):base(18,lbunc,ubunc), mgadsm(orbit_insertion,sequence,5,0,0,0,0.98531407996358,80330.0), tof(0), copy_of_x(18)
{
	if (probid < 1 || probid > 24) {
		pagmo_throw(value_error,"probid needs to be a number in [1,24]");
	}
	for (int i =0; i<5; i++) {
		mgadsm.sequence[i] = Data[probid-1][i];
	}
};

tandem::tandem(const int probid, const double tof_):base(18,lbcon,ubcon), mgadsm(orbit_insertion,sequence,5,0,0,0,0.98531407996358,80330.0), tof(tof_),copy_of_x(18)
{
	if (probid < 1 || probid > 24) {
		pagmo_throw(value_error,"probid needs to be an integer in [1,24]");
	}
	if (tof_ > 20 && tof_ <5) {
		pagmo_throw(value_error,"time of flight constraint needs to be a number in [5,20] (years)");
	}
	for (int i =0; i<5; i++) {
		mgadsm.sequence[i] = Data[probid-1][i];
	}
};



double tandem::objfun_(const std::vector<double>& x) const
{
	double obj = 0;
	if (tof!=0) { //constrained problem
		//Here we copy the chromosome into a new vector and we transform its time percentages into days
		copy_of_x = x;
		copy_of_x[4] = x[4]*365.25*tof;
		copy_of_x[5] = x[5]*(365.25*tof-copy_of_x[4]);
		copy_of_x[6] = x[6]*(365.25*tof-copy_of_x[4]-copy_of_x[5]);
		copy_of_x[7] = x[7]*(365.25*tof-copy_of_x[4]-copy_of_x[5]-copy_of_x[6]);
		MGA_DSM(copy_of_x, mgadsm, obj);
	} else {	//unconstrained problem
		MGA_DSM(x, mgadsm, obj);
	}
	//evaluating the mass from the dvs
	double rE[3];
	double vE[3];
	Planet_Ephemerides_Analytical (x[0],3,rE,vE);
	double VINFE = x[1];
	double udir = x[2];
	double vdir = x[3];
	double vtemp[3];
	vtemp[0]= rE[1]*vE[2]-rE[2]*vE[1];
	vtemp[1]= rE[2]*vE[0]-rE[0]*vE[2];
	vtemp[2]= rE[0]*vE[1]-rE[1]*vE[0];
	double iP1[3];
	double normvE=sqrt(vE[0]*vE[0]+vE[1]*vE[1]+vE[2]*vE[2]);
	iP1[0]=	vE[0]/normvE;
	iP1[1]=	vE[1]/normvE;
	iP1[2]=	vE[2]/normvE;
	double zP1[3];
	double normvtemp=sqrt(vtemp[0]*vtemp[0]+vtemp[1]*vtemp[1]+vtemp[2]*vtemp[2]);
	zP1[0]= vtemp[0]/normvtemp;
	zP1[1]= vtemp[1]/normvtemp;
	zP1[2]= vtemp[2]/normvtemp;
	double jP1[3];
	jP1[0]= zP1[1]*iP1[2]-zP1[2]*iP1[1];
	jP1[1]= zP1[2]*iP1[0]-zP1[0]*iP1[2];
	jP1[2]= zP1[0]*iP1[1]-zP1[1]*iP1[0];
	double theta=2*M_PI*udir; 		//See Picking a Point on a Sphere
	double phi=acos(2*vdir-1)-M_PI/2; //In this way: -pi/2<phi<pi/2 so phi can be used as out-of-plane rotation
	double vinf[3];
	vinf[0]=VINFE*(cos(theta)*cos(phi)*iP1[0]+sin(theta)*cos(phi)*jP1[0]+sin(phi)*zP1[0]);
	vinf[1]=VINFE*(cos(theta)*cos(phi)*iP1[1]+sin(theta)*cos(phi)*jP1[1]+sin(phi)*zP1[1]);
	vinf[2]=VINFE*(cos(theta)*cos(phi)*iP1[2]+sin(theta)*cos(phi)*jP1[2]+sin(phi)*zP1[2]);
	//We rotate it to the equatorial plane
	ecl2equ(vinf,vinf);
	//And we find the declination in degrees
	double normvinf=sqrt(vinf[0]*vinf[0]+vinf[1]*vinf[1]+vinf[2]*vinf[2]);
	double sindelta = vinf[2] / normvinf;
	double declination = asin(sindelta)/M_PI*180;

	//double m_initial = SoyuzFregat(VINFE,declination);
	double m_initial = Atlas501(VINFE,declination);

	//We evaluate the final mass
	double Isp = 312;
	double g0 = 9.80665;
	double sumDVvec=0;
	//double totaltime=x[4]+x[5]+x[6]+x[7];
	for (unsigned int i=1;i<=5;i++) {
		sumDVvec=sumDVvec+mgadsm.DV[i];
	}
	double m_final;
	sumDVvec=sumDVvec+0.165; //losses for 3 swgbys + insertion

	m_final = m_initial * exp(-sumDVvec/Isp/g0*1000);

	return -log(m_final);
}

std::ostream &tandem::print(std::ostream &s) const
{
	base::print(s);
	s << "Flyby sequence: " << '\n';
	for (size_t i = 0; i < 5; ++i) {
		s << mgadsm.sequence[i] << " ";
	}
	s << '\n';
	return s;
}


const double cassini1::lb[6] = {-1000, 30,100,30,400,1000};
const double cassini1::ub[6] = {0,400,470,400,2000,6000};

cassini1::cassini1():base(6,lb,ub) {}

double cassini1::objfun_(const std::vector<double>& x) const
{
	return cassini1f(x);
}

const double gtoc1::lb[8] = {3000,14,14,14,14,100,366,300};
const double gtoc1::ub[8] = {10000,2000,2000,2000,2000,9000,9000,9000};

gtoc1::gtoc1():base(8,lb,ub) {}

double gtoc1::objfun_(const std::vector<double>& x) const
{
	return gtoc1f(x);
}

const double cassini2::lb[22] = {-750, 3, 0, 0, 100, 100, 30, 400, 800, 0.01, 0.01, 0.01, 0.01, 0.01, 1.05, 1.05, 1.15, 1.7, -M_PI, -M_PI, -M_PI, -M_PI};
const double cassini2::ub[22] = {780,  5, 1, 1, 400, 500, 300, 1600, 2200, 0.9, 0.9, 0.9, 0.9, 0.9, 6, 6, 6.5, 291, M_PI, M_PI, M_PI, M_PI};
const int cassini2::sequence[6] = {3, 2, 2, 3, 5, 6};

cassini2::cassini2():base(22,lb,ub),mgadsm(total_DV_rndv,sequence,6,0,0,0,0,0) {}

double cassini2::objfun_(const std::vector<double>& x) const
{
	double obj = 0;
	MGA_DSM(x, mgadsm,
	        obj);
	return obj;
}

const double rosetta::lb[22] = {1460, 3, 0, 0, 300, 150, 150, 300, 700 , 0.01, 0.01, 0.01, 0.01, 0.01, 1.05, 1.05, 1.05, 1.05, -M_PI, -M_PI, -M_PI, -M_PI};
const double rosetta::ub[22] = {1825, 5, 1, 1, 500, 800, 800, 800, 1850, 0.9, 0.9 , 0.9 , 0.9 , 0.9 , 9   , 9   , 9   , 9   , M_PI , M_PI , M_PI , M_PI};
const int rosetta::sequence[6] = {3, 3, 4, 3, 3, 10};

rosetta::rosetta():base(22,lb,ub),mgadsm(rndv,sequence,6,0,0,0,0,0)
{
	//MGA_DSM stuff
	mgadsm.asteroid.keplerian[0] = 3.50294972836275;
	mgadsm.asteroid.keplerian[1] = 0.6319356;
	mgadsm.asteroid.keplerian[2] = 7.12723;
	mgadsm.asteroid.keplerian[3] = 50.92302;
	mgadsm.asteroid.keplerian[4] = 11.36788;
	mgadsm.asteroid.keplerian[5] = 0.0;
	mgadsm.asteroid.epoch = 52504.23754000012;
	mgadsm.asteroid.mu = 0.0;
};

double rosetta::objfun_(const std::vector<double>& x) const
{
	double obj = 0;
	MGA_DSM(x, mgadsm,
	        obj);
	return obj;
}

const double sagas::lb[12] = {7000, 0, 0, 0, 50, 300, 0.01, 0.01, 1.05, 8, -M_PI, -M_PI};
const double sagas::ub[12] = {9100, 7, 1, 1, 2000, 2000, 0.9, 0.9 ,7 ,500 ,M_PI,M_PI};
const int sagas::sequence[3] = {3,3,5};

sagas::sagas():base(12,lb,ub),mgadsm(time2AUs,sequence,3,50.0,6.782,1.782,0,0) {}

double sagas::objfun_(const std::vector<double>& x)  const
{
	double obj = 0;
	MGA_DSM(x, mgadsm,
	        obj);
	return obj;

}

laplace::laplace(const std::vector<int> &seq):base(4 * seq.size() - 2),mgadsm(0)
{
	const size_t seq_size = seq.size();
	if (seq_size < 2) {
		pagmo_throw(value_error,"flyby sequence size must be at least 2");
	}
	for (size_t i = 0; i < seq_size; ++i) {
		if (seq[i] < 2 || seq[i] > 5) {
			pagmo_throw(value_error,"invalid planet index in flyby sequence");
		}
	}
	mgadsm.reset(new mgadsmproblem(orbit_insertion,&seq[0],seq.size(),0,0,0,.97,4 * 71492));
	// Set bounds.
	LB[0] = 5475;
	UB[0] = 9132;
	LB[1] = 0.1;
	UB[1] = 3.5;
	LB[2] = LB[3] = 0;
	UB[2] = UB[3] = 1;
	for (size_t i = 0; i < seq_size - 1; ++i) {
		LB[i + 4] = 20;
		UB[i + 4] = 2500;
	}
	for (size_t i = 0; i < seq_size - 1; ++i) {
		LB[i + 3 + seq_size] = 0.01;
		UB[i + 3 + seq_size] = 0.99;
	}
	for (size_t i = 0; i < seq_size - 2; ++i) {
		LB[i + 2 * (seq_size + 1)] = 1.05;
		UB[i + 2 * (seq_size + 1)] = 100.0;
	}
	for (size_t i = 0; i < seq_size - 2; ++i) {
		LB[i + 3 * seq_size] = -M_PI;
		UB[i + 3 * seq_size] = M_PI;
	}
}

laplace::laplace(const laplace &l):base(l),mgadsm(new mgadsmproblem(*l.mgadsm)) {}

double laplace::objfun_(const std::vector<double> &x) const
{
	double obj = 0;
	MGA_DSM(x, *mgadsm, obj);
	const size_t sequence_size = (x.size() + 2) / 4;
	double totaltime = 0;
	for (size_t i = 0; i < sequence_size - 1 ; ++i) {
		totaltime += x[i + 4];
	}
	const double delta = totaltime - 8 * 365.25;
	// Penalise trajectory longer than 8 years by 200 meters/s per month.
	obj = std::max<double>(obj,obj + 0.2 / 30 * delta);
	return obj;
}

std::string laplace::solution(const std::vector<double> &x) const
{
	double obj = 0;
	MGA_DSM(x, *mgadsm, obj);
	std::ostringstream s;
	s.precision(15);
	s << std::scientific;
	const size_t seq_size = (x.size() + 2) / 4;
	pagmo_assert((x.size() + 2) % 4 == 0 && seq_size >= 2);
	pagmo_assert(mgadsm->sequence.size() == seq_size);
	s << "Flyby sequence:        ";
	for (size_t i = 0; i < seq_size; ++i) {
		s << mgadsm->sequence[i];
	}
	s << '\n';
	s << "Time of departure:     " << x[0] << '\n';
	s << "Vinf polar components: ";
	for (size_t i = 0; i < 3; ++i) {
		s << x[i + 1] << ' ';
	}
	s << '\n';
	double totaltime = 0;
	for (size_t i = 0; i < seq_size - 1; ++i) {
		s << "Leg time of flight:    " << x[i + 4] << '\n';
		totaltime += x[i + 4];
	}
	s << "Total time of flight:  " << totaltime << '\n';
	for (size_t i = 0; i < seq_size - 2; ++i) {
		s << "Flyby radius:          " << x[i + 2 * (seq_size + 1)] << '\n';
	}
	for (size_t i = 0; i < seq_size - 2; ++i) {
		s << "Vinf at flyby:         " << std::sqrt(mgadsm->vrelin_vec[i]) << '\n';
	}
	return s.str();
}

std::ostream &laplace::print(std::ostream &s) const
{
	base::print(s);
	s << "Flyby sequence: ";
	const size_t seq_size = mgadsm->sequence.size();
	for (size_t i = 0; i < seq_size; ++i) {
		s << mgadsm->sequence[i];
	}
	s << '\n';
	return s;
}

}
}
